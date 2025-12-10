library(dplyr)
library(tidyr)
library(stringr)

# Define a function that adds percentages of two populations, Path1 and Path2
# Then combine subpopulations within these two populations by weighted averaging
combinePercents <- function(df, resultingPath, Path1, Path2, constantPath,
                            subsets = c("TE", "DP", "EE", "MP", "Tcm", "Tpm", "Tem")){
    # First add percents to determine total
    resultingPath <- paste0(constantPath, resultingPath)
    Path1 <- paste0(constantPath, Path1)
    Path2 <- paste0(constantPath, Path2)
    df[resultingPath] <- df[Path1] + df[Path2]

    # Do weighted averaging of subsets
    for(sub in subsets){
        subPath1 <- paste0(Path1, "/", sub)
        subPath2 <- paste0(Path2, "/", sub)
        subResultingPath <- paste0(resultingPath, "/", sub)
        df[subResultingPath] <- ( df[subPath1] * df[Path1] + df[subPath2] * df[Path2] ) / ( df[Path1] + df[Path2] )
    }

    return(df)
}


#' Compute percent-of-total from FlowJo hierarchical gating
#'
#' @param df A data frame with columns:
#'   - ID: unique sample identifier
#'   - PopulationFullPath: hierarchical population name (e.g. "Transferred/Transduced/DP")
#'   - Percent: percent-of-parent from FlowJo (0–100)
#'   Other columns are preserved.
#' @param id_col Name of the ID column (default: "ID")
#' @param path_col Name of the population path column (default: "PopulationFullPath")
#' @param percent_col Name of the percent-of-parent column (default: "Percent")
#'
#' @return A tibble with a new column `PercentOfTotal`
compute_percent_of_total <- function(   df,
                                        id_col = "ID",
                                        path_col = "PopulationFullPath",
                                        percent_col = "Percent") {

    # Standardize column symbols
    id_col <- rlang::sym(id_col)
    path_col <- rlang::sym(path_col)
    percent_col <- rlang::sym(percent_col)

    df %>%
        group_by(!!id_col) %>%
        mutate(
        Depth = str_count(!!path_col, "/") + 1,
        ParentPath = if_else(
            str_detect(!!path_col, "/"),
            sub("/[^/]+$", "", !!path_col),
            NA_character_
        )
        ) %>%
        arrange(Depth, .by_group = TRUE) %>%
        group_modify(~ {
        d <- .x
        d$PercentOfTotal <- NA_real_

        # Depth 1 rows: already percent-of-total
        d$PercentOfTotal[d$Depth == 1] <- d[[rlang::as_name(percent_col)]][d$Depth == 1]

        # Propagate down hierarchy
        if (max(d$Depth, na.rm = TRUE) >= 2) {
            for (lev in sort(unique(d$Depth[d$Depth >= 2]))) {
            idx <- which(d$Depth == lev)
            if (length(idx)) {
                parent_pot <- d$PercentOfTotal[match(d$ParentPath[idx], d[[rlang::as_name(path_col)]])]
                d$PercentOfTotal[idx] <- (d[[rlang::as_name(percent_col)]][idx] / 100) * parent_pot
            }
            }
        }
        d
        }) %>%
        ungroup()
}

#' Derive population columns from a hierarchical FlowJo path
#'
#' @param df A data.frame/tibble containing a hierarchical path column.
#' @param path_col Name of the column holding the hierarchical path
#'        (default: "PopulationFullPath"). Paths use "/" as the delimiter.
#' @param max_depth Optional integer specifying the maximum depth to expand.
#'        If NULL (default), it's inferred as the maximum number of segments
#'        across all rows in `path_col`.
#'
#' @return The input `df` with additional columns:
#'   - Population
#'   - Population1Deriv, Population2Deriv, ..., Population{max_depth-1}Deriv
#'   Missing levels are filled with NA.
add_population_derivatives <- function( df,
                                        path_col = "PopulationFullPath",
                                        max_depth = NULL) {
    # Resolve the column programmatically
    path_col_sym <- rlang::sym(path_col)

    # Extract and clean path strings
    paths <- df[[path_col]]
    if (!is.character(paths)) {
        stop("`", path_col, "` must be a character column.")
    }

    # Split each path into its components by "/"
    # - Trim whitespace
    # - Remove accidental empty segments if any (e.g., leading/trailing '/')
    parts_list <- str_split(str_trim(paths), "/", simplify = FALSE) |>
        lapply(function(x) x[nzchar(x)])

    # Determine maximum depth if not provided
    depths <- vapply(parts_list, length, integer(1))
    if (is.null(max_depth)) {
        max_depth <- if (length(depths)) max(depths, na.rm = TRUE) else 1L
    }
    if (max_depth < 1L) max_depth <- 1L

    # Build output column names:
    # Level 1 -> "Population"
    # Level 2..max_depth -> "Population{level-1}Deriv"
    deriv_names <- if (max_depth >= 2) {
        paste0("Population", seq_len(max_depth - 1), "Deriv")
    } else character(0)
    out_names <- c("Population", deriv_names)

    # Row-wise pad each vector to `max_depth` with NA, then rbind into a matrix
    pad_to_depth <- function(v) {
        x <- rep(NA_character_, max_depth)
        x[seq_along(v)] <- v
        x
    }
    parts_mat <- do.call(rbind, lapply(parts_list, pad_to_depth))
    colnames(parts_mat) <- out_names

    # Bind derived columns to original data
    # (If you already have columns with the same names, they’ll be overwritten.)
    bind_cols(df, as_tibble(parts_mat))
}


#' Finalize FlowJo-derived table for arbitrary hierarchy depth
#'
#' - Removes helper columns (Depth, ParentPath, PopulationFullPath) if present
#' - Detects and orders Population derivatives: Population, Population1Deriv..PopulationNDeriv
#' - Relocates ID columns first, then Population columns, then everything else,
#'   ending with Percent/PercentOfTotal
#' - Sorts rows within ID by hierarchical path (Population -> Population1Deriv -> ... -> PopulationNDeriv)
#' - Optionally rounds Percent and PercentOfTotal
#'
#' @param df data.frame/tibble after compute_percent_of_total() and add_population_derivatives()
#' @param id_cols character vector of identifier columns to place first
#' @param percent_cols character vector of percent columns to place last
#' @param drop_cols helper columns to remove if present
#' @param round_digits integer digits to round Percent columns (set NULL to skip rounding)
#' @return tibble
finalize_flow_table <- function(df,
                                id_cols      = c("ID", "Group", "Replicate"),
                                percent_cols = c("Percent", "PercentOfTotal"),
                                drop_cols    = c("Depth", "ParentPath", "PopulationFullPath"),
                                round_digits = 6) {

    # ---- 1) Drop helper columns if they exist ----------------------------------------------
    out <- df %>% select(-any_of(drop_cols))

    # ---- 2) Identify "Population" + all "Population{n}Deriv" columns dynamically -----------
    nm <- names(out)

    has_root <- "Population" %in% nm
    deriv_idx <- str_which(nm, "^Population\\d+Deriv$")
    deriv_nums <- as.integer(str_match(nm[deriv_idx], "^Population(\\d+)Deriv$")[, 2])
    deriv_order <- deriv_idx[order(deriv_nums)]
    pop_cols <- c(if (has_root) "Population" else character(0), nm[deriv_order])

    # ---- 3) Ensure Percent columns are numeric ---------------------------------------------
    make_numeric <- intersect(percent_cols, nm)
    if (length(make_numeric)) {
        out <- out %>%
        mutate(across(all_of(make_numeric), ~ suppressWarnings(as.numeric(.x))))
    }

    # ---- 4) Create a hierarchy sort key ----------------------------------------------------
    if (length(pop_cols)) {
        out <- out %>%
        mutate(
            .sort_key = do.call(
            paste, c(lapply(out[pop_cols], function(x) ifelse(is.na(x), "", x)), sep = "\t")
            )
        )
    } else {
        out <- out %>% mutate(.sort_key = "")
    }

    # ---- 5) Relocate columns ---------------------------------------------------------------
    id_cols_keep      <- intersect(id_cols, names(out))
    percent_cols_keep <- intersect(percent_cols, names(out))

    out <- out %>%
        relocate(all_of(id_cols_keep), .before = 1) %>%
        relocate(all_of(pop_cols), .after = last_col()) %>%
        relocate(all_of(pop_cols), .after = all_of(id_cols_keep)) %>%
        relocate(all_of(percent_cols_keep), .after = last_col())

    # ---- 6) Sort rows ----------------------------------------------------------------------
    if (length(id_cols_keep)) {
        out <- out %>% arrange(across(all_of(id_cols_keep)), .sort_key)
    } else {
        out <- out %>% arrange(.sort_key)
    }

    # ---- 7) Optional rounding --------------------------------------------------------------
    if (!is.null(round_digits) && is.finite(round_digits) && round_digits >= 0) {
        out <- out %>%
        mutate(across(all_of(percent_cols_keep), ~ round(.x, round_digits)))
    }

    # ---- 8) Drop sort key ------------------------------------------------------------------
    out %>% select(-.sort_key)
}