library(dplyr)
library(tidyr)
library(stringr)
library(fcexpr)
library(rlang)
library(rstatix)
library(purrr)


# Function to load workspace with caching
load_workspace_cached <- function(wsp_file, rds_file, overwrite = FALSE) {
    # Compute the current MD5 checksum of the FlowJo workspace file.
    # tools::md5sum() returns a named character vector, so strip the name.
    current_md5 <- unname(tools::md5sum(wsp_file))

    # Check whether the cache file already exists.
    rds_exists <- file.exists(rds_file)

    # Default to rebuilding the cache unless a valid cached object is found.
    rebuild_cache <- overwrite || !rds_exists

    if (!rebuild_cache) {
        # Try to read the cached object. If it is malformed or unreadable,
        # fall back to rebuilding the cache.
        cached <- tryCatch(readRDS(rds_file), error = function(e) NULL)

        # A valid cache is expected to be a list with two elements:
        #   $workspace : the object returned by wsx_get_popstats()
        #   $md5       : the MD5 checksum of the WSP file at cache build time
        valid_cache <- is.list(cached) &&
            all(c("workspace", "md5") %in% names(cached)) &&
            is.character(cached$md5) &&
            length(cached$md5) == 1

        if (!valid_cache || !identical(cached$md5, current_md5)) {
            rebuild_cache <- TRUE
        }
    }

    if (rebuild_cache) {
        # Re-parse the FlowJo workspace because the cache is missing,
        # invalid, overwritten, or stale.
        workspace <- wsx_get_popstats(ws = wsp_file)

        # Store both the parsed workspace and the checksum in one RDS file.
        cached <- list(
            workspace = workspace,
            md5 = current_md5
        )

        saveRDS(cached, file = rds_file)
    } else {
        # Reuse the cached workspace object.
        workspace <- cached$workspace
    }

    return(workspace)
}

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

# -------------------------------------------------------------------
# plot_padj()
#
# Purpose:
#   Create a stat.test table compatible with ggpubr::stat_pvalue_manual()
#   for both:
#     1) faceted plots
#     2) non-faceted plots
#
# Features:
#   - pairwise t-tests for user-specified comparisons
#   - optional multiple-testing correction
#   - correction can be done within facet or across all tests
#   - manual y.positions so free_y facets behave properly
#   - adds comparison-level summary statistics
#   - preserves the USER-SUPPLIED comparison direction for:
#       * mean1, mean2
#       * sd1, sd2
#       * median1, median2
#       * cohen
#       * foldmean
#       * foldmedian
#   - explicitly records fold direction via:
#       * denominator
#       * numerator
#
# Important note on comparison direction:
#   rstatix::pairwise_t_test() may reorder the pair labels in its output
#   (for example, based on factor level order). That is fine for the p-value,
#   because the same two groups are being compared either way.
#
#   However, direction-dependent summary statistics ARE NOT symmetric.
#   For example:
#       c("10K", "25K") -> foldmean = mean(25K) / mean(10K)
#       c("25K", "10K") -> foldmean = mean(10K) / mean(25K)
#
#   Therefore this function stores the requested direction separately and uses
#   it when computing effect sizes and fold changes, even if the test output
#   displays group1/group2 in a different order.
#
# Arguments:
#   data            : data frame containing the variables of interest
#   comparisons     : list of pairwise comparisons, e.g.
#                     list(c("5K", "10K"), c("5K", "25K"))
#   p_adjust_method : p-value adjustment method. Examples:
#                     "BH", "holm", "bonferroni", "none"
#   y_col           : name of the numeric response column (string)
#   x_col           : name of the grouping column being compared (string)
#   facet_col       : name of the faceting column (string), or NULL
#                     for non-faceted analyses
#   id_col          : name of the subject/sample identifier column (string).
#                     REQUIRED when paired = TRUE. plot_padj() does NOT infer
#                     pairing from row order or any other implicit means --
#                     rstatix's paired t.test() pairs strictly by POSITION
#                     (the i-th row of group1 with the i-th row of group2),
#                     so plot_padj() uses id_col to both validate and enforce
#                     that positional pairing is correct:
#                       - errors out (printing the offending rows) if any
#                         id/facet does not have EXACTLY one row per compared
#                         group -- a missing sample or a duplicate would
#                         otherwise silently pair the wrong subjects together
#                         even when both groups end up the same length
#                       - sorts the data by id_col (within facet) before
#                         running the paired test/effect size, so the
#                         positional pairing lines up by subject
#   adjust_scope    : "intrafacet" or "total"
#                     - "intrafacet": adjust p-values separately within each facet
#                     - "total"     : adjust p-values across all comparisons
#   paired          : logical; whether the comparisons are paired
#   pool_sd         : logical; whether to pool SDs for the t-test and
#                     unpaired Cohen's d
#   step_fraction   : fraction of the panel maximum used to vertically
#                     separate significance brackets
#   min_step        : minimum absolute vertical spacing between brackets
#   hide_ns         : logical; if TRUE, remove non-significant rows after
#                     p-value adjustment (based on p.adj if available)
#   alpha           : significance threshold used when hide_ns = TRUE
#
# Returns:
#   A data frame suitable for ggpubr::stat_pvalue_manual(), containing
#   p-values, adjusted p-values, significance labels, annotation
#   y-positions, and summary/effect size statistics. Also records the
#   analysis provenance (p_adjust_method, adjust_scope, paired, pool_sd)
#   as constant columns so downstream consumers of the table know how it
#   was computed without needing to reread the call site.
#
# Notes:
#   - mean1 / median1 / sd1 correspond to the FIRST group supplied by the user
#   - mean2 / median2 / sd2 correspond to the SECOND group supplied by the user
#   - denominator = first group supplied by the user
#   - numerator   = second group supplied by the user
#   - foldmean   = mean(user_group2)   / mean(user_group1)
#   - foldmedian = median(user_group2) / median(user_group1)
#   - fold values are set to NA if the denominator is 0 or missing
#   - for paired Cohen's d, the effect size is based on:
#         mean(user_group2 - user_group1) / sd(user_group2 - user_group1)
#   - for unpaired Cohen's d:
#         if pool_sd = TRUE, pooled SD is used
#         if pool_sd = FALSE, the denominator is mean(sd1, sd2)
#   - the "topmost" comparison for y.position spacing is defined as the
#     last row remaining within each facet after filtering/adjustment
# -------------------------------------------------------------------

plot_padj <- function(
    data,
    comparisons,
    p_adjust_method = "BH",
    y_col,
    x_col,
    facet_col = NULL,
    id_col = NULL,
    adjust_scope = "intrafacet",
    paired = FALSE,
    pool_sd = FALSE,
    step_fraction = 0.08,
    min_step = 1,
    hide_ns = FALSE,
    alpha = 0.05
) {

    adjust_scope <- match.arg(adjust_scope, c("intrafacet", "total"))

    required_cols <- c(y_col, x_col)
    if (!is.null(facet_col)) {
        required_cols <- c(required_cols, facet_col)
    }
    if (paired && !is.null(id_col)) {
        required_cols <- c(required_cols, id_col)
    }

    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) {
        stop(paste0("Missing columns: ", paste(missing_cols, collapse = ", ")))
    }

    # ---- Pairing guard (paired = TRUE only) ------------------------------------------------
    # plot_padj() pairs samples POSITIONALLY under the hood (via rstatix's
    # paired t.test()), not by matching any column automatically. That means
    # incorrect row order -- or a missing/duplicate sample for one group --
    # can silently pair the wrong subjects together with no error at all.
    # id_col makes this explicit, checked, and enforced instead of assumed.
    if (paired) {
        if (is.null(id_col)) {
            stop(
                "id_col is required when paired = TRUE: plot_padj() pairs samples ",
                "positionally (the i-th row of group1 with the i-th row of group2) ",
                "and cannot verify or guarantee correct subject-level pairing without ",
                "a subject/sample identifier column."
            )
        }

        compared_groups <- unique(unlist(comparisons))
        pairing_keys <- c(facet_col, id_col)

        pairing_check <- data %>%
            filter(.data[[x_col]] %in% compared_groups) %>%
            dplyr::count(across(all_of(c(pairing_keys, x_col)))) %>%
            tidyr::pivot_wider(names_from = all_of(x_col), values_from = n, values_fill = 0)

        bad_pairing <- pairing_check %>%
            filter(if_any(all_of(intersect(compared_groups, colnames(pairing_check))), ~ . != 1))

        if (nrow(bad_pairing) > 0) {
            print(bad_pairing)
            stop(
                "paired = TRUE requires exactly one row per (", id_col,
                if (!is.null(facet_col)) paste0(", ", facet_col) else "",
                ", ", x_col, ") combination for the groups being compared -- see the ",
                "offending rows printed above (0 = missing sample, >1 = duplicate)."
            )
        }

        # Sort so the positional pairing rstatix/t.test() relies on lines up
        # by subject: within each facet, every id contributes exactly one row
        # per compared group (validated above), so sorting by id_col alone
        # guarantees each group's filtered vector comes out in the same id order.
        data <- data %>% arrange(across(all_of(pairing_keys)))
    }

    # ---- Sample-size guard -------------------------------------------------------------------
    # stats::t.test() needs >= 2 observations per compared group to compute a
    # variance (>= 2 pairs, i.e. subjects, when paired = TRUE). Without this
    # check, a facet that is missing -- or has only one row of -- a compared
    # group makes rstatix::pairwise_t_test() abort deep inside its call stack
    # with an opaque "not enough 'x'/'y' observations" error that does not
    # say which facet/group is at fault. Catch it here instead, with the
    # offending rows printed, before the test is ever run.
    {
        compared_groups <- unique(unlist(comparisons))
        size_keys <- if (!is.null(facet_col)) c(facet_col, x_col) else x_col

        group_sizes <- data %>%
            filter(.data[[x_col]] %in% compared_groups) %>%
            dplyr::count(across(all_of(size_keys)), name = "n")

        if (!is.null(facet_col)) {
            group_sizes <- tidyr::expand_grid(
                    !!facet_col := unique(data[[facet_col]]),
                    !!x_col := compared_groups
                ) %>%
                left_join(group_sizes, by = c(facet_col, x_col)) %>%
                mutate(n = ifelse(is.na(n), 0L, n))
        } else {
            group_sizes <- tibble(!!x_col := compared_groups) %>%
                left_join(group_sizes, by = x_col) %>%
                mutate(n = ifelse(is.na(n), 0L, n))
        }

        insufficient <- group_sizes %>% filter(n < 2)

        if (nrow(insufficient) > 0) {
            print(insufficient)
            stop(
                "Each group being compared needs at least 2 observations",
                if (!is.null(facet_col)) paste0(" per ", facet_col) else "",
                " (0 = group entirely missing) -- see the offending rows printed above."
            )
        }
    }

    test_formula <- as.formula(paste(y_col, "~", x_col))

    bad_length <- purrr::keep(comparisons, ~ length(.x) != 2)
    if (length(bad_length) > 0) {
        stop(
            "Each element of `comparisons` must contain exactly 2 groups; ",
            "found comparison(s) with a different length: ",
            paste(vapply(bad_length, paste, character(1), collapse = ", "), collapse = "; ")
        )
    }

    comparison_tbl <- purrr::map_dfr(seq_along(comparisons), function(i) {
        comp <- comparisons[[i]]
        tibble(
            comparison_id = i,
            requested_group1 = comp[1],
            requested_group2 = comp[2],
            pair_key = paste(sort(comp), collapse = "__")
        )
    })

    dup_pairs <- unique(comparison_tbl$pair_key[duplicated(comparison_tbl$pair_key)])
    if (length(dup_pairs) > 0) {
        stop(
            "`comparisons` contains duplicate/overlapping pairs (order-independent): ",
            paste(dup_pairs, collapse = ", ")
        )
    }

    p_to_signif <- function(pvals) {
        case_when(
            is.na(pvals)    ~ NA_character_,
            pvals <= 1e-04  ~ "****",
            pvals <= 1e-03  ~ "***",
            pvals <= 1e-02  ~ "**",
            pvals <= 5e-02  ~ "*",
            TRUE            ~ "ns"
        )
    }

    compute_cohen <- function(df_pair, g1, g2) {
        if (paired) {
            # Align by id_col (one row per subject, one column per group) so
            # that dropping missing values can't desync the two vectors the
            # way independently na.omit-ing each group's vector would.
            wide <- df_pair %>%
                filter(.data[[x_col]] %in% c(g1, g2)) %>%
                select(all_of(c(id_col, x_col, y_col))) %>%
                tidyr::pivot_wider(names_from = all_of(x_col), values_from = all_of(y_col))

            x1 <- wide[[g1]]
            x2 <- wide[[g2]]
            keep <- stats::complete.cases(x1, x2)
            x1 <- x1[keep]; x2 <- x2[keep]

            if (length(x1) < 2) return(NA_real_)
            d <- x2 - x1
            s <- sd(d)
            if (is.na(s) || s == 0) return(NA_real_)
            return(mean(d) / s)
        }

        x1 <- df_pair %>% filter(.data[[x_col]] == g1) %>% pull(.data[[y_col]]) %>% na.omit()
        x2 <- df_pair %>% filter(.data[[x_col]] == g2) %>% pull(.data[[y_col]]) %>% na.omit()

        if (length(x1) < 2 || length(x2) < 2) return(NA_real_)

        m1 <- mean(x1); m2 <- mean(x2)
        s1 <- sd(x1); s2 <- sd(x2)

        if (is.na(s1) || is.na(s2)) return(NA_real_)

        if (pool_sd) {
            pooled <- sqrt(((length(x1)-1)*s1^2 + (length(x2)-1)*s2^2)/(length(x1)+length(x2)-2))
            if (pooled == 0) return(NA_real_)
            return((m2 - m1)/pooled)
        }

        denom <- mean(c(s1, s2))
        if (denom == 0) return(NA_real_)
        (m2 - m1)/denom
    }

    get_extra_stats <- function(df_subset) {
        purrr::map_dfr(seq_len(nrow(comparison_tbl)), function(i) {
            g1 <- comparison_tbl$requested_group1[i]
            g2 <- comparison_tbl$requested_group2[i]

            df_pair <- df_subset %>% filter(.data[[x_col]] %in% c(g1, g2))

            # Left-joined onto a fixed two-row (g1, g2) skeleton so a group with
            # zero rows in this facet/subset yields NA summary stats instead of
            # a length-0 column that would break the tibble() below.
            summary_df <- tibble(!!x_col := c(g1, g2)) %>%
                left_join(
                    df_pair %>%
                        group_by(.data[[x_col]]) %>%
                        summarise(  mean_value = mean(.data[[y_col]], na.rm=TRUE),
                                    median_value = median(.data[[y_col]], na.rm=TRUE),
                                    sd_value = sd(.data[[y_col]], na.rm=TRUE),
                                    .groups="drop"),
                    by = x_col
                )

            m1 <- summary_df$mean_value[1]
            m2 <- summary_df$mean_value[2]
            med1 <- summary_df$median_value[1]
            med2 <- summary_df$median_value[2]
            s1 <- summary_df$sd_value[1]
            s2 <- summary_df$sd_value[2]

            tibble(
                comparison_id = i,
                requested_group1 = g1,
                requested_group2 = g2,
                pair_key = comparison_tbl$pair_key[i],
                denominator = g1,
                numerator = g2,
                mean1 = m1,
                mean2 = m2,
                sd1 = s1,
                sd2 = s2,
                median1 = med1,
                median2 = med2,
                cohen = compute_cohen(df_pair, g1, g2),
                foldmean = ifelse(is.na(m1) || m1==0, NA, m2/m1),
                foldmedian = ifelse(is.na(med1) || med1==0, NA, med2/med1)
            )
        })
    }

    add_adjusted_p <- function(stat_tbl) {
        if (p_adjust_method == "none") {
            stat_tbl %>% mutate(p.adj=p, p.adj.signif=p_to_signif(p))
        } else if (!is.null(facet_col) && adjust_scope=="intrafacet") {
            stat_tbl %>% group_by(.data[[facet_col]]) %>%
                mutate( p.adj=p.adjust(p, method=p_adjust_method),
                        p.adj.signif=p_to_signif(p.adj)) %>% ungroup()
        } else {
            stat_tbl %>% mutate(p.adj=p.adjust(p, method=p_adjust_method),
                                p.adj.signif=p_to_signif(p.adj))
        }
    }

    if (!is.null(facet_col)) {

        stat.test <- data %>%
            group_by(.data[[facet_col]]) %>%
            pairwise_t_test(formula=test_formula, comparisons=comparisons,
                            paired=paired, pool.sd=pool_sd, p.adjust.method="none") %>%
            ungroup() %>%
            mutate(pair_key=paste(pmin(group1,group2), pmax(group1,group2), sep="__")) %>%
            left_join(comparison_tbl, by="pair_key")

        stat.test <- add_adjusted_p(stat.test)

        extra.df <- data %>%
            group_split(.data[[facet_col]]) %>%
            purrr::map_dfr(~ get_extra_stats(.x) %>% mutate(!!facet_col := .x[[facet_col]][1]))

        facet_limits <- data %>%
            group_by(.data[[facet_col]]) %>%
            summarise(  y.max=max(.data[[y_col]], na.rm=TRUE),
                        y.range=diff(range(.data[[y_col]], na.rm=TRUE)),
                        .groups="drop")

        stat.test <- stat.test %>%
            left_join(extra.df, by=c(facet_col,"comparison_id","requested_group1","requested_group2","pair_key")) %>%
            left_join(facet_limits, by=facet_col)

        if (hide_ns) stat.test <- stat.test %>% filter(p.adj <= alpha)

        stat.test <- stat.test %>%
            group_by(.data[[facet_col]]) %>%
            mutate(
                step_size = max(step_fraction * first(y.range), min_step),
                y.position = first(y.max) + row_number() * step_size
            ) %>%
            ungroup() %>%
            select(-step_size)

    } else {

        stat.test <- data %>%
            pairwise_t_test(formula=test_formula, comparisons=comparisons,
                            paired=paired, pool.sd=pool_sd, p.adjust.method="none") %>%
            mutate(pair_key=paste(pmin(group1,group2), pmax(group1,group2), sep="__")) %>%
            left_join(comparison_tbl, by="pair_key")

        stat.test <- add_adjusted_p(stat.test)
        extra.df <- get_extra_stats(data)

        y.max <- max(data[[y_col]], na.rm=TRUE)
        y.range <- diff(range(data[[y_col]], na.rm=TRUE))

        stat.test <- stat.test %>% left_join(extra.df,
            by=c("comparison_id","requested_group1","requested_group2","pair_key"))

        if (hide_ns) stat.test <- stat.test %>% filter(p.adj <= alpha)

        step_size <- max(step_fraction * y.range, min_step)

        stat.test <- stat.test %>%
            mutate(
                y.position = y.max + row_number() * step_size
            )
    }

    stat.test <- stat.test %>%
        mutate(
            p_adjust_method = p_adjust_method,
            adjust_scope = adjust_scope,
            paired = paired,
            pool_sd = pool_sd
        )

    stat.test %>% select(-any_of(c("requested_group1","requested_group2","pair_key")))
}
