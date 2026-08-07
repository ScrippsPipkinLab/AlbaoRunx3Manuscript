#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hto_merge.py

Standalone API for pulling Multiplexing Capture (CMO / HTO) counts out of a
Cell Ranger 7 `cellranger multi` run and merging them onto an AnnData object
that has already been fully processed.

Nothing here re-runs, re-reads or depends on the preprocessing notebook. The
merge is a pure .obs / .obsm join keyed on barcode + batch, so .X, .raw,
.layers, embeddings, neighbour graphs, cluster labels and DE results are all
left bit-identical.

Expected keys on the processed object
-------------------------------------
    adata.obs["gem_id"]   per-cell Cell Ranger barcode
    adata.obs["batch"]    "gem1" or "gem2"

Both are configurable (barcode_key= / batch_key=). If the barcode column has
picked up extra suffixes over the course of the analysis, canonicalize_barcodes()
strips them back to the Cell Ranger form. If you are not certain which obs
column holds the barcode, run detect_barcode_key() first -- it scores every
candidate column against the real barcode set and reports match rates rather
than guessing.

Typical use
-----------
    import hto_merge as hto

    frames = hto.read_cmo_all({
        "gem1": "source_data/gem1/outs",
        "gem2": "source_data/gem2/outs",
    })

    tag_map = hto.parse_cmo_map("gem2.csv")        # HTO1..HTO9 -> Base/Null/WT/...

    hto.detect_barcode_key(adata, frames)          # optional sanity check
    hto.merge_hto(adata, frames, tag_map=tag_map)  # writes obs + obsm

Verified against scanpy 1.12 / anndata 0.13.
"""

import os
import re
import glob

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as pl


# Feature type string Cell Ranger writes for CMO features
CMO_FEATURE_TYPE = "Multiplexing Capture"

# Prefixes for what gets written onto the object
CMO_PREFIX = "cmo_"   # summary columns and obsm matrices
HTO_PREFIX = "hto_"   # flat per-tag obs columns, plot-ready

# Barcodes look like AAACCTGAGAAACCAT-1
BARCODE_RE = re.compile(r"[ACGTNacgtn]{8,}")


# ===========================================================================
# 1. Extraction
# ===========================================================================

def _cmo_frame(adata_full):
    """
    Pull the Multiplexing Capture features out of an AnnData that was read
    with gex_only=False. Returns a barcodes x tags DataFrame, or None if the
    matrix carries no CMO features.
    """

    if "feature_types" not in adata_full.var.columns:
        return None

    is_cmo = adata_full.var["feature_types"].astype(str) == CMO_FEATURE_TYPE

    if is_cmo.sum() == 0:
        return None

    subset = adata_full[:, is_cmo]
    counts = subset.X
    counts = counts.toarray() if hasattr(counts, "toarray") else np.asarray(counts)

    return pd.DataFrame(counts,
                        index=subset.obs_names.astype(str),
                        columns=subset.var_names.astype(str))
# End of function


def read_cmo_counts(gem_dir, source="auto", verbose=True):
    """
    Read the CMO counts for one GEM.

    gem_dir : path to that GEM's Cell Ranger `outs` directory, e.g.
              "source_data/gem1/outs"
    source  : "per_sample"  read outs/per_sample_outs/*/count/sample_filtered_feature_bc_matrix
              "raw"         read outs/multi/count/raw_feature_bc_matrix.h5
              "auto"        try per_sample, fall back to raw

    Returns a barcodes x tags DataFrame of raw CMO UMI counts. When read from
    the per-sample directories an extra "cellranger_call" column carries the
    sample_id each barcode was demultiplexed into.

    Note on why this is needed at all: `cellranger multi` already quantifies
    the CMO library into these same matrices, tagged
    feature_types == "Multiplexing Capture". Both sc.read_10x_mtx and
    sc.read_10x_h5 default to gex_only=True and drop those rows silently, which
    is why they are absent from the processed object. Nothing needs re-running
    in Cell Ranger.
    """

    if source in ("auto", "per_sample"):

        pattern = os.path.join(gem_dir, "per_sample_outs", "*",
                               "count", "sample_filtered_feature_bc_matrix")
        sample_dirs = sorted(glob.glob(pattern))

        frames = []

        for path in sample_dirs:

            sample_id = path.split(os.sep)[-3]

            full = sc.read_10x_mtx(path, gex_only=False)
            frame = _cmo_frame(full)
            del full  # Cleanup, these are large

            # This point release strips CMO rows from the per-sample matrices
            if frame is None:
                frames = None
                break

            frame["cellranger_call"] = sample_id
            frames.append(frame)
        # End loop

        if frames is not None and len(frames) > 0:
            out = pd.concat(frames, axis=0, join="outer")
            if verbose:
                print(gem_dir + ": " + str(out.shape[0]) +
                      " barcodes from " + str(len(sample_dirs)) +
                      " per-sample matrices")
            return out

        if source == "per_sample":
            raise ValueError("No Multiplexing Capture features under " + pattern)

        if verbose:
            print(gem_dir + ": per-sample matrices carry no CMO features, "
                  "falling back to the raw matrix.")

    # Fallback / explicit raw route
    path = os.path.join(gem_dir, "multi", "count", "raw_feature_bc_matrix.h5")

    full = sc.read_10x_h5(path, gex_only=False)
    full.var_names_make_unique()
    frame = _cmo_frame(full)
    del full  # Cleanup

    if frame is None:
        raise ValueError("No '" + CMO_FEATURE_TYPE + "' features in " + path +
                         ". Confirm the CMO library was in the multi config.")

    if verbose:
        print(gem_dir + ": " + str(frame.shape[0]) + " droplets from the raw matrix")

    return frame
# End of function


def read_cmo_all(gem_dirs, source="auto", verbose=True):
    """
    Read every GEM in one call.

    gem_dirs : {batch_label: outs_path}, where batch_label matches the values
               in adata.obs["batch"], e.g.
               {"gem1": "source_data/gem1/outs", "gem2": "source_data/gem2/outs"}
    """

    frames = {}

    for label in gem_dirs:
        frames[label] = read_cmo_counts(gem_dirs[label], source=source,
                                        verbose=verbose)
    # End loop

    return frames
# End of function


def parse_cmo_map(config_csv):
    """
    Parse the [samples] block of a cellranger multi config CSV and return
    {cmo_id: sample_id}, e.g. {"HTO1": "Base", "HTO2": "Null", ...}.

    Use it to rename the tag columns into the sample namespace the processed
    object already uses, so the HTO rows carry meaningful labels on a figure.
    """

    mapping = {}
    in_samples = False

    with open(config_csv, "r") as handle:
        for line in handle:
            line = line.strip()

            if line.startswith("["):
                in_samples = line.lower().startswith("[samples]")
                continue

            if not in_samples or line == "" or line.lower().startswith("sample_id"):
                continue

            fields = [f.strip() for f in line.split(",")]

            # cmo_ids is pipe-delimited when a sample carries more than one tag
            for cmo_id in fields[1].split("|"):
                mapping[cmo_id.strip()] = fields[0]
            # End loop
        # End loop

    return mapping
# End of function


# ===========================================================================
# 2. Key resolution
# ===========================================================================

def canonicalize_barcodes(values, suffix="-1", mode="auto", verbose=False):
    """
    Reduce whatever the processed object stores to the Cell Ranger barcode
    form. Handles the usual accretions:

        AAACCTGAGAAACCAT-1          -> AAACCTGAGAAACCAT-1
        AAACCTGAGAAACCAT-1-gem1     -> AAACCTGAGAAACCAT-1
        AAACCTGAGAAACCAT-1-0-1      -> AAACCTGAGAAACCAT-1
        AAACCTGAGAAACCAT            -> AAACCTGAGAAACCAT-1
        gem1_AAACCTGAGAAACCAT-1     -> AAACCTGAGAAACCAT-1

    mode : "nucleotide" match the longest ACGTN run, which is what a real
                        10x barcode is
           "token"      keep the first dash-separated token and re-append the
                        suffix; use for anything that is not literal sequence
           "auto"       nucleotide, falling back to token if nothing parses

    Values that cannot be parsed come back as NaN so they show up in the
    coverage report rather than silently failing to match.
    """

    series = pd.Series(values).astype(str)

    def _nucleotide(value):
        hit = BARCODE_RE.search(value)
        return hit.group(0).upper() + suffix if hit is not None else np.nan
    # End of function

    def _token(value):
        token = value.split("_")[-1].split("-")[0]
        return token + suffix if token != "" else np.nan
    # End of function

    if mode in ("auto", "nucleotide"):
        out = series.map(_nucleotide)
        if mode == "nucleotide" or out.notna().any():
            return pd.Series(out.values, index=series.index)
        if verbose:
            print("No nucleotide runs found; falling back to token parsing.")

    out = series.map(_token)

    return pd.Series(out.values, index=series.index)
# End of function


def detect_barcode_key(adata, frames, batch_key="batch", top=10, verbose=True):
    """
    Score every .obs column (plus obs_names) against the real Cell Ranger
    barcode set and report what fraction of cells each one would match.

    Run this once before merging if there is any doubt about which column
    holds the barcode. It returns a sorted DataFrame; a usable key sits at or
    near 1.00 match, and everything else will be far below.
    """

    # Union of barcodes across all GEMs
    barcodes = set()
    for label in frames:
        barcodes.update(canonicalize_barcodes(frames[label].index).dropna())
    # End loop

    candidates = {}

    # obs_names is worth testing too, in case no obs column survived
    candidates["<obs_names>"] = pd.Series(adata.obs_names, index=adata.obs_names)

    for column in adata.obs.columns:
        if column == batch_key:
            continue
        # Only string-like columns can hold a barcode
        if adata.obs[column].dtype.kind in "biufc":
            continue
        candidates[column] = adata.obs[column]
    # End loop

    rows = []

    for name in candidates:
        canon = canonicalize_barcodes(candidates[name])
        n_parsed = int(canon.notna().sum())
        n_match = int(canon.dropna().isin(barcodes).sum())
        rows.append({
            "key": name,
            "parsed_frac": n_parsed / max(adata.n_obs, 1),
            "match_frac": n_match / max(adata.n_obs, 1),
            "n_unique": int(canon.nunique()),
        })
    # End loop

    report = pd.DataFrame(rows).sort_values("match_frac", ascending=False)
    report = report.head(top).reset_index(drop=True)

    if verbose:
        print(report.to_string(index=False))

    return report
# End of function


def check_keys(adata, frames, barcode_key="gem_id", batch_key="batch"):
    """
    Verify the join can proceed. Checks that both keys exist, that the batch
    labels line up with the keys of `frames`, and that (batch, barcode) is
    unique. Returns True / False and prints what it found.
    """

    ok = True

    for key in [barcode_key, batch_key]:
        if key not in adata.obs.columns:
            print("MISSING obs['" + key + "'].")
            ok = False
    # End loop

    if not ok:
        print("Run detect_barcode_key() to find which column holds the barcode.")
        return False

    obs_batches = set(adata.obs[batch_key].astype(str).unique())
    frame_batches = set(str(k) for k in frames)

    print("batch values in obs : " + str(sorted(obs_batches)))
    print("batch keys in frames: " + str(sorted(frame_batches)))

    if not obs_batches.issubset(frame_batches):
        print("Batches with no matching frame: " +
              str(sorted(obs_batches - frame_batches)))
        ok = False

    canon = canonicalize_barcodes(adata.obs[barcode_key])
    print("barcodes parsed: " + str(int(canon.notna().sum())) + " / " + str(adata.n_obs))
    print("example: " + str(adata.obs[barcode_key].iloc[0]) + " -> " + str(canon.iloc[0]))

    pairs = pd.DataFrame({"b": adata.obs[batch_key].astype(str), "c": canon})
    n_dup = int(pairs.duplicated().sum())
    print("duplicate (batch, barcode) pairs: " + str(n_dup))
    if n_dup > 0:
        print("A duplicated key would fan out the join. Resolve before merging.")
        ok = False

    del canon, pairs  # Cleanup

    return ok
# End of function


# ===========================================================================
# 3. Normalization and per-cell summaries
# ===========================================================================

def clr_normalize(counts):
    """
    Centered log-ratio across tags within each cell, the standard treatment
    for hashing data (equivalent to Seurat NormalizeData(method="CLR",
    margin=2)):

        clr(x_i) = log1p( x_i / geometric_mean(x_1 .. x_n) )
    """

    values = counts.to_numpy(dtype=float)

    geo_mean = np.exp(np.log1p(values).sum(axis=1) / values.shape[1])
    geo_mean[geo_mean == 0] = 1.0

    return pd.DataFrame(np.log1p(values / geo_mean[:, None]),
                        index=counts.index, columns=counts.columns)
# End of function


def summarize_cmo(counts):
    """
    Per-cell summaries over the raw counts: total, top and runner-up tag, top
    fraction, and log2(top / second). The last is the useful one -- a low
    margin marks an inter-sample doublet or an ambiguous call that Cell Ranger
    assigned anyway.
    """

    values = counts.to_numpy(dtype=float)
    tags = np.asarray(counts.columns)

    order = np.argsort(-values, axis=1)
    rows = np.arange(values.shape[0])
    top_idx = order[:, 0]
    second_idx = order[:, 1] if values.shape[1] > 1 else order[:, 0]

    top_val = values[rows, top_idx]
    second_val = values[rows, second_idx]
    total = values.sum(axis=1)

    summary = pd.DataFrame(index=counts.index)
    summary[CMO_PREFIX + "total"] = total
    summary[CMO_PREFIX + "top"] = tags[top_idx]
    summary[CMO_PREFIX + "top_count"] = top_val
    summary[CMO_PREFIX + "top_frac"] = np.divide(
        top_val, total, out=np.zeros_like(top_val), where=total > 0)
    summary[CMO_PREFIX + "second"] = tags[second_idx]
    summary[CMO_PREFIX + "second_count"] = second_val
    summary[CMO_PREFIX + "margin"] = np.log2((top_val + 1.0) / (second_val + 1.0))

    return summary
# End of function


# ===========================================================================
# 4. Merge
# ===========================================================================

def merge_hto(adata, frames, tag_map=None, barcode_key="gem_id",
              batch_key="batch", add_obs_columns=True, group_key=None,
              copy=False):
    """
    Merge the CMO counts onto a processed object. Modifies in place (or
    returns a copy) and writes:

        adata.obsm["cmo_counts"]   cells x tags, raw UMI counts
        adata.obsm["cmo_clr"]      cells x tags, CLR normalized
        adata.obs["cmo_total"]     total CMO UMIs on the cell
        adata.obs["cmo_top"]       argmax tag
        adata.obs["cmo_top_frac"]  argmax / total
        adata.obs["cmo_second"]    runner-up tag
        adata.obs["cmo_margin"]    log2(top / second)
        adata.obs["cmo_cellranger_call"]   Cell Ranger's own assignment
        adata.obs["hto_<tag>"]     one flat column per tag (add_obs_columns)

    tag_map   : optional {cmo_id: sample_id} from parse_cmo_map()
    group_key : optional obs column holding the sample of origin; if given, a
                cmo_agree column is added comparing it to cmo_top

    Nothing on the feature axis is touched. Note that .obs columns are exactly
    what sc.pl.dotplot needs -- it resolves var_names through sc.get.obs_df,
    which accepts obs column names alongside genes, and reads them regardless
    of use_raw. Appending tags to .var would instead break plotting on an
    object whose .raw was frozen before the tags existed.
    """

    if copy:
        adata = adata.copy()

    for key in [barcode_key, batch_key]:
        if key not in adata.obs.columns:
            raise KeyError("adata.obs is missing '" + key +
                           "'. Run detect_barcode_key() to locate it.")
    # End loop

    # ---- Build one long frame keyed on (batch, barcode) --------------------
    stacked = []

    for label in frames:

        frame = frames[label].copy()
        frame.insert(0, "_barcode", canonicalize_barcodes(frame.index).values)
        frame.insert(0, "_batch", str(label))
        stacked.append(frame)
    # End loop

    stacked = pd.concat(stacked, axis=0, join="outer", ignore_index=True)
    stacked = stacked.dropna(subset=["_barcode"])

    duplicated = stacked.duplicated(subset=["_batch", "_barcode"])
    if duplicated.any():
        raise ValueError(str(int(duplicated.sum())) + " barcode(s) appear more "
                         "than once within a batch of the Cell Ranger output.")
    del duplicated  # Cleanup

    meta_cols = ["_batch", "_barcode", "cellranger_call"]
    tag_cols = [c for c in stacked.columns if c not in meta_cols]

    # ---- Align to the processed object, in its own cell order -------------
    keys = pd.MultiIndex.from_arrays([
        adata.obs[batch_key].astype(str).values,
        canonicalize_barcodes(adata.obs[barcode_key]).values,
    ])

    aligned = stacked.set_index(["_batch", "_barcode"]).reindex(keys)
    aligned.index = adata.obs_names

    n_missing = int(aligned[tag_cols].isna().all(axis=1).sum())
    print("Matched CMO counts for " + str(adata.n_obs - n_missing) + " / " +
          str(adata.n_obs) + " cells.")
    if n_missing == adata.n_obs:
        raise ValueError("No cell matched a CMO record. The barcode or batch "
                         "key is wrong -- run detect_barcode_key(adata, frames) "
                         "and check_keys(adata, frames) before merging.")

    if n_missing > 0:
        print("Warning: " + str(n_missing) + " cell(s) had no CMO record and "
              "are recorded as 0. If that fraction is large, the barcode or "
              "batch key is probably wrong -- check detect_barcode_key().")

    counts = aligned[tag_cols].fillna(0.0)

    # HTO1..HTO9 -> Base/Null/WT/... so the rows carry readable labels
    if tag_map is not None:
        counts = counts.rename(columns=tag_map)

    counts = counts[sorted(counts.columns)]

    # ---- Write ------------------------------------------------------------
    adata.obsm[CMO_PREFIX + "counts"] = counts
    adata.obsm[CMO_PREFIX + "clr"] = clr_normalize(counts)

    summary = summarize_cmo(counts)
    for column in summary.columns:
        adata.obs[column] = summary[column].values
    # End loop

    adata.obs[CMO_PREFIX + "top"] = pd.Categorical(adata.obs[CMO_PREFIX + "top"])
    adata.obs[CMO_PREFIX + "second"] = pd.Categorical(adata.obs[CMO_PREFIX + "second"])

    if add_obs_columns:
        for tag in counts.columns:
            adata.obs[HTO_PREFIX + str(tag)] = counts[tag].values
        # End loop

    if "cellranger_call" in aligned.columns:
        adata.obs[CMO_PREFIX + "cellranger_call"] = pd.Categorical(
            aligned["cellranger_call"].astype(object))

    # Optional cross-check against a sample-of-origin column
    if group_key is not None and group_key in adata.obs.columns:

        shared = set(counts.columns).intersection(
            set(adata.obs[group_key].astype(str).unique()))

        if len(shared) == 0:
            print("Skipping the agreement check: tag names and '" + group_key +
                  "' labels are in different namespaces. Pass tag_map=.")
        else:
            adata.obs[CMO_PREFIX + "agree"] = pd.Categorical(np.where(
                adata.obs[CMO_PREFIX + "top"].astype(str).values ==
                adata.obs[group_key].astype(str).values, "agree", "disagree"))
            print(adata.obs[CMO_PREFIX + "agree"].value_counts())
        del shared  # Cleanup

    del stacked, aligned, summary  # Cleanup

    return adata
# End of function


def hto_obs_keys(adata, restrict_to=None, order=None):
    """
    List the flat hto_<tag> columns written by merge_hto, optionally filtered
    and ordered for plotting.

    restrict_to : iterable of tag names to keep. Worth using on a subset
                  object -- a tag whose sample is absent carries ambient
                  signal only, and standard_scale="var" will stretch that row
                  across the full colour range as if it were real.
    order       : explicit tag order for the rows
    """

    keys = [c for c in adata.obs.columns if c.startswith(HTO_PREFIX)]
    tags = [k[len(HTO_PREFIX):] for k in keys]

    if restrict_to is not None:
        restrict_to = set(str(t) for t in restrict_to)
        dropped = [t for t in tags if t not in restrict_to]
        tags = [t for t in tags if t in restrict_to]
        if len(dropped) > 0:
            print("Excluding tags: " + str(dropped))
        del dropped  # Cleanup

    if order is not None:
        head = [t for t in order if t in tags]
        tags = head + [t for t in tags if t not in head]
        del head  # Cleanup

    return [HTO_PREFIX + str(t) for t in tags]
# End of function


def hto_composition(adata, groupby, call_key=CMO_PREFIX + "top"):
    """
    Percent composition of each cluster by assigned tag. Useful as the
    supplementary table behind the figure.
    """

    return (pd.crosstab(adata.obs[groupby], adata.obs[call_key].astype(str),
                        normalize="index") * 100.0).round(1)
# End of function


# ===========================================================================
# 5. Plotting
# ===========================================================================

def hto_dotplot(adata, groupby, hto_keys=None, categories_order=None,
                standard_scale="var", call_key=CMO_PREFIX + "top",
                swap_axes=True, figsize=None, ax=None):
    """
    HTO panel with dot semantics stated explicitly:

        colour = mean CLR signal for that tag in that cluster
        size   = percent of cells in that cluster assigned to that tag

    The second definition matters. scanpy's default dot size is the fraction
    of cells above expression_cutoff, and with hashing data ambient CMO leaves
    a nonzero count on nearly every barcode -- so on the default cutoff of 0
    every dot renders at ~100%. Composition is both more informative and
    easier to defend in a response letter.

    Built with DotPlot's dot_color_df / dot_size_df, so neither matrix is
    derived from .X.
    """

    if hto_keys is None:
        hto_keys = hto_obs_keys(adata)

    tags = [k[len(HTO_PREFIX):] for k in hto_keys]
    groups = adata.obs[groupby].astype(str)

    color_df = adata.obs[hto_keys].groupby(groups.values, observed=True).mean()

    size_df = pd.crosstab(groups, adata.obs[call_key].astype(str),
                          normalize="index")
    size_df = size_df.reindex(columns=tags, fill_value=0.0)
    size_df.columns = hto_keys

    color_df = color_df.reindex(size_df.index)
    for frame in (color_df, size_df):
        frame.index.name = None
        frame.columns.name = None
    # End loop

    if categories_order is not None:
        color_df = color_df.reindex(categories_order)
        size_df = size_df.reindex(categories_order)

    dp = sc.pl.DotPlot(adata, var_names=hto_keys, groupby=groupby,
                       categories_order=categories_order,
                       dot_color_df=color_df, dot_size_df=size_df,
                       standard_scale=standard_scale, figsize=figsize, ax=ax)

    if swap_axes:
        dp = dp.swap_axes()

    dp.legend(colorbar_title="mean HTO CLR\n(scaled per tag)",
              size_title="cells assigned\nto tag (%)")

    return dp
# End of function


def stacked_figure(adata, gene_names, hto_keys, groupby, categories_order=None,
                   standard_scale="var", height_ratios=(6, 2), figsize=(6, 10)):
    """
    Existing gene dotplot on top, HTO block beneath it sharing the same
    cluster columns, each panel keeping its own dot-size legend.

    Two panels rather than one because the blocks answer different questions
    -- transcript abundance versus sample provenance -- and a single dotplot
    forces one expression_cutoff and one size legend to describe both. The
    cutoff is global in scanpy; raising it to make the HTO dots meaningful
    would silently change the gene rows of an already-published panel.
    """

    fig = pl.figure(figsize=figsize)
    grid = fig.add_gridspec(2, 1, height_ratios=list(height_ratios), hspace=0.05)

    ax_genes = fig.add_subplot(grid[0])
    ax_hto = fig.add_subplot(grid[1])

    dp_genes = sc.pl.DotPlot(adata, var_names=list(gene_names), groupby=groupby,
                             categories_order=categories_order,
                             standard_scale=standard_scale, ax=ax_genes)
    dp_genes.swap_axes().make_figure()

    dp_hto = hto_dotplot(adata, groupby=groupby, hto_keys=hto_keys,
                         categories_order=categories_order,
                         standard_scale=standard_scale, ax=ax_hto)
    dp_hto.make_figure()

    return fig
# End of function


def single_panel_dotplot(adata, gene_names, hto_keys, groupby,
                         categories_order=None, standard_scale="var",
                         swap_axes=True, figsize=None, **kwds):
    """
    Alternative: one dotplot with the HTO rows bracketed under the genes.
    expression_cutoff stays at 0 so the published gene rows are untouched,
    which means the HTO dot sizes are uninformative -- colour still is.
    """

    var_names = list(gene_names) + list(hto_keys)

    return sc.pl.dotplot(
        adata, var_names, groupby=groupby, categories_order=categories_order,
        standard_scale=standard_scale, expression_cutoff=0.0,
        var_group_positions=[(len(gene_names), len(var_names) - 1)],
        var_group_labels=["HTO"], swap_axes=swap_axes, figsize=figsize,
        show=False, return_fig=True, **kwds)
# End of function


def plot_hto_distribution(adata, hto_keys, bins=60):
    """
    Per-tag histograms, for eyeballing the ambient versus carrying modes.
    """

    n = len(hto_keys)
    fig, axes = pl.subplots(n, 1, figsize=(4, 1.3 * n), sharex=True)
    axes = np.atleast_1d(axes)

    for index in range(0, n):
        axes[index].hist(adata.obs[hto_keys[index]].values, bins=bins)
        axes[index].set_ylabel(hto_keys[index], rotation=0, ha="right", fontsize=8)
    # End loop

    del index  # Cleanup
    fig.tight_layout()

    return fig
# End of function


# ===========================================================================
# Example
# ===========================================================================

if __name__ == "__main__":

    # ---- paths, edit to match wherever this is being run -------------------
    GEM_DIRS = {
        "gem1": "source_data/gem1/outs",
        "gem2": "source_data/gem2/outs",
    }
    CONFIG_CSV = "gem2.csv"
    PROCESSED = "h5ad/<processed>.h5ad"
    CLUSTER_KEY = "leiden"

    # ---- 1. extract --------------------------------------------------------
    frames = read_cmo_all(GEM_DIRS)
    tag_map = parse_cmo_map(CONFIG_CSV)
    print(tag_map)

    # ---- 2. merge ----------------------------------------------------------
    adata = sc.read_h5ad(PROCESSED)

    detect_barcode_key(adata, frames)          # confirm gem_id is the barcode
    check_keys(adata, frames)                  # confirm the join is sound

    merge_hto(adata, frames, tag_map=tag_map,
              barcode_key="gem_id", batch_key="batch")

    # ---- 3. figure ---------------------------------------------------------
    # Keep only tags whose samples are actually in this object
    keys = hto_obs_keys(adata, restrict_to=list(tag_map.values()))

    genes = ["Tcf7", "Slamf6", "Tox", "Bach2", "Jund", "Eomes", "Bcl6", "Cxcr3",
             "Cd28", "Cd27", "Il7r", "Pdcd1", "Ctla4", "Cx3cr1", "Klrg1",
             "Runx1", "Runx2", "Runx3", "Ets1", "Tbx21", "Zeb2", "Klf2",
             "Klf3", "Ifng", "Cxcr6"]
    genes = [g for g in genes if g in adata.var_names or
             (adata.raw is not None and g in adata.raw.var_names)]

    fig = stacked_figure(adata, genes, keys, groupby=CLUSTER_KEY)
    fig.savefig("figures/bubble_with_HTO.pdf", bbox_inches="tight")

    print(hto_composition(adata, CLUSTER_KEY))

    # ---- 4. optionally persist --------------------------------------------
    adata.write_h5ad(filename=PROCESSED.replace(".h5ad", "_HTO.h5ad"),
                     compression="gzip", compression_opts=9)

# End of script