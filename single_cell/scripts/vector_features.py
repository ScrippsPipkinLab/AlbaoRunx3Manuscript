#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vector_features.py

Extract counts for custom reference features (the retroviral construct --
MycRunx3 and anything else added to the reference for the re-alignment) from
a re-run `cellranger multi` output, and append them to the FEATURE axis of an
already-processed AnnData.

Companion to hto_merge.py, same join model: barcode + batch, two GEMs.

Why var is the right axis this time
-----------------------------------
Unlike the CMO tags, MycRunx3 is a genuine transcript quantified in the Gene
Expression library. It belongs on the var axis, and it needs to be there for
sc.pl.dotplot / sc.pl.umap / sc.tl.rank_genes_groups to treat it like any
other gene.

The complication is where exactly to put it on a processed object:

    .raw        log-normalized, full gene set, frozen before scaling.
                This is what sc.pl.* reads by default (use_raw=True), so a
                feature that is not here will not plot. PRIMARY TARGET.
    .X          scaled z-scores, possibly HVG-subset. Appending raw counts
                here mixes units within one matrix. Supported, but opt-in,
                and varm entries have to be padded to keep shapes valid.
    .obs        a plain per-cell column. Not the feature axis at all, but the
                lowest-risk option and enough for a dotplot row.

append_to_raw() is the default because it makes the feature behave like a
gene everywhere downstream without touching .X.

Read the note on multimapping at the bottom of this file before drawing
conclusions about endogenous versus vector-derived signal.

Verified against scanpy 1.12 / anndata 0.13.
"""

import os
import glob

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad

from hto_merge import canonicalize_barcodes, detect_barcode_key


GEX_FEATURE_TYPE = "Gene Expression"


# ===========================================================================
# 1. Finding the feature
# ===========================================================================

def list_features(gem_dir, pattern=None, source="raw"):
    """
    List features in the re-aligned reference, optionally filtered by a
    case-insensitive substring. Use this first to get the exact spelling of
    the construct as it appears in the new reference -- "MycRunx3",
    "Myc-Runx3", "MSCV-MycRunx3" and so on are all plausible, and the join
    silently returns nothing if the name is off by a character.

        list_features("source_data_realign/gem1/outs", "runx")
    """

    path = _matrix_path(gem_dir, source)
    full = _read_matrix(path)

    frame = full.var.copy()
    frame["feature_name"] = full.var_names.astype(str)
    del full  # Cleanup

    if pattern is not None:
        keep = frame["feature_name"].str.contains(pattern, case=False, regex=False)
        frame = frame[keep]
        del keep  # Cleanup

    cols = [c for c in ["feature_name", "gene_ids", "feature_types"]
            if c in frame.columns]

    return frame[cols].reset_index(drop=True)
# End of function


def diff_features(gem_dir, adata, source="raw"):
    """
    Features present in the re-aligned reference but absent from the
    processed object's gene list. On a re-run whose only change was adding the
    viral genome, this is exactly the set of construct features.
    """

    path = _matrix_path(gem_dir, source)
    full = _read_matrix(path)
    new_names = pd.Index(full.var_names.astype(str))
    del full  # Cleanup

    old_names = pd.Index(adata.var_names.astype(str))
    if adata.raw is not None:
        old_names = old_names.union(pd.Index(adata.raw.var_names.astype(str)))

    return list(new_names.difference(old_names))
# End of function


def _matrix_path(gem_dir, source):
    """Resolve which matrix to read for one GEM."""

    if source == "raw":
        return os.path.join(gem_dir, "multi", "count", "raw_feature_bc_matrix.h5")

    if source == "per_sample":
        return sorted(glob.glob(os.path.join(
            gem_dir, "per_sample_outs", "*", "count",
            "sample_filtered_feature_bc_matrix")))

    raise ValueError("source must be 'raw' or 'per_sample'.")
# End of function


def _read_matrix(path):
    """Read either an .h5 or an mtx directory, keeping all feature types."""

    if isinstance(path, list):
        raise ValueError("_read_matrix takes a single path.")

    if str(path).endswith(".h5"):
        full = sc.read_10x_h5(path, gex_only=False)
    else:
        full = sc.read_10x_mtx(path, gex_only=False)

    full.var_names_make_unique()

    return full
# End of function


# ===========================================================================
# 2. Extraction
# ===========================================================================

def read_feature_counts(gem_dir, features, source="raw", with_totals=True,
                        verbose=True):
    """
    Pull raw counts for named features from one GEM of the re-aligned run.

    gem_dir     : that GEM's `outs` directory from the RE-RUN
    features    : list of feature names, e.g. ["MycRunx3"]
    source      : "raw"        outs/multi/count/raw_feature_bc_matrix.h5
                  "per_sample" outs/per_sample_outs/*/count/sample_filtered...
    with_totals : also return a "_total_gex" column, the summed Gene
                  Expression counts per barcode in the new reference. Needed
                  if the values are going to be log-normalized to match .raw.

    "raw" is the default deliberately. Cell calling can shift slightly between
    Cell Ranger runs, so the re-run's filtered matrices may not contain every
    barcode the processed object kept. The raw matrix contains all of them,
    and the merge subsets to the cells that already exist.

    Returns a barcodes x features DataFrame of raw UMI counts.
    """

    if isinstance(features, str):
        features = [features]

    paths = _matrix_path(gem_dir, source)
    paths = paths if isinstance(paths, list) else [paths]

    frames = []

    for path in paths:

        full = _read_matrix(path)

        missing = [f for f in features if f not in full.var_names]
        if len(missing) > 0:
            raise KeyError("Not in the re-aligned reference: " + str(missing) +
                           ". Run list_features(gem_dir, '<substring>') to get "
                           "the exact feature name.")

        sub = full[:, features]
        counts = sub.X
        counts = counts.toarray() if hasattr(counts, "toarray") else np.asarray(counts)

        frame = pd.DataFrame(counts, index=sub.obs_names.astype(str),
                             columns=list(features))

        # Library size in the new reference, for optional normalization
        if with_totals:
            is_gex = full.var["feature_types"].astype(str) == GEX_FEATURE_TYPE
            totals = full[:, is_gex].X.sum(axis=1)
            frame["_total_gex"] = np.asarray(totals).ravel()
            del is_gex, totals  # Cleanup

        del full, sub  # Cleanup, these are large
        frames.append(frame)
    # End loop

    out = pd.concat(frames, axis=0, join="outer") if len(frames) > 1 else frames[0]
    del frames  # Cleanup

    if verbose:
        print(gem_dir + ": " + str(out.shape[0]) + " barcodes, features " +
              str(list(features)))
        print("  cells with nonzero signal: " +
              str({f: int((out[f] > 0).sum()) for f in features}))

    return out
# End of function


def read_feature_counts_all(gem_dirs, features, source="raw", with_totals=True,
                            verbose=True):
    """
    gem_dirs : {batch_label: outs_path}, batch_label matching adata.obs["batch"]
    """

    frames = {}

    for label in gem_dirs:
        frames[label] = read_feature_counts(gem_dirs[label], features,
                                            source=source,
                                            with_totals=with_totals,
                                            verbose=verbose)
    # End loop

    return frames
# End of function


# ===========================================================================
# 3. Alignment to the processed object
# ===========================================================================

def align_to_adata(adata, frames, features, barcode_key="gem_id",
                   batch_key="batch", verbose=True):
    """
    Reindex the per-GEM frames onto adata's cells, in adata's order. Returns a
    DataFrame of raw counts indexed by adata.obs_names, plus the matched
    "_total_gex" column when present.

    This is the same barcode + batch join used for the HTOs, so anything that
    worked there works here.
    """

    if isinstance(features, str):
        features = [features]

    for key in [barcode_key, batch_key]:
        if key not in adata.obs.columns:
            raise KeyError("adata.obs is missing '" + key + "'.")
    # End loop

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
        raise ValueError(str(int(duplicated.sum())) + " barcode(s) duplicated "
                         "within a batch of the re-aligned output.")
    del duplicated  # Cleanup

    keys = pd.MultiIndex.from_arrays([
        adata.obs[batch_key].astype(str).values,
        canonicalize_barcodes(adata.obs[barcode_key]).values,
    ])

    aligned = stacked.set_index(["_batch", "_barcode"]).reindex(keys)
    aligned.index = adata.obs_names

    n_missing = int(aligned[features].isna().all(axis=1).sum())

    if n_missing == adata.n_obs:
        raise ValueError("No cell matched the re-aligned run. Check the "
                         "barcode / batch keys with detect_barcode_key().")

    if verbose:
        print("Matched " + str(adata.n_obs - n_missing) + " / " +
              str(adata.n_obs) + " cells to the re-aligned run.")
        if n_missing > 0:
            print("Warning: " + str(n_missing) + " unmatched cell(s) recorded "
                  "as 0. A large number here means the two runs disagree on "
                  "barcodes, not that the construct is absent.")

    del stacked  # Cleanup

    return aligned.fillna(0.0)
# End of function


def _scale_values(aligned, features, adata, values, target_sum,
                  total_counts_key, verbose=True):
    """
    Convert raw counts to whatever scale the destination matrix is on.

    values : "counts"  raw UMI counts, unchanged
             "lognorm" log1p(count * target_sum / library_size), matching the
                       normalize_total + log1p that .raw holds
             "cpm"     counts per million, no log

    Library size is taken from adata.obs[total_counts_key] if that column
    exists, otherwise from the "_total_gex" column carried over from the
    re-aligned run. The two differ slightly -- the re-run's totals include the
    construct itself and any other added reference sequence.
    """

    raw = aligned[features].astype(float)

    if values == "counts":
        return raw

    if total_counts_key is not None and total_counts_key in adata.obs.columns:
        library = adata.obs[total_counts_key].astype(float).values
        if verbose:
            print("Library size from adata.obs['" + total_counts_key + "'].")
    elif "_total_gex" in aligned.columns:
        library = aligned["_total_gex"].astype(float).values
        if verbose:
            print("Library size from the re-aligned run's Gene Expression totals.")
    else:
        raise ValueError("No library size available. Pass total_counts_key= or "
                         "read the counts with with_totals=True.")

    library = np.where(library > 0, library, 1.0)

    scaled = raw.to_numpy() * (float(target_sum) / library[:, None])

    if values == "lognorm":
        scaled = np.log1p(scaled)
    elif values != "cpm":
        raise ValueError("values must be 'counts', 'lognorm' or 'cpm'.")

    return pd.DataFrame(scaled, index=raw.index, columns=raw.columns)
# End of function


# ===========================================================================
# 4. Writing onto the object
# ===========================================================================

def append_to_raw(adata, frames, features, values="lognorm", target_sum=2e4,
                  total_counts_key="total_counts", barcode_key="gem_id",
                  batch_key="batch", also_obs=True, copy=False, verbose=True):
    """
    Append the construct feature(s) to .raw, which is what sc.pl.dotplot,
    sc.pl.umap and friends read by default. RECOMMENDED ROUTE.

    .raw was frozen after normalize_total + log1p and before scaling, so
    values="lognorm" puts the new feature on the same scale as every other
    gene in that matrix. values="counts" would drop raw UMI counts into a
    log-normalized matrix -- available, but the colour scale on any plot mixing
    it with real genes then means two different things.

    also_obs writes a plain adata.obs["<feature>_counts"] column alongside, so
    the untransformed counts stay available regardless of what went into .raw.

    .X, .layers, .obsm, .varm and .obsp are untouched. Clusters, embeddings
    and existing DE results are unaffected.
    """

    if copy:
        adata = adata.copy()

    if isinstance(features, str):
        features = [features]

    if adata.raw is None:
        raise ValueError("This object has no .raw. Use append_to_X() or "
                         "add_to_obs() instead.")

    clash = [f for f in features if f in adata.raw.var_names]
    if len(clash) > 0:
        raise ValueError("Already in .raw: " + str(clash) +
                         ". Remove them first rather than duplicating.")

    aligned = align_to_adata(adata, frames, features, barcode_key=barcode_key,
                             batch_key=batch_key, verbose=verbose)

    scaled = _scale_values(aligned, features, adata, values, target_sum,
                           total_counts_key, verbose=verbose)

    # Rebuild .raw as old raw + the new columns
    raw_adata = adata.raw.to_adata()

    extra = ad.AnnData(
        X=sp.csr_matrix(scaled.to_numpy(dtype=np.float32)),
        obs=pd.DataFrame(index=raw_adata.obs_names),
        var=pd.DataFrame(index=list(scaled.columns)),
    )

    merged = ad.concat([raw_adata, extra], axis=1, merge="first")
    merged.obs = raw_adata.obs.copy()

    adata.raw = merged

    if verbose:
        print(".raw: " + str(raw_adata.n_vars) + " -> " + str(merged.n_vars) +
              " features (values='" + values + "')")

    # Keep the untransformed counts on obs as well
    if also_obs:
        for feature in features:
            adata.obs[str(feature) + "_counts"] = aligned[feature].values
        # End loop

    del raw_adata, extra, merged, aligned, scaled  # Cleanup

    return adata
# End of function


def append_to_X(adata, frames, features, values="counts", target_sum=2e4,
                total_counts_key="total_counts", barcode_key="gem_id",
                batch_key="batch", pad_varm=True, copy=False, verbose=True):
    """
    Append the construct feature(s) to .X and .var. Opt-in, because on a
    processed object .X is scaled z-scores over a HVG subset -- adding raw or
    log-normalized values puts two unit systems in one matrix, and anything
    that reads .X wholesale (PCA, HVG, DE on .X) will treat the new column as
    comparable when it is not.

    pad_varm : varm entries are indexed along the feature axis, so a stale
               varm["PCs"] would have the wrong shape after appending. With
               pad_varm=True the new feature gets a zero row in every varm
               matrix, which keeps shapes valid and gives it zero loading on
               the existing PCs -- correct, since it was not in the PCA.
               With False, varm is dropped.

    .layers is dropped if present, for the same shape reason -- a layer has no
    meaningful value for a feature that was not in the original run.
    """

    if copy:
        adata = adata.copy()

    if isinstance(features, str):
        features = [features]

    clash = [f for f in features if f in adata.var_names]
    if len(clash) > 0:
        raise ValueError("Already in .var: " + str(clash) + ".")

    aligned = align_to_adata(adata, frames, features, barcode_key=barcode_key,
                             batch_key=batch_key, verbose=verbose)

    scaled = _scale_values(aligned, features, adata, values, target_sum,
                           total_counts_key, verbose=verbose)

    if verbose and values == "counts":
        print("Note: appending raw counts to a scaled .X. Colour scales and "
              "any .X-wide operation will mix units.")

    new_values = scaled.to_numpy(dtype=np.float32)

    if sp.issparse(adata.X):
        new_X = sp.hstack([adata.X, sp.csr_matrix(new_values)], format="csr")
    else:
        new_X = np.hstack([np.asarray(adata.X), new_values])

    new_var = pd.concat(
        [adata.var, pd.DataFrame(index=list(scaled.columns))],
        axis=0)
    new_var["is_construct"] = new_var.index.isin(list(scaled.columns))

    out = ad.AnnData(X=new_X, obs=adata.obs.copy(), var=new_var)
    out.obsm = adata.obsm.copy()
    out.obsp = adata.obsp.copy()
    out.uns = adata.uns.copy()
    if adata.raw is not None:
        out.raw = adata.raw.to_adata()

    if pad_varm:
        n_new = len(scaled.columns)
        for key in adata.varm:
            block = np.asarray(adata.varm[key])
            out.varm[key] = np.vstack([block, np.zeros((n_new, block.shape[1]),
                                                       dtype=block.dtype)])
        # End loop
        del n_new  # Cleanup
    elif len(adata.varm) > 0 and verbose:
        print("Dropping varm: " + str(list(adata.varm)))

    layer_keys = [k for k in adata.layers if k is not None]
    if len(layer_keys) > 0 and verbose:
        print("Dropping layers: " + str(layer_keys))
    del layer_keys  # Cleanup

    if verbose:
        print(".X: " + str(adata.n_vars) + " -> " + str(out.n_vars) +
              " features (values='" + values + "')")

    del aligned, scaled, new_values, new_X, new_var  # Cleanup

    return out
# End of function


def add_to_obs(adata, frames, features, values="counts", target_sum=2e4,
               total_counts_key="total_counts", barcode_key="gem_id",
               batch_key="batch", suffix="_counts", copy=False, verbose=True):
    """
    Lowest-risk route: per-cell columns only, nothing on the feature axis.
    Enough for a dotplot row, a UMAP overlay or a per-cluster table, since
    sc.pl.* resolves names through sc.get.obs_df.
    """

    if copy:
        adata = adata.copy()

    if isinstance(features, str):
        features = [features]

    aligned = align_to_adata(adata, frames, features, barcode_key=barcode_key,
                             batch_key=batch_key, verbose=verbose)

    scaled = _scale_values(aligned, features, adata, values, target_sum,
                           total_counts_key, verbose=verbose)

    for feature in features:
        adata.obs[str(feature) + suffix] = scaled[feature].values
    # End loop

    del aligned, scaled  # Cleanup

    return adata
# End of function


# ===========================================================================
# 5. Readouts
# ===========================================================================

def construct_summary(adata, feature, groupby, counts_key=None,
                      endogenous=None):
    """
    Per-cluster table: percent of cells with any construct signal, mean count
    among positive cells, and -- if `endogenous` names the endogenous gene --
    the mean of that gene in construct-positive versus construct-negative
    cells.

    That last comparison is the one the reviewer question turns on.
    """

    if counts_key is None:
        counts_key = str(feature) + "_counts"

    counts = adata.obs[counts_key].astype(float)
    positive = counts > 0

    table = pd.DataFrame({
        "n_cells": adata.obs.groupby(groupby, observed=True).size(),
        "pct_positive": positive.groupby(
            adata.obs[groupby], observed=True).mean() * 100.0,
        "mean_count_in_positive": counts.where(positive).groupby(
            adata.obs[groupby], observed=True).mean(),
    })

    if endogenous is not None:
        # Read the endogenous gene from .raw when it exists, so the numbers
        # are log-normalized expression rather than scaled z-scores
        values = sc.get.obs_df(adata, keys=[endogenous],
                               use_raw=adata.raw is not None)[endogenous]
        table[endogenous + "_in_pos"] = values.where(positive).groupby(
            adata.obs[groupby], observed=True).mean()
        table[endogenous + "_in_neg"] = values.where(~positive).groupby(
            adata.obs[groupby], observed=True).mean()
        del values  # Cleanup

    return table.round(3)
# End of function


# ===========================================================================
# Example
# ===========================================================================

if __name__ == "__main__":

    GEM_DIRS = {
        "gem1": "source_data_realign/gem1/outs",
        "gem2": "source_data_realign/gem2/outs",
    }
    PROCESSED = "h5ad/<processed>.h5ad"

    adata = sc.read_h5ad(PROCESSED)

    # ---- 1. confirm the exact feature name in the new reference -----------
    print(list_features(GEM_DIRS["gem1"], "runx"))
    print("features unique to the re-run: " +
          str(diff_features(GEM_DIRS["gem1"], adata)))

    FEATURES = ["MycRunx3"]

    # ---- 2. extract from both GEMs ---------------------------------------
    frames = read_feature_counts_all(GEM_DIRS, FEATURES)

    # ---- 3. confirm the join keys ----------------------------------------
    detect_barcode_key(adata, frames)

    # ---- 4. append ---------------------------------------------------------
    append_to_raw(adata, frames, FEATURES, values="lognorm",
                  barcode_key="gem_id", batch_key="batch")

    # MycRunx3 now behaves like any other gene
    sc.pl.dotplot(adata, ["Runx3", "MycRunx3"], groupby="RPE",
                  standard_scale="var", swap_axes=True)

    print(construct_summary(adata, "MycRunx3", groupby="RPE",
                            endogenous="Runx3"))

    adata.write_h5ad(filename=PROCESSED.replace(".h5ad", "_construct.h5ad"),
                     compression="gzip", compression_opts=9)


# ---------------------------------------------------------------------------
# Multimapping caveat
# ---------------------------------------------------------------------------
#
# If the construct sequence in the new reference contains the Runx3 CDS, then
# reads from that region align equally well to the transgene and to the
# endogenous locus. Cell Ranger discards reads that map to more than one gene,
# so those reads are counted for NEITHER feature. The practical consequences:
#
#   - MycRunx3 counts reflect only the construct-unique portion of the
#     transcript (the Myc tag, vector UTRs, the LTR-proximal sequence), not
#     total construct expression. Low counts are expected and are not evidence
#     of low transgene expression.
#   - Endogenous Runx3 counts can be DEPRESSED in transduced cells relative to
#     the original alignment, because shared reads that previously counted
#     toward Runx3 now multimap and are dropped. Comparing endogenous Runx3
#     between the two alignments is therefore not apples-to-apples.
#
# Worth checking how the construct was added to the reference before the
# numbers go into a response letter. Whether the CDS is present in the custom
# contig, and whether the 3' end of the construct transcript is distinguishable
# within a 10x 3' read's reach, both determine what these counts can support.
# ---------------------------------------------------------------------------