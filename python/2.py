#!/usr/bin/env python3
# coding: utf-8

"""
Optimized Sequence-group aware phylogenetic label checker with dynamic margin.
Produces multiple CSV/XLSX outputs and an optionally precomputed patristic distance pickle.

Usage (example):
 python phylo_label_checker_full.py --tree "path/to/tree.tre" --meta "path/to/meta.xlsx" --outdir "out/folder" --max_cpus 4 --dist_matrix_pickle None
"""

# ------------------------- CONFIG (edit here) -------------------------
DEFAULT_TREE = r"G:\Paper\nema-Nanopore-Sequencing\pylogenetic\data\Ok-Selected\Seq\tree\18S_final.tre"
DEFAULT_META = r"G:\Paper\nema-Nanopore-Sequencing\pylogenetic\data\Ok-Selected\by_cp_f-h\18S.xlsx"
DEFAULT_OUTDIR = r"G:\Paper\nema-Nanopore-Sequencing\pylogenetic\data\Ok-Selected\Seq\tree\phylocheck1\18S"
# Column names
DEFAULT_ID_COL = "ID"
DEFAULT_GENUS_COL = "Main_Organism"
DEFAULT_FAMILY_COL = "Family"
DEFAULT_SEQ_COL = "RNA Sequence"
# Threshold / tuning
DEFAULT_FAMILY_MULTIPLIER = 2.0
DEFAULT_GENUS_MULTIPLIER = 2.0
# dynamic_margin = max(0.05 * min_other, 0.01 * tree_diameter)
# percent identity threshold for mapping non-tree reps
DEFAULT_MIN_PID_FLAG = 80.0
# Pairwise mode: if True and Biopython available uses globalxx alignment; otherwise k-mer fallback
DEFAULT_USE_PAIRWISE2 = True
# Show progress bars if tqdm available
USE_TQDM = True
# Maximum CPUs to use (limits parallel processing to this number of cores; set to 1 for single-threaded)
import os
DEFAULT_MAX_CPUS = os.cpu_count() or 1
# Path to precomputed distance matrix pickle file (set to None to compute on-the-fly)
DEFAULT_DIST_MATRIX_PICKLE = None  # e.g., r"G:\path\to\dist_matrix.pkl"
# ---------------------------------------------------------------------

import sys
import argparse
import hashlib
import math
import pickle
import copy
from collections import defaultdict
from typing import List, Dict, Tuple

# dependencies
try:
    import pandas as pd
    import numpy as np
except Exception:
    print("Please install required Python packages: pandas, numpy")
    raise

try:
    import dendropy
except Exception:
    print("Please install dendropy (pip install dendropy).")
    raise

# optional: biopython for pairwise alignment (globalxx)
USE_PAIRWISE2 = False
if DEFAULT_USE_PAIRWISE2:
    try:
        from Bio import pairwise2  # type: ignore
        import warnings
        from Bio import BiopythonDeprecationWarning
        warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)
        USE_PAIRWISE2 = True
    except Exception:
        USE_PAIRWISE2 = False
        print("Warning: Biopython not available; falling back to k-mer similarity approximation.")

# optional tqdm
if USE_TQDM:
    try:
        from tqdm import tqdm
        TQDM = tqdm
    except Exception:
        TQDM = lambda x, **k: x
else:
    TQDM = lambda x, **k: x

import concurrent.futures
import logging
# Set up logging for better error handling
logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')

# -------------------- Helper functions --------------------
def md5_hash(s: str) -> str:
    s2 = (s or "").upper().strip()
    return hashlib.md5(s2.encode("utf-8")).hexdigest()

def read_metadata(path: str, id_col: str, genus_col: str, family_col: str, seq_col: str) -> pd.DataFrame:
    df = pd.read_excel(path, dtype=str)
    missing_cols = [c for c in (id_col, genus_col, family_col, seq_col) if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in metadata Excel: {missing_cols}")
    df = df[[id_col, genus_col, family_col, seq_col]].copy()
    df.columns = ['ID', 'Genus', 'Family', 'Sequence']
    df['Sequence'] = df['Sequence'].fillna('').astype(str)
    return df

def percent_identity_pairwise(s1: str, s2: str) -> float:
    if not s1 or not s2:
        return 0.0
    if USE_PAIRWISE2:
        try:
            score = pairwise2.align.globalxx(s1, s2, one_alignment_only=True, score_only=True)
            pid = 100.0 * score / max(len(s1), len(s2))
            return float(pid)
        except Exception as ve:
            logging.warning(f"Pairwise alignment failed: {ve}. Falling back to k-mer.")
    # k-mer Jaccard fallback (k=6); note: this is an approximation of similarity, not true percent identity
    k = 6
    def kmers(s):
        sU = s.upper()
        if len(sU) < k:
            return {sU}
        return set(sU[i:i+k] for i in range(len(sU) - k + 1))
    a = kmers(s1)
    b = kmers(s2)
    if not a or not b:
        return 0.0
    inter = len(a & b)
    uni = len(a | b)
    return 100.0 * inter / uni

# Double-sweep to estimate tree diameter (fast)
def estimate_tree_diameter(tips_labels: List[str], cached_distance_func) -> float:
    """
    Double-sweep method:
     - pick an arbitrary tip A (tips[0])
     - find tip B farthest from A
     - find tip C farthest from B
     - diameter approximated by dist(B,C)
    This requires O(n) distance computations (twice).
    """
    if not tips_labels:
        return 0.0
    if len(tips_labels) == 1:
        return 0.0
    a = tips_labels[0]
    # distances from a to all tips
    max_d = -1.0
    far_b = a
    for t in tips_labels:
        d = cached_distance_func(a, t)
        if math.isnan(d):
            continue
        if d > max_d:
            max_d = d
            far_b = t
    # now from far_b find farthest
    max_d2 = -1.0
    far_c = far_b
    for t in tips_labels:
        d = cached_distance_func(far_b, t)
        if math.isnan(d):
            continue
        if d > max_d2:
            max_d2 = d
            far_c = t
    if max_d2 < 0:
        return 0.0
    return float(max_d2)

def cached_distance(a_label: str, b_label: str, dist_cache: Dict[Tuple[str, str], float]) -> float:
    if a_label == b_label:
        return 0.0
    key = (a_label, b_label) if a_label <= b_label else (b_label, a_label)
    return dist_cache.get(key, float('nan'))

def mean_dist_to_group_on_the_fly(id_: str, group_ids: List[str], tree_tips_all: set, dist_cache: Dict[Tuple[str, str], float]) -> float:
    vals = []
    for other in group_ids:
        if other == id_:
            continue
        if other not in tree_tips_all:
            continue
        d = cached_distance(id_, other, dist_cache)
        if not math.isnan(d):
            vals.append(d)
    if not vals:
        return float('nan')
    return float(sum(vals) / len(vals))

def compute_rep_metrics(rep, meta_rep_in_tree, rep_family, rep_genus, tree_tips_all, dist_cache):
    rows = meta_rep_in_tree[meta_rep_in_tree['Representative'] == rep]
    if rows.empty:
        fam = ''
        gen = ''
    else:
        fam = rows.iloc[0]['Family']
        gen = rows.iloc[0]['Genus']
    mean_dist_family = mean_dist_to_group_on_the_fly(rep, rep_family.get(fam, []), tree_tips_all, dist_cache)
    mean_dist_genus  = mean_dist_to_group_on_the_fly(rep, rep_genus.get(gen, []), tree_tips_all, dist_cache)
    return {
        'Representative': rep,
        'Family': fam,
        'Genus': gen,
        'mean_dist_family': mean_dist_family,
        'mean_dist_genus': mean_dist_genus
    }

def compute_family_threshold(fam_reps, args, tree_tips_all, dist_cache):
    fam, reps = fam_reps
    vals = [mean_dist_to_group_on_the_fly(r, reps, tree_tips_all, dist_cache) for r in reps]
    arr = np.asarray([v for v in vals if not math.isnan(v)])
    if arr.size == 0:
        return fam, float('nan'), float('nan'), float('nan'), len(reps)
    else:
        median = float(np.median(arr))
        std = float(np.std(arr))
        threshold = median + args.family_multiplier * std
        return fam, threshold, median, std, len(reps)

def compute_genus_threshold(gen_reps, args, tree_tips_all, dist_cache):
    gen, reps = gen_reps
    vals = [mean_dist_to_group_on_the_fly(r, reps, tree_tips_all, dist_cache) for r in reps]
    arr = np.asarray([v for v in vals if not math.isnan(v)])
    if arr.size == 0:
        return gen, float('nan'), float('nan'), float('nan'), len(reps)
    else:
        median = float(np.median(arr))
        std = float(np.std(arr))
        threshold = median + args.genus_multiplier * std
        return gen, threshold, median, std, len(reps)

def flag_rep(row, family_thresholds, genus_thresholds, families_list, rep_family, tree_diameter, tree_tips_all, dist_cache):
    rep = row['Representative']
    fam = row['Family']
    gen = row['Genus']
    md_f = row['mean_dist_family']
    md_g = row['mean_dist_genus']
    reasons = []
    fam_thr = family_thresholds.get(fam, float('nan'))
    gen_thr = genus_thresholds.get(gen, float('nan'))
    if not math.isnan(md_f) and not math.isnan(fam_thr) and md_f > fam_thr:
        reasons.append(f"high_mean_dist_to_family (>{fam_thr:.4g})")
    if not math.isnan(md_g) and not math.isnan(gen_thr) and md_g > gen_thr:
        reasons.append(f"high_mean_dist_to_genus (>{gen_thr:.4g})")
    # compute min_other (closest other family)
    min_other = None
    min_other_fam = None
    for other_f in families_list:
        if other_f == fam:
            continue
        other_reps = rep_family.get(other_f, [])
        if not other_reps:
            continue
        val = mean_dist_to_group_on_the_fly(rep, other_reps, tree_tips_all, dist_cache)
        if math.isnan(val):
            continue
        if (min_other is None) or (val < min_other):
            min_other = val
            min_other_fam = other_f
    # dynamic margin per representative:
    # use max(0.05 * min_other, 0.01 * tree_diameter)
    if min_other is None or math.isnan(min_other):
        dynamic_margin = 0.01 * tree_diameter
    else:
        dynamic_margin = max(0.05 * min_other, 0.01 * tree_diameter)
    # decide "closer to other family" using dynamic margin
    if (min_other is not None) and (not math.isnan(md_f)) and (min_other + dynamic_margin < md_f):
        reasons.append(f"closer_to_family_{min_other_fam} (dist {min_other:.4g} < own {md_f:.4g}; margin {dynamic_margin:.4g})")
    if reasons:
        return {
            'Representative': rep,
            'Family': fam,
            'Genus': gen,
            'reasons': "; ".join(reasons),
            'mean_dist_family': md_f,
            'mean_dist_genus': md_g,
            'min_other_dist': (min_other if min_other is not None else float('nan')),
            'dynamic_margin_used': dynamic_margin
        }
    return None

def map_rep(rep, rep_to_seq, in_tree_rep_seqs, args, suspects_rep_df, rep_to_members, meta):
    seq = rep_to_seq.get(rep, '')
    best_pid = -1.0
    best_inrep = None
    for inrep, inseq in in_tree_rep_seqs.items():
        pid = percent_identity_pairwise(seq, inseq)
        if pid > best_pid:
            best_pid = pid
            best_inrep = inrep
    if best_inrep and best_pid >= args.min_pid_flag:
        row = suspects_rep_df[suspects_rep_df['Representative'] == best_inrep]
        if not row.empty:
            reasons = row.iloc[0]['reasons']
            fam = row.iloc[0]['Family']
            gen = row.iloc[0]['Genus']
            members = rep_to_members.get(rep, [])
            mapped = []
            for m in members:
                seq_group_row = meta.loc[meta['ID'] == m, 'SequenceGroup']
                seq_group = seq_group_row.values[0] if not seq_group_row.empty else ''
                mapped.append({
                    'ID': m,
                    'SequenceGroup': seq_group,
                    'Representative': rep,
                    'MappedToRepresentative': best_inrep,
                    'MappedPID': best_pid,
                    'Family': fam,
                    'Genus': gen,
                    'reasons': f"mapped_to_suspicious_rep {best_inrep} (pid {best_pid:.1f}%) -> {reasons}"
                })
            return mapped
    return []

def compute_pairwise_family_stats(fam_a, fam_b, rep_family, tree_tips_all, dist_cache):
    reps_a = rep_family.get(fam_a, [])
    reps_b = rep_family.get(fam_b, [])
    if not reps_a or not reps_b:
        return None
    vals = []
    for r_a in reps_a:
        for r_b in reps_b:
            if r_a == r_b:
                continue
            d = cached_distance(r_a, r_b, dist_cache)
            if not math.isnan(d):
                vals.append(d)
    if not vals:
        return None
    arr = np.asarray(vals)
    return {
        'FamilyA': fam_a,
        'FamilyB': fam_b,
        'Mean_distance': np.mean(arr),
        'Min_distance': np.min(arr),
        'Max_distance': np.max(arr)
    }

# -------------------- Main optimized flow --------------------
def main(args):
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    meta_basename = os.path.basename(args.meta)
    meta_file = os.path.splitext(meta_basename)[0]
    print("Reading metadata...")
    meta = read_metadata(args.meta, args.id_col, args.genus_col, args.family_col, args.seq_col)
    total_rows = len(meta)
    print("Grouping identical sequences by MD5 hash...")
    meta['SeqHash'] = meta['Sequence'].apply(md5_hash)
    hash2group = {}
    group_counter = 0
    for h in meta['SeqHash'].unique():
        group_counter += 1
        label = f"SEQGRP_{group_counter:06d}"
        hash2group[h] = label
    meta['SequenceGroup'] = meta['SeqHash'].map(hash2group)
    n_groups = meta['SequenceGroup'].nunique()
    print(f"Metadata rows: {total_rows}, unique identical-sequence groups: {n_groups}")
    seqgrp_to_ids = meta.groupby('SequenceGroup')['ID'].apply(list).to_dict()
    id_to_seq = meta.set_index('ID')['Sequence'].to_dict()
    print("Loading tree (dendropy)...")
    tree = dendropy.Tree.get(path=args.tree, schema="newick")
    tree_tips_all = {leaf.taxon.label for leaf in tree.leaf_node_iter() if getattr(leaf, "taxon", None) and getattr(leaf.taxon, "label", None)}
    print(f"Tree tips: {len(tree_tips_all)}")
    meta_ids = set(meta['ID'].tolist())
    tips_in_meta = sorted(list(tree_tips_all & meta_ids))
    missing_in_tree = sorted(list(meta_ids - tree_tips_all))
    missing_in_meta = sorted(list(tree_tips_all - meta_ids))
    print(f"IDs in metadata & tree intersection: {len(tips_in_meta)}; missing in tree: {len(missing_in_tree)}")
    print("Selecting representative per SequenceGroup (prefer tree tips)...")
    seqgroup_to_rep: Dict[str, str] = {}
    for grp, ids in seqgrp_to_ids.items():
        rep = None
        for _id in ids:
            if _id in tree_tips_all:
                rep = _id
                break
        if rep is None:
            rep = max(ids, key=lambda x: len(id_to_seq.get(x, '')))
        seqgroup_to_rep[grp] = rep
    id_to_rep = {i: seqgroup_to_rep[meta.loc[meta['ID'] == i, 'SequenceGroup'].values[0]] for i in meta['ID']}
    meta['Representative'] = meta['ID'].map(id_to_rep)
    reps_all = set(seqgroup_to_rep.values())
    reps_in_tree = sorted([r for r in reps_all if r in tree_tips_all])
    reps_not_in_tree = sorted([r for r in reps_all if r not in tree_tips_all])
    print(f"Representatives total: {len(reps_all)}, reps-in-tree: {len(reps_in_tree)}, reps-not-in-tree: {len(reps_not_in_tree)}")
    rep_to_seq = {rep: id_to_seq.get(rep, '') for rep in reps_all}
    meta_rep_in_tree = meta[meta['Representative'].isin(reps_in_tree)].copy()
    rep_family = meta_rep_in_tree.groupby('Family')['Representative'].unique().to_dict()
    rep_genus  = meta_rep_in_tree.groupby('Genus')['Representative'].unique().to_dict()
    rep_family = {k: list(v) for k, v in rep_family.items()}
    rep_genus  = {k: list(v) for k, v in rep_genus.items()}
    if len(reps_in_tree) < 2:
        print("Not enough representatives present in the tree to compute patristic statistics. Exiting.")
        meta.to_csv(os.path.join(outdir, meta_file + '_metadata_with_seqgroups.csv'), index=False)
        return

    print("Preparing patristic distance matrix...")
    # determine pickle path
    dist_matrix_pickle = args.dist_matrix_pickle if args.dist_matrix_pickle and args.dist_matrix_pickle.lower() != 'none' else os.path.join(outdir, meta_file + '_patristic_dist_matrix.pkl')
    taxon_label_to_taxon = {t.label: t for t in tree.taxon_namespace if getattr(t, "label", None)}
    dist_cache: Dict[Tuple[str, str], float] = {}
    if os.path.exists(dist_matrix_pickle):
        print(f"Loading precomputed distance matrix from {dist_matrix_pickle}...")
        with open(dist_matrix_pickle, 'rb') as f:
            dist_cache = pickle.load(f)
    else:
        print("Computing phylogenetic distance matrix (dendropy)...")
        pdm = tree.phylogenetic_distance_matrix()
        print("Precomputing all pairwise distances...")
        taxon_objects = [taxon_label_to_taxon[label] for label in reps_in_tree if label in taxon_label_to_taxon]
        # compute upper triangle
        for i, taxon_a in enumerate(TQDM(taxon_objects, desc="precompute dists")):
            for taxon_b in taxon_objects[i+1:]:
                try:
                    d = pdm.patristic_distance(taxon_a, taxon_b)
                    key = (taxon_a.label, taxon_b.label) if taxon_a.label <= taxon_b.label else (taxon_b.label, taxon_a.label)
                    dist_cache[key] = float(d)
                except Exception as e:
                    logging.warning(f"Error computing distance {getattr(taxon_a,'label',None)}-{getattr(taxon_b,'label',None)}: {e}")
        try:
            print(f"Saving precomputed distance matrix to {dist_matrix_pickle}...")
            with open(dist_matrix_pickle, 'wb') as f:
                pickle.dump(dist_cache, f)
        except Exception as e:
            logging.warning(f"Could not save dist matrix pickle: {e}")

    print("Preparing on-demand patristic distance cache (using precomputed dict)...")
    # estimate tree diameter using double-sweep (efficient)
    print("Estimating tree diameter (double-sweep)...")
    tree_diameter = estimate_tree_diameter(reps_in_tree, lambda a, b: cached_distance(a, b, dist_cache))
    print(f"Estimated tree diameter (patristic): {tree_diameter:.6g}")

    print("Computing mean patristic distances for representatives (on-demand, cached, parallelized)...")
    rep_metrics = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_cpus) as executor:
        futures = [executor.submit(compute_rep_metrics, rep, meta_rep_in_tree, rep_family, rep_genus, tree_tips_all, dist_cache) for rep in reps_in_tree]
        for future in TQDM(concurrent.futures.as_completed(futures), total=len(reps_in_tree), desc="reps"):
            rep_metrics.append(future.result())
    rep_metrics_df = pd.DataFrame(rep_metrics)

    # compute dynamic thresholds per family/genus (parallelize if many groups)
    print("Computing dynamic thresholds...")
    family_thresholds = {}
    family_summary_rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_cpus) as executor:
        futures = [executor.submit(compute_family_threshold, (fam, reps), args, tree_tips_all, dist_cache) for fam, reps in rep_family.items()]
        for future in TQDM(concurrent.futures.as_completed(futures), total=len(rep_family), desc="family thresholds"):
            fam, thr, med, std, count = future.result()
            family_thresholds[fam] = thr
            family_summary_rows.append({
                'Family': fam,
                'Family_threshold': thr,
                'Median_dist': med,
                'Std_dev': std,
                'Member_count': count
            })

    genus_thresholds = {}
    genus_summary_rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_cpus) as executor:
        futures = [executor.submit(compute_genus_threshold, (gen, reps), args, tree_tips_all, dist_cache) for gen, reps in rep_genus.items()]
        for future in TQDM(concurrent.futures.as_completed(futures), total=len(rep_genus), desc="genus thresholds"):
            gen, thr, med, std, count = future.result()
            genus_thresholds[gen] = thr
            genus_summary_rows.append({
                'Genus': gen,
                'Genus_threshold': thr,
                'Median_dist': med,
                'Std_dev': std,
                'Member_count': count
            })

    family_summary_df = pd.DataFrame(family_summary_rows)
    genus_summary_df = pd.DataFrame(genus_summary_rows)
    rep_thresholds_df = pd.concat([family_summary_df, genus_summary_df], axis=1, ignore_index=False)

    print("Flagging suspicious representatives using dynamic margin (parallelized)...")
    suspect_reps = []
    families_list = list(rep_family.keys())
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_cpus) as executor:
        futures = [executor.submit(flag_rep, row, family_thresholds, genus_thresholds, families_list, rep_family, tree_diameter, tree_tips_all, dist_cache) for _, row in rep_metrics_df.iterrows()]
        for future in TQDM(concurrent.futures.as_completed(futures), total=len(rep_metrics_df), desc="flagging reps"):
            result = future.result()
            if result:
                suspect_reps.append(result)
    suspects_rep_df = pd.DataFrame(suspect_reps)

    # Propagate flags to members
    print("Propagating flags from representatives to all group members...")
    rep_to_members = defaultdict(list)
    for grp, rep in seqgroup_to_rep.items():
        rep_to_members[rep].extend(seqgrp_to_ids[grp])
    member_suspect_rows = []
    for _, r in suspects_rep_df.iterrows():
        rep = r['Representative']
        members = rep_to_members.get(rep, [])
        for m in members:
            seq_group_row = meta.loc[meta['ID'] == m, 'SequenceGroup']
            seq_group = seq_group_row.values[0] if not seq_group_row.empty else ''
            member_suspect_rows.append({
                'ID': m,
                'SequenceGroup': seq_group,
                'Representative': rep,
                'Family': r['Family'],
                'Genus': r['Genus'],
                'reasons': r['reasons'],
                'mean_dist_family': r['mean_dist_family'],
                'mean_dist_genus': r['mean_dist_genus'],
                'min_other_dist': r.get('min_other_dist', float('nan')),
                'dynamic_margin_used': r.get('dynamic_margin_used', float('nan'))
            })
    member_suspects_df = pd.DataFrame(member_suspect_rows)

    # Map non-tree representatives to in-tree reps by sequence similarity
    print("Mapping representatives NOT in tree to nearest in-tree rep by sequence similarity (pid)...")
    in_tree_rep_seqs = {rep: rep_to_seq.get(rep, '') for rep in reps_in_tree}
    mapped_rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_cpus) as executor:
        futures = [executor.submit(map_rep, rep, rep_to_seq, in_tree_rep_seqs, args, suspects_rep_df, rep_to_members, meta) for rep in reps_not_in_tree]
        for future in TQDM(concurrent.futures.as_completed(futures), total=len(reps_not_in_tree), desc="map_reps"):
            results = future.result()
            mapped_rows.extend(results)
    mapped_flags_df = pd.DataFrame(mapped_rows)

    # Combine
    if not member_suspects_df.empty and not mapped_flags_df.empty:
        mapped_df_norm = mapped_flags_df.rename(columns={'MappedToRepresentative': 'MappedTo', 'MappedPID': 'MappedPID'})
        suspects_final_df = pd.concat([member_suspects_df, mapped_df_norm], sort=False, ignore_index=True)
    elif not member_suspects_df.empty:
        suspects_final_df = member_suspects_df.copy()
    elif not mapped_flags_df.empty:
        suspects_final_df = mapped_flags_df.copy()
    else:
        suspects_final_df = pd.DataFrame(columns=['ID', 'SequenceGroup', 'Representative', 'Family', 'Genus', 'reasons'])

    # inconsistent groups (identical sequence but different labels)
    print("Checking identical-sequence groups for inconsistent labels...")
    inconsistent_groups = []
    for grp, ids in seqgrp_to_ids.items():
        if len(ids) <= 1:
            continue
        labels = meta[meta['ID'].isin(ids)][['Family', 'Genus']].drop_duplicates()
        if len(labels) > 1:
            inconsistent_groups.append({
                'SequenceGroup': grp,
                'member_count': len(ids),
                'members': ";".join(ids),
                'distinct_labels': labels.to_dict(orient='records')
            })
    inconsistent_df = pd.DataFrame(inconsistent_groups)

    # Create clean tree without suspects
    print("Creating clean tree without suspicious representatives...")
    sys.setrecursionlimit(10000)
    try:
        clean_tree = copy.deepcopy(tree)
        suspect_labels = suspects_rep_df['Representative'].tolist()
        # prune_taxa_with_labels removes taxa with provided labels (keeps others)
        clean_tree.prune_taxa_with_labels(suspect_labels)
        clean_tree_path = os.path.join(outdir, meta_file + '_clean.tre')
        clean_tree.write(path=clean_tree_path, schema="newick")
    except Exception as e:
        logging.warning(f"Could not produce clean tree: {e}")
        clean_tree_path = None

    # Pairwise similarity stats
    print("Computing pairwise family stats...")
    pairwise_rows = []
    families_sorted = sorted(families_list)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_cpus) as executor:
        futures = []
        for i, fam_a in enumerate(families_sorted):
            for fam_b in families_sorted[i:]:  # include within family
                futures.append(executor.submit(compute_pairwise_family_stats, fam_a, fam_b, rep_family, tree_tips_all, dist_cache))
        for future in TQDM(concurrent.futures.as_completed(futures), total=len(futures), desc="pairwise family stats"):
            result = future.result()
            if result:
                pairwise_rows.append(result)
    pairwise_stats_df = pd.DataFrame(pairwise_rows)

    # Tree tip presence summary
    print("Creating tree tip presence summary...")
    tip_presence_rows = []
    tip_presence_rows.append({'Category': 'Intersection', 'Count': len(tips_in_meta), 'IDs': ";".join(tips_in_meta)})
    tip_presence_rows.append({'Category': 'Missing in tree', 'Count': len(missing_in_tree), 'IDs': ";".join(missing_in_tree)})
    tip_presence_rows.append({'Category': 'Missing in metadata', 'Count': len(missing_in_meta), 'IDs': ";".join(missing_in_meta)})
    tip_presence_df = pd.DataFrame(tip_presence_rows)

    # Suspect statistics
    print("Creating suspect statistics...")
    if not suspects_rep_df.empty:
        suspect_stats = suspects_rep_df.groupby('Family').agg(
            Num_suspects=('Representative', 'count'),
            Mean_distance_over_threshold=('mean_dist_family', 'mean')
        ).reset_index()
    else:
        suspect_stats = pd.DataFrame(columns=['Family', 'Num_suspects', 'Mean_distance_over_threshold'])
    family_counts = pd.DataFrame([{'Family': fam, 'Total': len(reps)} for fam, reps in rep_family.items()])
    suspect_stats = pd.merge(family_counts, suspect_stats, on='Family', how='left').fillna(0)
    suspect_stats['Num_suspects'] = suspect_stats['Num_suspects'].astype(int)
    suspect_stats['Percent_suspect'] = (suspect_stats['Num_suspects'] / suspect_stats['Total']) * 100

    # Save outputs
    print("Saving outputs to:", outdir)
    try:
        meta.to_csv(os.path.join(outdir, meta_file + '_metadata_with_seqgroups.csv'), index=False)
    except Exception as e:
        logging.warning(f"Failed to write metadata CSV: {e}")

    try:
        suspects_final_df.to_csv(os.path.join(outdir, meta_file + '_suspect_sequences_propagated.csv'), index=False)
    except Exception as e:
        logging.warning(f"Failed to write suspects CSV: {e}")

    try:
        seqgroups_df = pd.DataFrame({
            'SequenceGroup': list(seqgrp_to_ids.keys()),
            'members': [";".join(seqgrp_to_ids[g]) for g in seqgrp_to_ids.keys()],
            'representative': [seqgroup_to_rep[g] for g in seqgrp_to_ids.keys()]
        })
        seqgroups_df.to_csv(os.path.join(outdir, meta_file + '_sequence_groups_summary.csv'), index=False)
    except Exception as e:
        logging.warning(f"Failed to write sequence groups CSV: {e}")

    try:
        inconsistent_df.to_csv(os.path.join(outdir, meta_file + '_inconsistent_sequence_groups.csv'), index=False)
    except Exception as e:
        logging.warning(f"Failed to write inconsistent groups CSV: {e}")

    try:
        rep_thresholds_df.to_csv(os.path.join(outdir, meta_file + '_representative_thresholds_summary.csv'), index=False)
    except Exception as e:
        logging.warning(f"Failed to write representative thresholds CSV: {e}")

    try:
        pairwise_stats_df.to_csv(os.path.join(outdir, meta_file + '_pairwise_similarity_stats.csv'), index=False)
    except Exception as e:
        logging.warning(f"Failed to write pairwise stats CSV: {e}")

    try:
        tip_presence_df.to_csv(os.path.join(outdir, meta_file + '_tree_tip_presence_summary.csv'), index=False)
    except Exception as e:
        logging.warning(f"Failed to write tip presence CSV: {e}")

    try:
        suspect_stats.to_csv(os.path.join(outdir, meta_file + '_suspect_statistics.csv'), index=False)
    except Exception as e:
        logging.warning(f"Failed to write suspect statistics CSV: {e}")

    excel_path = os.path.join(outdir, meta_file + '_phylo_seqgroup_check_summary_dynamic_margin.xlsx')
    try:
        with pd.ExcelWriter(excel_path) as writer:
            meta.to_excel(writer, sheet_name='metadata_with_seqgroups', index=False)
            seqgroups_df.to_excel(writer, sheet_name='sequence_groups', index=False)
            suspects_final_df.to_excel(writer, sheet_name='suspects_propagated', index=False)
            inconsistent_df.to_excel(writer, sheet_name='inconsistent_groups', index=False)
            rep_metrics_df.to_excel(writer, sheet_name='rep_metrics', index=False)
            pd.DataFrame([{
                'total_meta_rows': total_rows,
                'unique_sequence_groups': n_groups,
                'representatives_total': len(reps_all),
                'reps_in_tree': len(reps_in_tree),
                'reps_not_in_tree': len(reps_not_in_tree),
                'suspect_records_count': len(suspects_final_df),
                'inconsistent_groups_count': len(inconsistent_df),
                'estimated_tree_diameter': tree_diameter
            }]).to_excel(writer, sheet_name='summary', index=False)
    except Exception as e:
        logging.warning(f"Failed to write Excel summary: {e}")

    print("Done. Files written to:", os.path.abspath(outdir))
    print("Unique distance cache entries computed:", len(dist_cache))
    if clean_tree_path:
        print("Clean tree path:", clean_tree_path)

# -------------------- CLI --------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Optimized Sequence-group aware phylogenetic label checker with dynamic margin.")
    parser.add_argument('--tree', default=DEFAULT_TREE, help='Path to tree file (Newick).')
    parser.add_argument('--meta', default=DEFAULT_META, help='Path to metadata Excel (columns: ID, Main_Organism, Family, RNA Sequence).')
    parser.add_argument('--outdir', default=DEFAULT_OUTDIR, help='Output folder.')
    parser.add_argument('--id_col', default=DEFAULT_ID_COL, help='ID column in metadata')
    parser.add_argument('--genus_col', default=DEFAULT_GENUS_COL, help='Genus column in metadata')
    parser.add_argument('--family_col', default=DEFAULT_FAMILY_COL, help='Family column in metadata')
    parser.add_argument('--seq_col', default=DEFAULT_SEQ_COL, help='Sequence column in metadata')
    parser.add_argument('--family_multiplier', type=float, default=DEFAULT_FAMILY_MULTIPLIER)
    parser.add_argument('--genus_multiplier', type=float, default=DEFAULT_GENUS_MULTIPLIER)
    parser.add_argument('--min_pid_flag', type=float, default=DEFAULT_MIN_PID_FLAG)
    parser.add_argument('--max_cpus', type=int, default=DEFAULT_MAX_CPUS, help='Maximum number of CPU cores to use for parallel processing.')
    parser.add_argument('--dist_matrix_pickle', default=DEFAULT_DIST_MATRIX_PICKLE, help='Path to precomputed distance matrix pickle file ("None" or omitted to store under outdir).')
    args = parser.parse_args()
    # propagate CLI into args (to keep original variable names)
    args.family_multiplier = args.family_multiplier
    args.genus_multiplier = args.genus_multiplier
    args.min_pid_flag = args.min_pid_flag
    # run
    main(args)
