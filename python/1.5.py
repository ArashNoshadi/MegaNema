
from __future__ import annotations
import argparse
import os
import json
import re
import math
import numpy as np
import pandas as pd
from Bio.Seq import Seq

# -------------------------
# Defaults (edit if needed)
# -------------------------
DEFAULT_INPUT = r'G:\Paper\nema-Nanopore-Sequencing\pylogenetic\data\Ok-Selected\by_cp_f-h\COX1.xlsx'
DEFAULT_OUTPUT_DIR = r'G:\Paper\nema-Nanopore-Sequencing\pylogenetic\data\Ok-Selected\by_cp_f-h\NUMTs'

_non_dna_re = re.compile(r'[^ATCGUatcgu]')

# -------------------------
# Helpers
# -------------------------
def clean_sequence_for_translation(seq: str) -> str:
    if not isinstance(seq, str):
        return ""
    s = seq.replace(" ", "").upper().replace("U", "T")
    s = _non_dna_re.sub("", s)
    return s

def compute_gc(seq: str) -> float:
    if not seq:
        return float('nan')
    s = seq.upper()
    a = s.count('A'); t = s.count('T'); g = s.count('G'); c = s.count('C')
    atgc = a + t + g + c
    if atgc == 0:
        return float('nan')
    return (g + c) / atgc * 100.0

def detect_gc_outliers(df: pd.DataFrame, col="gc"):
    """
    Return (mean, std, lower, upper, mask_series)
    mask_series boolean: True => outlier
    """
    if col not in df.columns:
        return (float('nan'), float('nan'), float('nan'), float('nan'),
                pd.Series(False, index=df.index))
    try:
        gc_series = pd.to_numeric(df[col], errors='coerce').dropna()
        if gc_series.empty:
            return (float('nan'), float('nan'), float('nan'), float('nan'),
                    pd.Series(False, index=df.index))
        mean_gc = float(gc_series.mean())
        std_gc = float(gc_series.std(ddof=0))
        lower = mean_gc - 2 * std_gc
        upper = mean_gc + 2 * std_gc
        mask = (pd.to_numeric(df[col], errors='coerce') < lower) | (pd.to_numeric(df[col], errors='coerce') > upper)
        return (mean_gc, std_gc, lower, upper, mask.fillna(False))
    except Exception:
        return (float('nan'), float('nan'), float('nan'), float('nan'),
                pd.Series(False, index=df.index))

def translate_frame(subseq: str, table: int) -> str:
    if len(subseq) < 3:
        return ""
    try:
        return str(Seq(subseq).translate(table=table, to_stop=False))
    except Exception:
        return ""

def translate_six_frames(seq: str, table: int = 5):
    seq = seq or ""
    out = []
    for f in range(3):
        sub = seq[f:]
        trim = (len(sub) // 3) * 3
        prot = translate_frame(sub[:trim], table) if trim >= 3 else ""
        out.append({"frame": f"+{f+1}", "protein": prot})
    try:
        rc = str(Seq(seq).reverse_complement())
    except Exception:
        rc = ""
    for f in range(3):
        sub = rc[f:]
        trim = (len(sub) // 3) * 3
        prot = translate_frame(sub[:trim], table) if trim >= 3 else ""
        out.append({"frame": f"-{f+1}", "protein": prot})
    return out

def orf_len_before_first_stop(prot: str) -> int:
    if not prot:
        return 0
    i = prot.find('*')
    return len(prot) if i == -1 else i

def has_internal_stop(prot: str) -> bool:
    if prot is None:
        return True
    return '*' in prot[:-1]

def detect_ipsc_from_frames(frames_json_str):
    """Return 'Present' / 'Absent' / 'Unknown' based on internal stops in frames_summary JSON."""
    try:
        frames = json.loads(frames_json_str)
        if not isinstance(frames, list) or len(frames) == 0:
            return "Unknown"
        for f in frames:
            prot = f.get("protein", "")
            if prot and ('*' in prot[:-1] or f.get("internal_stop") is True):
                return "Present"
        return "Absent"
    except Exception:
        return "Unknown"

def evaluate_sequence_sixframes(
    seq_raw,
    min_aa_abs: int = 30,
    min_aa_pct: float = 0.80,
    table: int = 5
):
    seq = clean_sequence_for_translation(seq_raw)
    gc = compute_gc(seq)
    len_nt = len(seq)
    if len_nt < 3:
        return dict(
            clean_seq=seq, gc=gc, frames_summary=json.dumps([]),
            best_frame=None, best_prot="", best_prot_len=0, best_orf_len=0,
            decision_threshold=max(min_aa_abs, 0),
            numt_flag=True, final_reason="too_short_after_cleaning"
        )
    aa_pct_thresh = math.ceil((len_nt / 3.0) * min_aa_pct)
    decision_thresh = max(min_aa_abs, aa_pct_thresh)
    frames = translate_six_frames(seq, table=table)
    frames_info = []
    for f in frames:
        prot = f["protein"]
        info = {
            "frame": f["frame"],
            "protein_len": len(prot),
            "orf_len": orf_len_before_first_stop(prot),
            "internal_stop": has_internal_stop(prot),
            "protein": prot
        }
        frames_info.append(info)
    clean_hits = [fi for fi in frames_info if (not fi["internal_stop"]) and (fi["orf_len"] >= decision_thresh)]
    if clean_hits:
        best = max(clean_hits, key=lambda x: x["orf_len"])
        numt_flag = False
        reason = "clean_frame_meets_absolute_or_relative_threshold"
    else:
        best = max(frames_info, key=lambda x: x["orf_len"])
        numt_flag = True
        reason = "no_frame_meets_threshold"
    return dict(
        clean_seq=seq,
        gc=gc,
        frames_summary=json.dumps(frames_info, ensure_ascii=False),
        best_frame=best["frame"] if best else None,
        best_prot=best["protein"] if best else "",
        best_prot_len=best["protein_len"] if best else 0,
        best_orf_len=best["orf_len"] if best else 0,
        decision_threshold=decision_thresh,
        numt_flag=numt_flag,
        final_reason=reason
    )

# -------------------------
# I/O helpers
# -------------------------
def load_table(path: str) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xls", ".xlsx"):
        return pd.read_excel(path, engine="openpyxl")
    else:
        try:
            return pd.read_csv(path)
        except Exception:
            return pd.read_csv(path, sep="\t")

def detect_sequence_column(df: pd.DataFrame):
    candidates = ['RNA Sequence', 'Sequence', 'sequence', 'rna sequence',
                  'RNA_Sequence', 'RNAsequence', 'seq', 'Seq', 'clean_seq']
    for c in candidates:
        if c in df.columns:
            return c
    # fallback: any column with average length > 20
    for c in df.columns:
        try:
            if df[c].astype(str).map(len).mean() > 20:
                return c
        except Exception:
            continue
    return None

# default bins (C1..C5)
DEFAULT_SIZE_BINS = [
    (100, 150, "C1(100–150)"),
    (151, 300, "C2(151–300)"),
    (301, 450, "C3(301–450)"),
    (451, 600, "C4(451–600)"),
    (601, 661, "C5(601–661)")
]

# -------------------------
# Main
# -------------------------
def main():
    p = argparse.ArgumentParser(description="NUMT detection + plotting-ready length tables (no BLAST).")
    p.add_argument("--input", default=DEFAULT_INPUT, help="input xlsx/csv")
    p.add_argument("--output_dir", default=DEFAULT_OUTPUT_DIR, help="output directory")
    p.add_argument("--seq_col", default="", help="sequence column name (auto-detect if empty)")
    p.add_argument("--min_aa_abs", type=int, default=30)
    p.add_argument("--min_aa_pct", type=float, default=0.80)
    p.add_argument("--trans_table", type=int, default=5)
    p.add_argument("--gc_outliers", action="store_true")
    p.add_argument("--drop_protein_strings", action="store_true")
    p.add_argument("--size_bins_json", default="", help="optional JSON file defining bins [[low,high,label],...]")
    p.add_argument("--bin_width", type=int, default=0, help="if >0, create uniform bins of this width (bp)")
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = load_table(args.input)

    seq_col = args.seq_col if args.seq_col else detect_sequence_column(df)
    if not seq_col or seq_col not in df.columns:
        raise ValueError(f"Sequence column not found. Use --seq_col to set it. Columns: {list(df.columns)}")

    # run evaluation
    results = df[seq_col].apply(lambda s: evaluate_sequence_sixframes(
        s, min_aa_abs=args.min_aa_abs, min_aa_pct=args.min_aa_pct, table=args.trans_table
    ))
    res_df = pd.DataFrame(results.tolist())

    # optionally drop long protein strings in frames_summary to save space
    if args.drop_protein_strings and "frames_summary" in res_df:
        slim = []
        for js in res_df["frames_summary"]:
            try:
                frames = json.loads(js)
                for f in frames:
                    f.pop("protein", None)
                slim.append(json.dumps(frames, ensure_ascii=False))
            except Exception:
                slim.append(js)
        res_df["frames_summary"] = slim

    out_df = pd.concat([df.reset_index(drop=True), res_df], axis=1)

    # GC outliers
    mean_gc = std_gc = lower = upper = float('nan')
    if args.gc_outliers:
        mean_gc, std_gc, lower, upper, mask = detect_gc_outliers(out_df, col="gc")
        out_df["GC_outlier"] = mask
    else:
        out_df["GC_outlier"] = False

    out_df["Final_remove_flag"] = out_df["numt_flag"] | out_df["GC_outlier"]

    # ensure Length_Check present
    if 'Length_Check' not in out_df.columns:
        out_df['Length_Check'] = -1
    out_df['Coverage_label'] = out_df['Length_Check'].apply(lambda x: "High" if x == 1 else ("Low" if x == 0 else "Unknown"))

    # infer or normalize species/genus/family columns for summaries
    def infer_species(df):
        if 'Species' in df.columns:
            return df['Species'].astype(str)
        if 'Organism' in df.columns:
            return df['Organism'].astype(str).map(lambda s: " ".join(s.split()[0:2]) if isinstance(s, str) and len(s.split())>=2 else s)
        return pd.Series([""] * len(df), index=df.index)

    def infer_genus(df):
        if 'Genus' in df.columns:
            return df['Genus'].astype(str)
        if 'Main_Organism' in df.columns:
            return df['Main_Organism'].astype(str).map(lambda s: s.split()[0] if isinstance(s, str) and s.strip() else "")
        if 'Organism' in df.columns:
            return df['Organism'].astype(str).map(lambda s: s.split()[0] if isinstance(s, str) and s.strip() else "")
        return pd.Series([""] * len(df), index=df.index)

    out_df['Species_norm'] = infer_species(out_df)
    out_df['Genus_norm'] = infer_genus(out_df)
    if 'Family' not in out_df.columns and 'family' in out_df.columns:
        out_df['Family'] = out_df['family']
    if 'Family' not in out_df.columns:
        out_df['Family'] = ""

    # IPSC status
    out_df['IPSC_status'] = out_df['frames_summary'].apply(lambda js: detect_ipsc_from_frames(js) if pd.notna(js) else "Unknown")

    # compute seq length in bp (prefer best_orf_len*3 else len(clean_seq))
    out_df['best_orf_len'] = pd.to_numeric(out_df.get('best_orf_len', 0), errors='coerce').fillna(0).astype(int)
    out_df['seq_len_bp'] = out_df.apply(
        lambda r: int(r['best_orf_len'] * 3) if r['best_orf_len'] > 0 else len(str(r.get('clean_seq',''))),
        axis=1
    )

    # Save per-sequence report & kept set
    base = os.path.splitext(os.path.basename(args.input))[0]
    report = os.path.join(args.output_dir, f"{base}_numt_report.xlsx")
    kept = os.path.join(args.output_dir, f"{base}_numt_filtered.xlsx")
    out_df.to_excel(report, index=False, engine="openpyxl")
    out_df.loc[~out_df["Final_remove_flag"]].to_excel(kept, index=False, engine="openpyxl")

    # Prepare NUMT-only dataframe
    numt_df = out_df[out_df['numt_flag'] == True].copy()

    # -------------------------
    # 1) Per-length exact counts by IPSC status
    # -------------------------
    if numt_df.empty:
        per_length = pd.DataFrame(columns=['seq_len_bp','IPSC_status','numt_count'])
        pivot_per_length = pd.DataFrame()
    else:
        per_length = numt_df.groupby(['seq_len_bp','IPSC_status']).size().reset_index(name='numt_count').sort_values(['seq_len_bp','IPSC_status'])
        pivot_per_length = per_length.pivot_table(index='seq_len_bp', columns='IPSC_status', values='numt_count', fill_value=0).reset_index()
        # ensure order
        for col in ['Present','Absent','Unknown']:
            if col not in pivot_per_length.columns:
                pivot_per_length[col] = 0
        pivot_per_length = pivot_per_length[['seq_len_bp','Present','Absent','Unknown']]

    per_length_csv = os.path.join(args.output_dir, f"{base}_lengthcounts_per_length.csv")
    pivot_csv = os.path.join(args.output_dir, f"{base}_lengthcounts_pivot_per_length.csv")
    per_length.to_csv(per_length_csv, index=False)
    pivot_per_length.to_csv(pivot_csv, index=False)

    # -------------------------
    # 2) Binned counts by IPSC status
    # -------------------------
    size_bins = DEFAULT_SIZE_BINS
    if args.size_bins_json:
        try:
            with open(args.size_bins_json, 'r') as fh:
                jd = json.load(fh)
                size_bins = [(int(x[0]), int(x[1]), x[2]) for x in jd]
        except Exception:
            pass

    # if bin_width specified, create uniform bins
    if args.bin_width and args.bin_width > 0 and not numt_df.empty:
        minlen = int(numt_df['seq_len_bp'].min())
        maxlen = int(numt_df['seq_len_bp'].max())
        bins = []
        start = max(0, (minlen // args.bin_width) * args.bin_width)
        while start <= maxlen:
            end = start + args.bin_width - 1
            label = f"{start+1}-{end+1}"
            bins.append((start+1, end+1, label))
            start += args.bin_width
        if bins:
            size_bins = bins

    binned_records = []
    for low, high, label in size_bins:
        for ips in ['Present','Absent','Unknown']:
            if numt_df.empty:
                cnt = 0
                mean_len = np.nan
            else:
                sel = numt_df[(numt_df['seq_len_bp'] >= low) & (numt_df['seq_len_bp'] <= high) & (numt_df['IPSC_status']==ips)]
                cnt = len(sel)
                mean_len = int(sel['seq_len_bp'].mean()) if cnt>0 else np.nan
            binned_records.append({
                'Size_bin_label': label,
                'Size_low_bp': low,
                'Size_high_bp': high,
                'IPSC_status': ips,
                'numt_count': cnt,
                'mean_length_bp': mean_len
            })
    binned_df = pd.DataFrame(binned_records)

    # pivot binned to wide
    if not binned_df.empty:
        binned_pivot = binned_df.pivot_table(index=['Size_bin_label','Size_low_bp','Size_high_bp'], columns='IPSC_status', values='numt_count', fill_value=0).reset_index()
        for col in ['Present','Absent','Unknown']:
            if col not in binned_pivot.columns:
                binned_pivot[col] = 0
        binned_pivot = binned_pivot[['Size_bin_label','Size_low_bp','Size_high_bp','Present','Absent','Unknown']]
    else:
        binned_pivot = pd.DataFrame()

    binned_csv = os.path.join(args.output_dir, f"{base}_lengthcounts_binned.csv")
    binned_pivot_csv = os.path.join(args.output_dir, f"{base}_lengthcounts_binned_pivot.csv")
    binned_df.to_csv(binned_csv, index=False)
    binned_pivot.to_csv(binned_pivot_csv, index=False)

    # -------------------------
    # 3) Save everything to additional_summaries.xlsx
    # -------------------------
    summary_xlsx = os.path.join(args.output_dir, f"{base}_additional_summaries.xlsx")
    with pd.ExcelWriter(summary_xlsx, engine="openpyxl") as writer:
        per_length.to_excel(writer, sheet_name='LengthDist_per_length', index=False)
        pivot_per_length.to_excel(writer, sheet_name='LengthDist_per_length_pivot', index=False)
        binned_df.to_excel(writer, sheet_name='LengthDist_binned_long', index=False)
        binned_pivot.to_excel(writer, sheet_name='LengthDist_binned_pivot', index=False)
        out_df.to_excel(writer, sheet_name='per_sequence_full', index=False)
        # quick summary
        quick = {
            'Total_sequences': [len(out_df)],
            'Total_NUMTs_flagged': [int(out_df['numt_flag'].sum())],
            'GC_outliers_flagged': [int(out_df['GC_outlier'].sum())],
            'Final_removed_total': [int(out_df['Final_remove_flag'].sum())]
        }
        pd.DataFrame(quick).to_excel(writer, sheet_name='quick_summary', index=False)

    # console summary
    print("\n=== DONE: plotting-ready tables produced ===")
    print(f"Per-length CSV:       {per_length_csv}")
    print(f"Per-length pivot CSV: {pivot_csv}")
    print(f"Binned CSV:           {binned_csv}")
    print(f"Binned pivot CSV:     {binned_pivot_csv}")
    print(f"Additional summaries: {summary_xlsx}")
    print(f"Per-sequence report:  {report}")
    print(f"Kept sequences:       {kept}")

if __name__ == "__main__":
    main()
