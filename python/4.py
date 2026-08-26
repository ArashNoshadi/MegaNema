import pandas as pd
import subprocess
import os
import tempfile
import re
import platform
from Bio import AlignIO
from Bio.Align import PairwiseAligner
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
# ===========================
# CONFIGURATION & INPUTS
# ===========================
# --- SHUTDOWN CONFIGURATION ---
# اگر می‌خواهید کامپیوتر بعد از اتمام کار خاموش شود، این مقدار را True کنید
SHUTDOWN_AFTER_DONE = False  # Set to True to shut down PC after completion


INPUT_FILES_LIST = [
   # r"C:\Users\pc\Desktop\clear suspect\ITS1.xlsx",
   # r"C:\Users\pc\Desktop\clear suspect\ITS2.xlsx",
   # r"C:\Users\pc\Desktop\clear suspect\18S.xlsx",
  #  r"C:\Users\pc\Desktop\clear suspect\28S.xlsx",
    r"C:\Users\pc\Desktop\clear suspect\5.8S.xlsx",
   # r"C:\Users\pc\Desktop\clear suspect\COX1.xlsx",


]

MAFFT_BAT = r"C:\Users\pc\Desktop\mafft-win\mafft.bat"
OUTPUT_DIR = r"C:\Users\pc\Desktop\clear suspect\phylogen by counsensus"
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
# ===========================
# HELPER FUNCTIONS
# ===========================
def perform_shutdown():
    """
    Executes system shutdown command based on the operating system.
    """
    system_name = platform.system()
    print(f"\n[System] Initiating shutdown sequence for {system_name}...")
   
    try:
        if system_name == "Windows":
            # /s = shutdown, /t 5 = wait 5 seconds
            os.system("shutdown /s /t 5")
        elif system_name == "Linux" or system_name == "Darwin":
            # Linux or Mac
            os.system("shutdown -h now")
        else:
            print("[System] OS not recognized. Shutdown aborted.")
    except Exception as e:
        print(f"[Error] Could not shut down: {e}")
def clean_sequence(seq):
    """
    Cleans the sequence:
    1. Converts to Upper Case.
    2. Converts Uracil (U) to Thymine (T) for DNA compatibility.
    3. Replaces non-standard characters with 'N', but preserves IUPAC ambiguity codes and gaps (-).
    """
    if not isinstance(seq, str):
        return ""
   
    seq = seq.strip().upper()
    if not seq:
        return ""
   
    seq = seq.replace('U', 'T')
   
    # Preserve ACGTRYSWKMBDHVN-
    # Replace anything else with 'N'
    return re.sub(r'[^ACGTRYSWKMBDHVN-]', 'N', seq)
def dumb_consensus(alignment, threshold=0.5, ambiguous='N', gap='-'):
    """
    Improved manual consensus: accounts for gaps in denominator, handles ties by setting to ambiguous.
    - Denominator is total sequences (including gaps).
    - If gap frequency >= threshold, output gap.
    - Else, if max non-gap frequency >= threshold, output that base (if tie, ambiguous).
    - Else, ambiguous.
    """
    cons = ''
    aln_len = alignment.get_alignment_length()
    num_seqs = len(alignment)
    for i in range(aln_len):
        col = [rec.seq[i] for rec in alignment]
        counts = Counter(col)
        gap_count = counts.get(gap, 0)
        if gap_count / num_seqs >= threshold:
            cons += gap
            continue
        non_gap_counts = {k: v for k, v in counts.items() if k != gap}
        if not non_gap_counts:
            cons += ambiguous
            continue
        max_count = max(non_gap_counts.values())
        max_bases = [k for k, v in non_gap_counts.items() if v == max_count]
        if len(max_bases) > 1:
            cons += ambiguous  # Tie: ambiguous
        elif max_count / num_seqs >= threshold:
            cons += max_bases[0]
        else:
            cons += ambiguous
    return cons.replace(gap, 'N').upper()  # Replace gaps with N and ensure uppercase
def get_consensus_sequence(sequences, mafft_path):
    """
    Hybrid: For len==2, use PairwiseAligner + dumb_consensus; for >2, MAFFT + dumb_consensus.
    """
    # 1. Clean Sequences
    clean_seqs = [clean_sequence(s) for s in sequences]
    clean_seqs = [s for s in clean_seqs if s]
   
    if not clean_seqs:
        return ""
    unique_seqs = list(set(clean_seqs))
    if len(unique_seqs) == 1:
        return unique_seqs[0].upper()
   
    # For len==2: Pairwise + dumb_consensus
    if len(clean_seqs) == 2:
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aln = aligner.align(clean_seqs[0], clean_seqs[1])[0]
        # Convert to MultipleSeqAlignment for dumb_consensus
        msa = MultipleSeqAlignment([
            SeqRecord(Seq(aln[0]), id="seq0"),
            SeqRecord(Seq(aln[1]), id="seq1")
        ])
        return dumb_consensus(msa, threshold=0.5, ambiguous='N', gap='-')
   
    # For >2, use MAFFT
    tmp_in_path = tempfile.mktemp(suffix='.fasta')
    tmp_out_path = tempfile.mktemp(suffix='.aligned.fasta')
   
    try:
        # Write input file and close it
        with open(tmp_in_path, 'w') as tmp_in:
            for i, seq in enumerate(clean_seqs):
                tmp_in.write(f">seq_{i}\n{seq}\n")
       
        cmd = f'"{mafft_path}" --auto "{tmp_in_path}"'
       
        with open(tmp_out_path, 'w') as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True, shell=True)
       
        if result.returncode != 0:
            stderr_msg = result.stderr.strip() if result.stderr else "No stderr"
            print(f"[MAFFT ERROR] RC={result.returncode}, STDERR={stderr_msg}")
            return ""
       
        if not os.path.exists(tmp_out_path) or os.path.getsize(tmp_out_path) <= 0:
            print(f"[MAFFT OUTPUT ERROR] {tmp_out_path} missing or empty")
            return ""
       
        try:
            alignment = AlignIO.read(tmp_out_path, "fasta")
            if len(alignment) == 0:
                print("[ALIGNIO ERROR] Empty alignment")
                return ""
           
            cons = dumb_consensus(alignment, threshold=0.5, ambiguous='N', gap='-')
            cons = cons.strip()
            if not cons:
                return ""
            return cons.upper()
        except Exception as e:
            print(f"[ALIGNIO/CONSENSUS ERROR] {str(e)}")
            try:
                with open(tmp_out_path, 'r') as f:
                    preview = f.read(500)
                    print(f"[FILE PREVIEW] {repr(preview)[:300]}...")
            except Exception as preview_e:
                print(f"[PREVIEW ERROR] {preview_e}")
            return ""
       
    except Exception as e:
        print(f"[GENERAL ERROR in get_consensus] {str(e)}")
        return ""
    finally:
        for path in [tmp_in_path, tmp_out_path]:
            if os.path.exists(path):
                try:
                    os.remove(path)
                except OSError as e:
                    print(f"[TEMP CLEANUP ERROR] {path}: {e}")
def process_species_group(args):
    key, seq_list, mafft_path = args
    cons = get_consensus_sequence(seq_list, mafft_path)
    if cons:
        main_org, family, species = key
        return {
            'Main_Organism': main_org, 'Family': family, 'ID_Name': species,
            'Level': 'Species', 'Consensus_Sequence': cons
        }
    return None
def process_genus_group(args):
    main_org, family, seq_list, mafft_path = args
    cons = get_consensus_sequence(seq_list, mafft_path)
    if cons:
        return {
            'Main_Organism': main_org, 'Family': family, 'ID_Name': f"{main_org} sp.",
            'Level': 'Genus', 'Consensus_Sequence': cons
        }
    return None
def save_results(final_df, out_excel, out_fasta):
    """
    Save the current results to Excel and FASTA.
    Also clean all-N sequences before saving.
    Add Sequence_Length column.
    """
    if final_df.empty:
        return
    # Clean: remove rows where Consensus_Sequence is all N's
    mask_all_n = final_df['Consensus_Sequence'].str.match(r'^N+$', na=False)
    final_df = final_df[~mask_all_n]
    # Add length column
    final_df['Sequence_Length'] = final_df['Consensus_Sequence'].apply(len)
    final_df.to_excel(out_excel, index=False)
    
    with open(out_fasta, "w") as f:
        for _, row in final_df.iterrows():
            header = str(row['ID_Name']).replace(" ", "_")
            family = str(row['Family']).replace(" ", "_")
            f.write(f">{header}|{family}\n{row['Consensus_Sequence']}\n")

def analyze_and_save(df, base_name, mafft_path, output_dir, suffix_label):
    """
    Processes dataframe groupings (Species + Genus levels) and saves results.
    Uses multiprocessing with limited workers for parallel processing to speed up without freezing.
    Extracts seq_list before submit to reduce pickling overhead.
    Supports resume: if output exists, clean all-N sequences, skip processed groups, append new results.
    Live save: save after each new result.
    Re-process groups that were all-N and removed.
    """
    out_excel = os.path.join(output_dir, f"Consensus_{base_name}_{suffix_label}.xlsx")
    out_fasta = os.path.join(output_dir, f"Consensus_{base_name}_{suffix_label}.fasta")
    
    # Load existing if exists
    if os.path.exists(out_excel):
        try:
            existing_df = pd.read_excel(out_excel)
            # Identify all-N rows to re-process
            mask_all_n = existing_df['Consensus_Sequence'].str.match(r'^N+$', na=False)
            all_n_df = existing_df[mask_all_n]
            # Clean: remove all-N rows
            existing_df = existing_df[~mask_all_n]
            # Save cleaned back immediately
            save_results(existing_df, out_excel, out_fasta)
            # Processed species: tuples (Main_Organism, Family, ID_Name) for Species level (non all-N)
            processed_species = set(
                existing_df[existing_df['Level'] == 'Species'][['Main_Organism', 'Family', 'ID_Name']]
                .itertuples(index=False, name=None)
            )
            # Processed genus: set of Main_Organism for Genus level (non all-N)
            processed_genus = set(
                existing_df[existing_df['Level'] == 'Genus']['Main_Organism'].unique()
            )
            # Collect groups to re-process from all-N rows
            reprocess_species = set(
                all_n_df[all_n_df['Level'] == 'Species'][['Main_Organism', 'Family', 'ID_Name']]
                .itertuples(index=False, name=None)
            )
            reprocess_genus = set(
                all_n_df[all_n_df['Level'] == 'Genus']['Main_Organism'].unique()
            )
            tqdm.write(f"  > Resuming from existing file: {len(existing_df)} entries after cleaning. Re-processing {len(reprocess_species)} species and {len(reprocess_genus)} genus groups that were all-N.")
        except Exception as e:
            tqdm.write(f"  > [ERROR] Failed to load/clean existing file: {e}. Starting fresh.")
            processed_species = set()
            processed_genus = set()
            reprocess_species = set()
            reprocess_genus = set()
            existing_df = pd.DataFrame()
    else:
        processed_species = set()
        processed_genus = set()
        reprocess_species = set()
        reprocess_genus = set()
        existing_df = pd.DataFrame()
   
    final_df = existing_df.copy()
   
    # --- 1. Species Level Analysis ---
    try:
        df_sp = df.dropna(subset=['Main_Organism', 'Family', 'Species'])
        species_groups = {key: group for key, group in df_sp.groupby(['Main_Organism', 'Family', 'Species'])}
        # Unprocessed groups: not in processed_species or in reprocess_species
        unprocessed_species_groups = [
            (key, species_groups[key]) for key in species_groups 
            if key not in processed_species or key in reprocess_species
        ]
        total_species = len(species_groups)
        unprocessed_count = len(unprocessed_species_groups)
       
        desc_text = f"> {suffix_label:<15} | Species Level"
        if unprocessed_count > 0:
            with ProcessPoolExecutor(max_workers=8) as executor:
                futures = []
                for key, group in unprocessed_species_groups:
                    seq_list = group['RNA Sequence'].tolist()
                    futures.append(executor.submit(process_species_group, (key, seq_list, mafft_path)))
                for future in tqdm(as_completed(futures), total=unprocessed_count, desc=desc_text, unit="grp", leave=False):
                    result = future.result()
                    if result:
                        new_df = pd.DataFrame([result])
                        final_df = pd.concat([final_df, new_df], ignore_index=True)
                        final_df = final_df.drop_duplicates(subset=['Main_Organism', 'Family', 'ID_Name', 'Level'], keep='last')
                        save_results(final_df, out_excel, out_fasta)  # Live save
        num_species_cons = len(final_df[final_df['Level']=='Species']) - len(existing_df[existing_df['Level']=='Species'])
        tqdm.write(f"  > Species: {num_species_cons} new / {unprocessed_count} unprocessed (total groups: {total_species})")
    except KeyError as e:
        tqdm.write(f"  > [SKIP Species] Missing column: {e}")
    except Exception as e:
        tqdm.write(f"  > [ERROR] Species processing failed: {e}")
   
    # --- 2. Genus Level Analysis ---
    try:
        df_gen = df.dropna(subset=['Main_Organism'])
        genus_groups = {key: group for key, group in df_gen.groupby(['Main_Organism'])}
        # Unprocessed genus: not in processed_genus or in reprocess_genus
        unprocessed_genus_groups = []
        for main_org in genus_groups:
            if main_org in processed_genus and main_org not in reprocess_genus:
                continue
            group = genus_groups[main_org]
            mode_families = group['Family'].mode()
            if mode_families.empty:
                continue
            family = mode_families.iloc[0]
            unprocessed_genus_groups.append((main_org, family, group))
        total_genus = len(genus_groups)
        unprocessed_count = len(unprocessed_genus_groups)
       
        desc_text = f"> {suffix_label:<15} | Genus Level"
        if unprocessed_count > 0:
            with ProcessPoolExecutor(max_workers=4) as executor:
                futures = []
                for main_org, family, group in unprocessed_genus_groups:
                    seq_list = group['RNA Sequence'].tolist()
                    futures.append(executor.submit(process_genus_group, (main_org, family, seq_list, mafft_path)))
                for future in tqdm(as_completed(futures), total=unprocessed_count, desc=desc_text, unit="grp", leave=False):
                    result = future.result()
                    if result:
                        new_df = pd.DataFrame([result])
                        final_df = pd.concat([final_df, new_df], ignore_index=True)
                        final_df = final_df.drop_duplicates(subset=['Main_Organism', 'Family', 'ID_Name', 'Level'], keep='last')
                        save_results(final_df, out_excel, out_fasta)  # Live save
        num_genus_cons = len(final_df[final_df['Level']=='Genus']) - len(existing_df[existing_df['Level']=='Genus'])
        tqdm.write(f"  > Genus: {num_genus_cons} new / {unprocessed_count} unprocessed (total groups: {total_genus})")
    except KeyError as e:
        tqdm.write(f"  > [SKIP Genus] Missing column: {e}")
    except Exception as e:
        tqdm.write(f"  > [ERROR] Genus processing failed: {e}")
   
    # --- 3. Final Save (though live saves happened) ---
    if not final_df.empty:
        tqdm.write(f"  > Final saved: {out_excel} ({len(final_df)} entries)")
    else:
        tqdm.write(f"  > [WARNING] No results to save for {suffix_label}.")
def process_single_file(file_path):
    """Processes a single file: creates dedicated subfolder, runs ALL_DATA and FILTERED analysis."""
    file_name = os.path.basename(file_path)
    base_name = os.path.splitext(file_name)[0]
   
    subdir = os.path.join(OUTPUT_DIR, base_name)
    os.makedirs(subdir, exist_ok=True)
   
    tqdm.write(f"\n>>> File: {file_name}")
    tqdm.write(f"    -> Output folder: {subdir}")
   
    if not os.path.exists(file_path):
        tqdm.write(f"    !!! File not found: {file_path}")
        return
    try:
        df = pd.read_excel(file_path)
        df.columns = [c.strip() for c in df.columns]
        tqdm.write(f"    Loaded: {len(df)} rows")
    except Exception as e:
        tqdm.write(f"    !!! Error reading Excel: {e}")
        return
   
    # RUN 1: ALL DATA
    analyze_and_save(df, base_name, MAFFT_BAT, subdir, "ALL_DATA")
   
    # RUN 2: FILTERED DATA (if column exists)
    if 'RNA Length Check' in df.columns:
        try:
            rna_check = pd.to_numeric(df['RNA Length Check'], errors='coerce')
            df_filtered = df[rna_check == 1].copy()
            if len(df_filtered) > 0:
                analyze_and_save(df_filtered, base_name, MAFFT_BAT, subdir, "FILTERED_CHECK")
            else:
                tqdm.write("    > [Skip] Filtered dataset empty")
        except Exception as e:
            tqdm.write(f"    > [Error] Filtering failed: {e}")
    else:
        tqdm.write("    > [Skip] 'RNA Length Check' column missing")
# ===========================
# MAIN
# ===========================
if __name__ == "__main__":
    print(f"--- Batch Processing Started ---")
    print(f"Output base: {OUTPUT_DIR}")
    print(f"MAFFT: {MAFFT_BAT}")
   
    if not os.path.exists(MAFFT_BAT):
        print(f"CRITICAL ERROR: MAFFT executable not found at: {MAFFT_BAT}")
        exit(1)
   
    with tqdm(INPUT_FILES_LIST, desc="Overall Progress", unit="file") as pbar:
        for fpath in pbar:
            pbar.set_postfix({"File": os.path.basename(fpath)})
            process_single_file(fpath)
           
    print("\n--- All files processed ---")
   
    # CHECK FOR SHUTDOWN
    if SHUTDOWN_AFTER_DONE:
        print("Shutdown requested. Powering off in 5 seconds...")
        perform_shutdown()