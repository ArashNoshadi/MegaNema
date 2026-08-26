import os
import ssl
import certifi
import socket
import time
import logging
import io
import requests
import sqlite3
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from Bio import Entrez, SeqIO
from Bio.Seq import UndefinedSequenceError
import re
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# ─── تنظیم SSL با Certifi ────────────────────────────────────────────────────
os.environ['SSL_CERT_FILE'] = certifi.where()

# ─── GLOBAL TIMEOUT FOR SOCKETS ───────────────────────────────────────────────
socket.setdefaulttimeout(5)

# ─── NCBI SETTINGS ─────────────────────────────────────────────────────────────
Entrez.email   = "@gmail.com"
Entrez.tool    = "n***"
Entrez.api_key = "////7609"
Entrez.base    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# ─── PATHS ────────────────────────────────────────────────────────────────────
INPUT_FILE     = Path(r"C:/Users/pc/Desktop/nema-itss.txt")
PARQUET_PATH   = Path(r"C:/Users/pc/Desktop/out-RNA/nema-itss.parquet")
DB_PATH        = Path(r"C:/Users/pc/Desktop/out-RNA/nema-itss-ids.db")
EXCEL_PATH     = Path(r"C:/Users/pc/Desktop/out-RNA/nema-itss.xlsx")

# ─── PARAMETERS ───────────────────────────────────────────────────────────────
BATCH_SIZE       = 200
BULK_FETCH_SIZE  = 50
MAX_WORKERS      = 3
EFETCH_RETRIES   = 3
DELAY            = 0.1
CHECK_INTERVAL   = 3
VACUUM_INTERVAL  = 5000

# ─── LOGGING ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger()

# ─── FIELDS & MARKERS ─────────────────────────────────────────────────────────
FIELDNAMES = [
    "ID","DEFINITION","marker","RNA Sequence","Length (RNA)",
    "Host","Geo Loc Name","Collection Date","Organism","Mol Type","Product"
]
marker_map = [
    # ── 5S ریبوزومی ─────────────────────────────────────────────────────────────
    (r"\b5[\s\-]?s[\s\-]?rRNA\b",                          "5S"),

    # ── External Transcribed Spacer (ETS) ──────────────────────────────────────
    (r"\bexternal[\s\-]?transcribed[\s\-]?spacer\b",       "ETS"),
    (r"\bets\b",                                           "ETS"),

    # ── Internal Transcribed Spacer (ITS) ──────────────────────────────────────
    (r"\bITS[\s\-]?region\b",                              "ITS"),
    (r"\bits[\s\-]?1\b",                                   "ITS1"),
    (r"\bits[\s\-]?2\b",                                   "ITS2"),
    (r"\binternal[\s\-]?transcribed[\s\-]?spacer[\s\-]?1\b","ITS1"),
    (r"\binternal[\s\-]?transcribed[\s\-]?spacer[\s\-]?2\b","ITS2"),

    # ── 18S (SSU) ────────────────────────────────────────────────────────────────
    (r"\b18[\s\-]?s\b",                                     "18S"),
    (r"\bssu[\s\-]?rrna\b",                                 "18S"),
    (r"\brRNA[\s\-]?small[\s\-]?subunit\b",                 "18S"),
    (r"\bsmall[\s\-]?subunit[\s\-]?ribosomal\b",            "18S"),
    (r"\b18[\s\-]?s[-\s]?rDNA\b",                           "18S"),

    # ── 28S (LSU) ───────────────────────────────────────────────────────────────
    (r"\b28[\s\-]?s\b",                                     "28S"),
    (r"\blsu[\s\-]?rrna\b",                                 "28S"),
    (r"\brRNA[\s\-]?large[\s\-]?subunit\b",                 "28S"),
    (r"\blarge[\s\-]?subunit[\s\-]?ribosomal\b",            "28S"),
    (r"\b28[\s\-]?s[-\s]?rDNA\b",                           "28S"),

    # ── 12S (mito SSU) ──────────────────────────────────────────────────────────
    (r"\b12[\s\-]?s\b",                                     "12S"),
    (r"\brrn[sS]\b",                                        "12S"),

    # ── 16S (mito LSU) ──────────────────────────────────────────────────────────
    (r"\b16[\s\-]?s\b",                                     "16S"),
    (r"\brrn[lL]\b",                                        "16S"),

    # ── 5.8S ─────────────────────────────────────────────────────────────────────
    (r"\b5[.\s\-]?8[\s\-]?s\b",                             "5.8S"),
    (r"\b5[.\s\-]?8[\s\-]?s[\s\-]?rRNA\b",                  "5.8S"),
    (r"\b5[.\s\-]?8[\s\-]?s[\s\-]?subunit\b",               "5.8S"),
    (r"\b5[.\s\-]?8[\s\-]?s[\s\-]?ribosomal\b",             "5.8S"),
    (r"\b5[.\s\-]?8[\s\-]?s[\s\-]?ribosomal[\s\-]?RNA\b",   "5.8S"),

    # ── 23S ─────────────────────────────────────────────────────────────────────
    (r"\b23[\s\-]?s\b",                                     "23S"),
    (r"\b23[\s\-]?s[\s\-]?rRNA\b",                          "23S"),

    # ── ژن‌های میتوکندریایی: COX1, COX2, COX3, CYTB ────────────────────────────
    #   COX1
    (r"\bcox[\s\-]?1\b",                                    "COX1"),
    (r"\bco[iI]\b",                                         "COX1"),
    (r"\bcytochrome[\s\-]?c[\s\-]?oxidase[\s\-]?subunit[\s]?1\b","COX1"),
    (r"\bcytochrome[\s\-]?oxidase[\s\-]?subunit\s?I\b",     "COX1"),
    #   COX2
    (r"\bcox[\s\-]?2\b",                                    "COX2"),
    (r"\bcytochrome[\s\-]?c[\s\-]?oxidase[\s\-]?subunit[\s]?2\b","COX2"),
    (r"\bcytochrome[\s\-]?oxidase[\s\-]?subunit\s?II\b",    "COX2"),
    #   COX3
    (r"\bcox[\s\-]?3\b",                                    "COX3"),
    (r"\bcytochrome[\s\-]?c[\s\-]?oxidase[\s\-]?subunit[\s]?3\b","COX3"),
    (r"\bcytochrome[\s\-]?oxidase[\s\-]?subunit\s?III\b",   "COX3"),
    #   CYTB
    (r"\bcytb\b",                                           "CYTB"),
    (r"\bcytochrome[\s\-]?b\b",                             "CYTB"),

    # ── NADH Dehydrogenase (ND1–ND6) ────────────────────────────────────────────
    (r"\bnd([1-6])\b",                                      lambda m: f"ND{m.group(1)}"),
    (r"\bnad[hH]?[\s\-]?dehydrogenase[\s\-]?subunit[\s]?([1-6])\b",
                                                           lambda m: f"ND{m.group(1)}"),

    # ── ژن‌های ATP ───────────────────────────────────────────────────────────────
    (r"\batp[fF][\s\-]?atp[hH]\b",                          "atpF-atpH"),
    (r"\batp[\s\-]?ase[\s\-]?subunit\s?6\b",                "ATP6"),
    (r"\batp6\b",                                           "ATP6"),
    (r"\batp[\s\-]?synthetase\s?subunit\s?alpha\b",         "ATPα"),

    # ── گیاهی (rbcL, matK) ───────────────────────────────────────────────────────
    (r"\brbc[\s\-]?l\b",                                    "rbcL"),
    (r"\bmat[\s\-]?k\b",                                    "matK"),

    # ── ترو-psb و psb-k/psb-i ────────────────────────────────────────────────────
    (r"\btrn[hH][\s\-]?psb[aA]\b",                          "trnH-psbA"),
    (r"\bpsb[kK][\s\-]?psb[iI]\b",                          "psbK-psbI"),

    # ── متغیرهای منطقه‌ای 18S/28S ───────────────────────────────────────────────
    (r"\bssu[-\s]?v4\b",                                    "18S-V4"),
    (r"\bssu[-\s]?v9\b",                                    "18S-V9"),
    (r"\blsu[-\s]?d2\b",                                    "28S-D2"),
    (r"\b28[\s\-]?s[-\s]?d2\b",                             "28S-D2"),

    # ── اضافه: تشخیص 'mitochondrion' و 'mitochondrial' ─────────────────────────
    (r"\bmitochondrion\b",                                  "mitochondrion"),
    (r"\bmitochondrial\b",                                  "mitochondrion"),
]


compiled_markers = [(re.compile(p, re.IGNORECASE), v) for p, v in marker_map]

# ─── HELPERS ──────────────────────────────────────────────────────────────────
def ensure_dir(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)

def load_ids(path: Path):
    xml = path.read_text()
    logger.info("Loaded input XML for IDs.")
    return re.findall(r"<Id>(\d+)</Id>", xml)

def check_internet(host="8.8.8.8", port=53, timeout=3):
    try:
        socket.setdefaulttimeout(timeout)
        socket.socket().connect((host, port))
        return True
    except OSError:
        return False

def parse_record(rec):
    desc = rec.description or ""
    found = []
    for regex, canon in compiled_markers:
        for m in regex.finditer(desc):
            val = canon(m) if callable(canon) else canon
            found.append(val)
    # حذف تکراری‌ها و حفظ ترتیب
    seen = set()
    markers = [x for x in found if not (x in seen or seen.add(x))]
    marker_str = "|".join(markers) if markers else "N/A"

    try:
        seq_str = str(rec.seq)
        length = len(rec.seq)
    except (UndefinedSequenceError, TypeError):
        seq_str = "N/A"
        length = 0

    entry = dict.fromkeys(FIELDNAMES, "N/A")
    entry.update({
        "ID": rec.id,
        "DEFINITION": desc,
        "marker": marker_str,
        "RNA Sequence": seq_str,
        "Length (RNA)": length
    })

    for feat in rec.features:
        q = feat.qualifiers
        if feat.type == "source":
            entry.update({
                "Host":           q.get("host", ["N/A"])[0],
                "Geo Loc Name":   q.get("geo_loc_name", ["N/A"])[0],
                "Collection Date":q.get("collection_date", ["N/A"])[0],
                "Organism":       q.get("organism", ["N/A"])[0],
                "Mol Type":       q.get("mol_type", ["N/A"])[0]
            })
        elif feat.type == "CDS":
            entry["Product"] = q.get("product", ["N/A"])[0]

    return entry

def fetch_bulk(id_list):
    time.sleep(DELAY)
    ids_str = ",".join(id_list)
    try:
        post = Entrez.epost(db="nucleotide", id=ids_str)
        data = Entrez.read(post); post.close()
        we, qk = data["WebEnv"], data["QueryKey"]
    except Exception as e:
        logger.error(f"EPost FAILED: {e}")
        return []

    recs = []
    for attempt in range(1, EFETCH_RETRIES+1):
        try:
            resp = requests.get(
                Entrez.base + "efetch.fcgi",
                params={"db":"nucleotide","WebEnv":we,"query_key":qk,
                        "rettype":"gb","retmode":"text"},
                timeout=10, verify=True
            )
            resp.raise_for_status()
            recs = list(SeqIO.parse(io.StringIO(resp.text), "gb"))
            break
        except Exception as e:
            logger.warning(f"EFetch attempt {attempt} failed: {e}")
            time.sleep(DELAY * (2**(attempt-1)))
    if not recs:
        logger.error("All EFetch retries failed")
        return []

    return [parse_record(r) for r in recs]

class ParquetAppender:
    def __init__(self, path):
        self.path = path
        self.writer = None
        self.schema = None

    def write(self, df: pd.DataFrame):
        table = pa.Table.from_pandas(df, schema=self.schema, preserve_index=False)
        if self.writer is None:
            self.schema = table.schema
            self.writer = pq.ParquetWriter(self.path, self.schema, compression='SNAPPY')
        self.writer.write_table(table)

    def close(self):
        if self.writer:
            self.writer.close()

def init_db():
    conn = sqlite3.connect(DB_PATH)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS processed (
            id TEXT PRIMARY KEY,
            ts TEXT
        );
    """)
    conn.commit()
    return conn

def log_ids(conn, records, counter):
    now = time.strftime("%Y-%m-%d %H:%M:%S")
    data = [(r["ID"], now) for r in records]
    conn.executemany("INSERT OR IGNORE INTO processed(id, ts) VALUES(?,?)", data)
    conn.commit()
    if counter % VACUUM_INTERVAL == 0:
        conn.execute("VACUUM;")
        conn.commit()

def main():
    logger.info("Pipeline STARTED")
    ensure_dir(PARQUET_PATH)

    ids = load_ids(INPUT_FILE)
    batches = [ids[i:i+BATCH_SIZE] for i in range(0, len(ids), BATCH_SIZE)]
    logger.info(f"{len(batches)} batches to process")

    conn = init_db()
    appender = ParquetAppender(PARQUET_PATH)
    total = 0
    log_counter = 0

    for bi, batch in enumerate(batches, 1):
        logger.info(f"Batch {bi}/{len(batches)} start")
        while not check_internet():
            time.sleep(CHECK_INTERVAL)

        recs_all = []
        chunks = [batch[i:i+BULK_FETCH_SIZE] for i in range(0, len(batch), BULK_FETCH_SIZE)]
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
            for fut in as_completed([pool.submit(fetch_bulk, ch) for ch in chunks]):
                recs_all.extend(fut.result() or [])

        if recs_all:
            df = pd.DataFrame.from_records(recs_all, columns=FIELDNAMES)
            appender.write(df)
            total += len(recs_all)

            log_ids(conn, recs_all, log_counter)
            log_counter += len(recs_all)

            logger.info(f"Wrote {len(recs_all)} rows (total {total})")
        else:
            logger.warning("Empty batch")

        if bi % 20 == 0:
            logger.info("Pausing 5 minutes")
            time.sleep(50)

    appender.close()
    conn.close()

    logger.info("Converting Parquet to Excel")
    df_all = pd.read_parquet(PARQUET_PATH, engine="pyarrow")
    df_all.to_excel(EXCEL_PATH, index=False, engine="xlsxwriter")
    logger.info("Excel saved to %s", EXCEL_PATH)
    logger.info("Pipeline FINISHED")

if __name__ == '__main__':
    main()  


