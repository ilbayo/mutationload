# app/mutation_load.py

import io
import re
import pandas as pd

# flexible header groups
_COL_GROUPS = {
    "Chromosome": ["chromosome", "chrom", "chr", "#chrom", "#chromosome"],
    "Position": ["position", "pos", "start", "bp"],
    "AltAlleleFreq": [
        "altallelefreq",
        "alt_allele_freq",
        "alt_freq",
        "allele_freq",
        "af",
        "altfreq",
    ],
}


def _pick_col(df: pd.DataFrame, wanted: list[str], label: str) -> str:
    """Find a column in df that matches any of wanted (case-insensitive)."""
    lower_map = {c.lower(): c for c in df.columns}
    for w in wanted:
        lw = w.lower()
        if lw in lower_map:
            return lower_map[lw]
    raise ValueError(f"Missing required column: {label}")


def load_variant_file(raw) -> pd.DataFrame:
    """
    Accepts: bytes, str, file-like.
    Returns: dataframe with columns Chromosome, Position, AltAlleleFreq
    """
    # 1) normalize to text
    if isinstance(raw, (bytes, bytearray)):
        text = raw.decode("utf-8", errors="ignore")
    elif hasattr(raw, "read"):
        # file-like
        text = raw.read()
        if isinstance(text, bytes):
            text = text.decode("utf-8", errors="ignore")
    elif isinstance(raw, str):
        text = raw
    else:
        raise ValueError(f"Unsupported upload type: {type(raw)}")

    text = text.strip()
    if not text:
        raise ValueError("Uploaded file is empty")

    # Sometimes users paste TSV but with multiple spaces — normalize
    # We'll first try auto-sep
    df = None
    try:
        df = pd.read_csv(io.StringIO(text), sep=None, engine="python", comment="#")
    except Exception:
        # fallback: collapse 2+ spaces to tab
        cleaned = re.sub(r"[ ]{2,}", "\t", text)
        df = pd.read_csv(io.StringIO(cleaned), sep="\t", comment="#")

    if df is None or df.empty:
        raise ValueError("No columns to parse from file")

    # 2) map columns
    chrom_col = _pick_col(df, _COL_GROUPS["Chromosome"], "Chromosome")
    pos_col = _pick_col(df, _COL_GROUPS["Position"], "Position")
    af_col = _pick_col(df, _COL_GROUPS["AltAlleleFreq"], "AltAlleleFreq")

    df = df.rename(
        columns={
            chrom_col: "Chromosome",
            pos_col: "Position",
            af_col: "AltAlleleFreq",
        }
    )

    # 3) coerce numeric
    df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
    df["AltAlleleFreq"] = pd.to_numeric(df["AltAlleleFreq"], errors="coerce")
    df = df.dropna(subset=["Position", "AltAlleleFreq"])

    if df.empty:
        raise ValueError("After cleaning, no valid rows remained")

    df["Position"] = df["Position"].astype(int)

    return df


def compute_mutation_load(
    df: pd.DataFrame,
    chrom: str,
    start: int,
    end: int,
    bin_size: int = 30,
):
    """
    Reproduce your sliding-window “event_count * avg_alt_freq” idea.
    Filters AltAlleleFreq <= 0.35 like your local code.
    """
    # filter by region and AF
    region = df[
        (df["Chromosome"] == chrom)
        & (df["Position"] >= start)
        & (df["Position"] <= end)
        & (df["AltAlleleFreq"] <= 0.35)
    ]

    xs, ys = [], []
    for s in range(start, end - bin_size + 2):
        e = s + bin_size - 1
        bin_events = region[(region["Position"] >= s) & (region["Position"] <= e)]
        if len(bin_events) == 0:
            continue
        event_count = len(bin_events)
        avg_af = bin_events["AltAlleleFreq"].mean()
        score = event_count * avg_af
        xs.append(s)
        ys.append(score)

    return xs, ys
