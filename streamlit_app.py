# app.py — DBsimilarity (Streamlit port)
# ------------------------------------------------
# Author: Ricardo M. Borges & LAABio-IPPN-UFRJ
# Date: 2025-10-24
#
# What you get (all original computations preserved, modular, multi-file):
# 1) Data intake (multi-upload): target structures (CSV), annotated table(s), metadata (CSV/TXT ; or TAB; comma ok)
# 2) Cleaning: drop empty/duplicate SMILES; robust column mapping; InChI/InChIKey/Formula/Exact masses + [M±H], [M+Na], [M+K], [M+NH4]
# 3) Merge target × annotated with tolerant key preprocessing (ignores _0/0_)
# 4) MassQL builders: MS1 + MS2 (with adducts, ppm, intensity%, optional fragment cardinality)
# 5) Fingerprints + similarity matrix → 3-column edges CSV (threshold) + network (Plotly) + isolated nodes CSV
# 6) Mordred descriptors (2D/3D) → PCA (Plotly), dendrogram (Matplotlib), Agglomerative clusters CSV, t-SNE
# 7) Exports: CSVs + HTML (network and t-SNE). Buttons for every artifact.
# 8) Performance knobs: sample limits, descriptor mode, similarity metric (Dice/Tanimoto), tsne params, etc.
# ------------------------------------------------

from __future__ import annotations

# Standard lib (trim this list to what you actually use)
import io
import os
from pathlib import Path
from typing import List, Dict, Tuple, Optional

import re
NUM_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")

import tempfile
from io import BytesIO

from rdkit import DataStructs
from rdkit.Chem import PandasTools
from rdkit.Chem.rdMolDescriptors import CalcMolFormula


import numpy as np
import pandas as pd
import streamlit as st

# --- Chemistry deps (guarded) ---
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    from rdkit.Chem.Draw import MolToImage
except Exception as e:
    st.stop()  # hard stop if RDKit isn't present




# Mordred depends on RDKit; import it only after RDKit succeeded
try:
    from mordred import Calculator, descriptors
except Exception:
    st.error(
        "The 'mordred' package isn’t installed. Add it to your conda/pip deps "
        "(e.g., environment.yml -> pip: [mordred])."
    )
    st.stop()

# --- Viz & analysis ---
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import networkx as nx
import matplotlib.pyplot as plt
from scipy.stats import chi2

# ML / stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.cluster import AgglomerativeClustering

try:
    import mpld3
    HAS_MPLD3 = True
except Exception:
    HAS_MPLD3 = False
# ...
if HAS_MPLD3:
    html_str = mpld3.fig_to_html(fig)
    st.download_button("⬇️ Download dendrogram (HTML, interactive)", data=html_str,
                       file_name="dendrogram.html", mime="text/html")



st.set_page_config(page_title="DBsimilarity — Streamlit", layout="wide")
st.title("DBsimilarity — Structure Mining, MassQL & Similarity (Streamlit)")

# --- persistent state for merged table ---
if "merged_df" not in st.session_state:
    st.session_state["merged_df"] = pd.DataFrame()

# -------- Arrow-free table renderer (no pyarrow dependency) --------
def show_df(df: pd.DataFrame, caption: str | None = None, rows: int | None = None):
    if caption:
        st.caption(caption)
    if rows is not None:
        df = df.head(rows)
    try:
        # avoid st.dataframe to keep pyarrow optional
        st.markdown(df.to_html(index=False), unsafe_allow_html=True)
    except Exception:
        st.write(df)

# Persistent downloads container
if "exports" not in st.session_state:
    st.session_state["exports"] = {}
exports: Dict[str, bytes] = st.session_state["exports"]


# ---------- Helpers ----------
H_MASS   = 1.007276466
NA_MASS  = 22.989218
K_MASS   = 38.963157
NH4_MASS = 18.033823

# (multiplier_of_M, adduct_mass, charge)
ADDUCT_SPECS = {
    "[M+H]+":     (1,  +H_MASS,  +1),
    "[M+Na]+":    (1,  +NA_MASS, +1),
    "[M+K]+":     (1,  +K_MASS,  +1),
    "[M+NH4]+":   (1,  +NH4_MASS,+1),
    "[M-H]-":     (1,  -H_MASS,  -1),
    # NEW:
    "[2M+H]+":    (2,  +H_MASS,  +1),   # dimer protonated
    "[M+2H]+2":   (1,  +2*H_MASS, +2),  # doubly charged
}


NUM_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")

def _to_float(val):
    if val is None:
        return None
    s = str(val).strip().replace(",", ".")
    m = NUM_RE.search(s)
    return float(m.group(0)) if m else None

@st.cache_data(show_spinner=False)
def smart_read_table(file: bytes, filename: str) -> pd.DataFrame:
    """Robust reader: tries ; then TAB, then comma."""
    name = filename.lower()
    content = io.BytesIO(file)
    encodings = ["utf-8", "ISO-8859-1", "latin1"]
    seps = [";", "\t", ","]
    for enc in encodings:
        for sep in seps:
            try:
                content.seek(0)
                df = pd.read_csv(content, sep=sep, encoding=enc)
                # Heuristic: if only 1 column and sep wrong, keep trying
                if df.shape[1] == 1 and sep != ",":
                    continue
                return df
            except Exception:
                continue
    # Fallback pandas auto
    content.seek(0)
    return pd.read_csv(content, engine="python")

@st.cache_data(show_spinner=False)
def read_sdf_to_df(file: bytes, filename: str) -> pd.DataFrame:
    tmp = Path(tempfile.gettempdir()) / filename
    with open(tmp, "wb") as f:
        f.write(file)
    df = PandasTools.LoadSDF(str(tmp), embedProps=True, molColName="ROMol")
    if "SMILES" not in df.columns:
        df["SMILES"] = df["ROMol"].apply(lambda m: Chem.MolToSmiles(m) if m is not None else None)
    return df.drop(columns=["ROMol"], errors="ignore")


# Column resolver (case-insensitive)
def resolve_col(df: pd.DataFrame, wanted: str) -> str:
    if wanted in df.columns: return wanted
    lowmap = {c.lower(): c for c in df.columns}
    key = wanted.lower()
    if key in lowmap: return lowmap[key]
    raise KeyError(f"Column '{wanted}' not found. Available: {', '.join(map(str, df.columns))}")

# Molecule utilities
def mol_from_smiles(smiles: str):
    try:
        m = Chem.MolFromSmiles(smiles)
        return m
    except Exception:
        return None

# Add calculated chem columns
def enrich_chem(df: pd.DataFrame, smiles_col: str) -> pd.DataFrame:
    df = df.copy()
    mols = [mol_from_smiles(s) for s in df[smiles_col].astype(str)]

    def _safe_inchi(m):
        try:
            return Chem.MolToInchi(m) if m else None
        except Exception:
            return None

    def _safe_inchikey(m):
        try:
            return Chem.MolToInchiKey(m) if m else None
        except Exception:
            return None

    df["MolFormula"] = [CalcMolFormula(m) if m else None for m in mols]
    df["Inchi"]      = [_safe_inchi(m) for m in mols]
    df["InchiKey"]   = [_safe_inchikey(m) for m in mols]
    exact = [Descriptors.ExactMolWt(m) if m else None for m in mols]
    df["ExactMass"]  = exact
    for ad, (mult, add_mass, z) in ADDUCT_SPECS.items():
        denom = abs(z) if z != 0 else 1
        df[f"mz {ad}"] = [(mult * mw + add_mass) / denom if mw is not None else None for mw in exact]
    return df



# MassQL builders (MS1 / MS2)
def _fmt(x: float, decimals: int=5) -> str: return f"{x:.{decimals}f}"

def compute_adduct_mzs(neutral_mass: float, adducts: List[str], decimals: int = 5) -> List[str]:
    mzs = []
    for a in adducts:
        mult, add_mass, z = ADDUCT_SPECS[a]
        denom = abs(z) if z != 0 else 1
        mz = (mult * neutral_mass + add_mass) / denom
        mzs.append(_fmt(mz, decimals))
    return mzs

@st.cache_data(show_spinner=False)
def build_ms1_queries(
    df_meta: pd.DataFrame, name_col: str, mass_col: str,
    adducts: List[str], ppm: int, intensity: int, decimals: int,
    separate: bool
) -> Dict:
    name_c = resolve_col(df_meta, name_col)
    mass_c = resolve_col(df_meta, mass_col)
    out = {}

    for _, row in df_meta.iterrows():
        name = str(row[name_c]).strip()
        mono = _to_float(row[mass_c])
        if mono is None:
            continue

        mzs = compute_adduct_mzs(mono, adducts, decimals)

        if not separate:
            mz_list = " OR ".join(mzs)
            q = (
                f"# {name}\n"
                "QUERY scaninfo(MS2DATA) WHERE\n"
                f"MS2PREC=(\n\t{mz_list}\n):TOLERANCEPPM={ppm}:INTENSITYPERCENT={intensity}"
            )
            out[name] = q
        else:
            per_adduct = {}
            for a, mz in zip(adducts, mzs):
                q = (
                    f"# {name} {a}\n"
                    "QUERY scaninfo(MS2DATA) WHERE\n"
                    f"MS2PREC=(\n\t{mz}\n):TOLERANCEPPM={ppm}:INTENSITYPERCENT={intensity}"
                )
                per_adduct[a] = q
            out[name] = per_adduct

    return out


@st.cache_data(show_spinner=False)
def build_ms2_queries(df_meta: pd.DataFrame, name_col: str, mass_col: str, frag_col: str,
                      adducts: List[str], ppm_prec:int, ppm_prod:int, intensity:int,
                      decimals:int, cmin:int, cmax:int, clamp:bool, ignore_zero:bool) -> Dict:
    name_c = resolve_col(df_meta, name_col)
    mass_c = resolve_col(df_meta, mass_col)
    frag_c = resolve_col(df_meta, frag_col)
    out = {}
    for _, row in df_meta.iterrows():
        name = str(row[name_c]).strip()
        mono = _to_float(row[mass_c])
        if mono is None:
            continue
        prec_list = " OR ".join(compute_adduct_mzs(mono, adducts, decimals))
        # parse fragments
        raw = row.get(frag_c, None)
        frags = []
        if pd.notna(raw) and str(raw).strip() != "":
            tokens = re.split(r"[;,\s]+", str(raw))
            seen = set()
            for t in tokens:
                v = _to_float(t)
                if v is None: continue
                if ignore_zero and v <= 0: continue
                vv = float(_fmt(v, decimals))
                if vv not in seen:
                    frags.append(vv); seen.add(vv)
        if frags:
            n = len(frags)
            cmi = min(cmin, n) if clamp else cmin
            cma = min(cmax, n) if clamp else cmax
            if cmi < 1: cmi = 1
            if cma < cmi: cma = cmi
            prod_list = " OR ".join(_fmt(v, decimals) for v in frags)
            q = (
                f"# {name}\n"
                "QUERY scaninfo(MS2DATA) WHERE\n"
                f"MS2PREC=({prec_list}):TOLERANCEPPM={ppm_prec} AND\n"
                f"MS2PROD=({prod_list}):CARDINALITY=range(min={cmi},max={cma}):"
                f"TOLERANCEPPM={ppm_prod}:INTENSITYPERCENT={intensity}"
            )
        else:
            q = (
                f"# {name}\n"
                "QUERY scaninfo(MS2DATA) WHERE\n"
                f"MS2PREC=({prec_list}):TOLERANCEPPM={ppm_prec}"
            )
        out[name] = q
    return out

# Merge utilities
def preprocess_key_column(df: pd.DataFrame, key: str) -> pd.DataFrame:
    d = df.copy()
    if key in d.columns:
        d[key] = d[key].astype(str).str.replace(r"_0|0_", "", regex=True)
    return d

def merge_by_key(dfs: List[pd.DataFrame], key: str) -> pd.DataFrame:
    dfs2 = [preprocess_key_column(x, key) for x in dfs]
    merged = pd.concat(dfs2, axis=0, ignore_index=True)
    merged = merged.groupby(key, as_index=False).first()
    for i, df in enumerate(dfs2, start=1):
        merged[f"in_df{i}"] = merged[key].isin(df[key]).astype(int)
    return merged.fillna(0)

# Similarity graph
def fingerprints(smiles: List[str], radius:int=2, nBits:int=2048):
    mols = [mol_from_smiles(s) for s in smiles]
    return [AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=nBits) if m else None for m in mols]

@st.cache_data(show_spinner=False)
def pairwise_similarity(ikeys: List[str], fps: List, metric: str = "dice") -> pd.DataFrame:
    n = len(fps)
    rows = []
    for i in range(n):
        if fps[i] is None: continue
        # Bulk comparisons for speed
        sims = [ (DataStructs.DiceSimilarity(fps[i], fps[j]) if metric=="dice" else DataStructs.TanimotoSimilarity(fps[i], fps[j])) if fps[j] is not None else 0.0 for j in range(i+1, n)]
        for j, s in enumerate(sims, start=i+1):
            rows.append((ikeys[i], ikeys[j], float(s)))
    return pd.DataFrame(rows, columns=["SOURCE","TARGET","CORRELATION"])

def _set_state_df_if_changed(key: str, new_df: pd.DataFrame):
    """Update session_state[key] only if it is missing or differs by content."""
    cur = st.session_state.get(key, None)
    if cur is None or not isinstance(cur, pd.DataFrame) or not cur.equals(new_df):
        st.session_state[key] = new_df


from PIL import Image, UnidentifiedImageError

STATIC_DIR = Path(__file__).parent / "static"
LOGO_PATH = STATIC_DIR / "LAABio.png"
try:
    logo = Image.open(LOGO_PATH)  # raises if missing
    st.sidebar.image(logo, use_container_width=True)
except FileNotFoundError:
    st.sidebar.warning("Logo not found at static/LAABio.png")

LOGO_DBsimilarity_PATH = STATIC_DIR / "DBsimilarity.png"
try:
    logo = Image.open(LOGO_DBsimilarity_PATH)  # raises if missing
    st.sidebar.image(logo, use_container_width=True)
except FileNotFoundError:
    st.sidebar.warning("Logo not found at static/DBsimilarity.png")

st.markdown(
    """
    Developed by **Ricardo M Borges** and **LAABio-IPPN-UFRJ**  
    contact: ricardo_mborges@yahoo.com.br  

    🔗 Details: [GitHub repository](https://github.com/RicardoMBorges/DBsimilarity_st)

    Check also: [DAFdiscovery](https://dafdiscovery.streamlit.app/)
    
    Check also: [TLC2Chrom](https://tlc2chrom.streamlit.app/)
    """
)

# PayPal donate button
st.sidebar.markdown("""
<hr>
<center>
<p>To support the app development:</p>
<a href="https://www.paypal.com/donate/?business=2FYTFNDV4F2D4&no_recurring=0&item_name=Support+with+%245+→+Send+receipt+to+tlc2chrom.app@gmail.com+with+your+login+email+→+Access+within+24h!&currency_code=USD" target="_blank">
    <img src="https://www.paypalobjects.com/en_US/i/btn/btn_donate_SM.gif" alt="Donate with PayPal button" border="0">
</a>
</center>
""", unsafe_allow_html=True)

st.sidebar.markdown("""---""")

TUTORIAL_URL = "https://github.com/RicardoMBorges/DBsimilarity_st"
try:
    st.sidebar.link_button("📘 Tutorial", TUTORIAL_URL)
except Exception:
    st.sidebar.markdown(
        f'<a href="{TUTORIAL_URL}" target="_blank">'
        '<button style="padding:0.6rem 1rem; border-radius:8px; border:1px solid #ddd; cursor:pointer;">📘 Tutorial</button>'
        '</a>',
        unsafe_allow_html=True,
    )


MockData_URL = "https://github.com/RicardoMBorges/DBsimilarity_st"
try:
    st.sidebar.link_button("Mock Data", MockData_URL)
except Exception:
    st.sidebar.markdown(
        f'<a href="{MockData_URL}" target="_blank">'
        '<button style="padding:0.6rem 1rem; border-radius:8px; border:1px solid #ddd; cursor:pointer;">Mock Data</button>'
        '</a>',
        unsafe_allow_html=True,
    )
    
VIDEO_URL = "https://github.com/RicardoMBorges/new_figure_metaboAnalyst"
try:
    st.sidebar.link_button("Video", VIDEO_URL)
except Exception:
    st.sidebar.markdown(
        f'<a href="{VIDEO_URL}" target="_blank">'
        '<button style="padding:0.6rem 1rem; border-radius:8px; border:1px solid #ddd; cursor:pointer;">📘 Tutorial</button>'
        '</a>',
        unsafe_allow_html=True,
    )

st.sidebar.markdown("""---""")

st.sidebar.markdown("Please Cite:")

st.sidebar.markdown("""Borges, R. M., Ferreira, G. de A., Campos, M. M., Teixeira, A. M., Costa, F. das N., & Chagas, F. O. (2024). Data Base similarity (DBsimilarity) of natural products to aid compound identification on MS and NMR pipelines, similarity networking, and more. Phytochemical Analysis, 35(1), 93–101. https://doi.org/10.1002/pca.3277

""")

st.sidebar.markdown("""---""")


# ---------- Sidebar: uploads & options ----------
st.sidebar.header("1) Upload your data")
up_targets = st.sidebar.file_uploader("Target structures (CSV with SMILES or SDF) — multiple allowed", type=["csv","txt","tsv","sdf"], accept_multiple_files=True)
up_annot   = st.sidebar.file_uploader("Annotated table(s) (CSV/TXT) — multiple allowed", type=["csv","txt","tsv"], accept_multiple_files=True)

st.sidebar.header("2) Column mapping")
col_smiles = st.sidebar.text_input("SMILES column name (targets/annotated)", value="SMILES")
col_key    = st.sidebar.text_input("Key column to merge (e.g., file_path / name)", value="file_path")

st.sidebar.header("3) Feature toggles")
run_massql   = st.sidebar.checkbox("Build MassQL (MS1/MS2)", value=True)
run_merge    = st.sidebar.checkbox("Merge targets × annotated (co-occurrence)", value=True)
run_similarity = st.sidebar.checkbox("Compute similarity + network", value=True)
run_descriptors = st.sidebar.checkbox("Mordred descriptors + PCA/Dendro/t-SNE", value=False)

st.sidebar.divider()
st.sidebar.header("Performance / Limits")
max_rows = st.sidebar.number_input("Max rows to process (0 = no cap)", min_value=0, value=0, step=100)
sim_metric = st.sidebar.selectbox("Similarity metric", ["dice","tanimoto"], index=0)


st.sidebar.markdown("---")

# MassQL params (sidebar)
if run_massql:
    st.sidebar.subheader("MassQL params")
    adducts = st.sidebar.multiselect(
        "Adducts",
        list(ADDUCT_SPECS.keys()),
        default=["[M+H]+","[M+Na]+","[M+K]+","[M+NH4]+","[2M+H]+","[M+2H]+2"]
    )
    ppm_ms1 = st.sidebar.number_input("MS1 ppm", 1, 100, 10)
    inten_ms1 = st.sidebar.number_input("MS1 intensity %", 0, 100, 1)
    separate_adducts = st.sidebar.checkbox("Separate query per adduct (MS1)", value=False)

    name_col = st.sidebar.text_input("Name column (in your table)", value="Compound name")
    mass_col = st.sidebar.text_input("Mass column (neutral)", value="MolWeight")
    frag_col = st.sidebar.text_input("Fragments column (m/z list)", value="Fragments")

    ppm_prec  = st.sidebar.number_input("MS2 ppm (precursor)", 1, 100, 10)
    ppm_prod  = st.sidebar.number_input("MS2 ppm (product)",   1, 100, 10)
    inten_ms2 = st.sidebar.number_input("MS2 intensity %",     0, 100, 5)
    card_min  = st.sidebar.number_input("MS2 cardinality min", 1, 10, 1)
    card_max  = st.sidebar.number_input("MS2 cardinality max", 1, 20, 5)
    clamp_card = st.sidebar.checkbox("Clamp cardinality to #frags", value=True)
    # in the sidebar (anywhere after the feature toggles)
    use_merged_for_massql = st.sidebar.checkbox(
        "Use merged table for MassQL input (if available)", value=True
    )




# Descriptor params
if run_descriptors:
    st.sidebar.markdown("---")
    st.sidebar.subheader("Descriptors")
    use_3d = st.sidebar.checkbox("Use 3D descriptors (requires 3D conf gen)", value=False)
    tsne_perp = st.sidebar.slider("**t-SNE perplexity**", 5, 50, 20, 1)
    st.sidebar.markdown("""
    ### 
    **Perplexity** is a smooth measure of the *effective number of neighbors* that t-SNE uses to define the local density around each point.
    - **Low perplexity (e.g. 5–20)** → emphasizes **local structure**  
    - **High perplexity (e.g. 30–100)** → emphasizes **global structure**  
    """)
    tsne_seed = st.sidebar.number_input("t-SNE random_state", 0, 9999, 0, 1)
    st.sidebar.markdown("""
    **`random_state`** controls the *random initialization* of points in the t-SNE algorithm.
    - Each run of t-SNE starts from a random layout.  
      If the `random_state` is **not fixed**, your clusters might appear in different positions or even shapes each time.
    - Setting a fixed `random_state` (e.g. `random_state=42`) ensures **reproducible results** —  
      meaning the same input data always produces the same embedding.
    """)
    
    dendro_thresh = st.sidebar.number_input("**Dendrogram color threshold**", 1, 1_000_000, 250_000, 1_000)

# ---------- Load & assemble targets / annotated ----------
def concat_frames(frames: List[pd.DataFrame]) -> pd.DataFrame:
    frames = [f for f in frames if isinstance(f, pd.DataFrame) and len(f)]
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)

frames_target: List[pd.DataFrame] = []
frames_annot:  List[pd.DataFrame] = []

# Read all target files
if up_targets:
    for f in up_targets:
        if f.type == "chemical/x-mdl-sdfile" or f.name.lower().endswith(".sdf"):
            frames_target.append(read_sdf_to_df(f.read(), f.name))
        else:
            frames_target.append(smart_read_table(f.read(), f.name))

# Read all annotated files
if up_annot:
    for f in up_annot:
        frames_annot.append(smart_read_table(f.read(), f.name))

# Clean + enrich EACH frame (per-file) to get InchiKey early
def _clean_enrich_frames(frames: List[pd.DataFrame], smiles_col: str) -> List[pd.DataFrame]:
    out = []
    for df in frames:
        if not isinstance(df, pd.DataFrame) or df.empty:
            continue
        d = df.copy()
        if smiles_col in d.columns:
            d[smiles_col] = d[smiles_col].replace('', np.nan)
            d = d.dropna(subset=[smiles_col]).drop_duplicates(subset=[smiles_col])
            d = enrich_chem(d, smiles_col)  # <- adds InchiKey
        out.append(d)
    return out

frames_target = _clean_enrich_frames(frames_target, col_smiles)
frames_annot  = _clean_enrich_frames(frames_annot,  col_smiles)

# Concatenated previews (nice for UI)
df_target = concat_frames(frames_target)
df_annot  = concat_frames(frames_annot)

# Apply caps for previews only
if max_rows and len(df_target) > max_rows:
    df_target = df_target.head(max_rows)
if max_rows and len(df_annot) > max_rows:
    df_annot = df_annot.head(max_rows)

st.markdown("### Uploaded tables")
col1, col2 = st.columns(2)
with col1:
    show_df(df_target.head(1), caption="Targets (preview)")
with col2:
    show_df(df_annot.head(1), caption="Annotated (preview)")

# ---------- Auto-merge ASAP using InchiKey when ≥2 files are present ----------
def _frames_with_key(frames: List[pd.DataFrame], key: str) -> List[pd.DataFrame]:
    ok = []
    for d in frames:
        if isinstance(d, pd.DataFrame) and not d.empty and (key in d.columns):
            ok.append(d[[key] + [c for c in d.columns if c != key]])  # keep key + rest
    return ok

# Collect all enriched frames that have InchiKey
key_col_for_merge = "InchiKey"  # <- force InchiKey as the merging key
frames_all = _frames_with_key(frames_target, key_col_for_merge) + _frames_with_key(frames_annot, key_col_for_merge)

if len(frames_all) >= 2:
    merged_auto = merge_by_key(frames_all, key_col_for_merge)
    _set_state_df_if_changed("merged_df", merged_auto)
else:
    if isinstance(st.session_state.get("merged_df"), pd.DataFrame) and not st.session_state["merged_df"].empty:
        _set_state_df_if_changed("merged_df", pd.DataFrame())


# Apply caps
if max_rows and len(df_target) > max_rows:
    df_target = df_target.head(max_rows)
if max_rows and len(df_annot) > max_rows:
    df_annot = df_annot.head(max_rows)

# ---------- Cleaning: drop empty/duplicates; enrich chem ----------
if not df_target.empty and col_smiles in df_target.columns:
    df0 = df_target.copy()
    df0[col_smiles] = df0[col_smiles].replace('', np.nan)
    df0 = df0.dropna(subset=[col_smiles])
    df0 = df0.drop_duplicates(subset=[col_smiles])
    df_target = enrich_chem(df0, col_smiles)
    st.success(f"Targets cleaned: {len(df0)} rows → {len(df_target)} unique SMILES.")

# ---------- Build merged_df early (even if only one table has the key) ----------
dfs_to_merge_early = []
for _df in [df_target, df_annot]:
    if isinstance(_df, pd.DataFrame) and not _df.empty and (col_key in _df.columns):
        dfs_to_merge_early.append(_df)

if dfs_to_merge_early:
    merged_early = merge_by_key(dfs_to_merge_early, col_key)
    _set_state_df_if_changed("merged_df", merged_early)



st.markdown("---")

# ---------- Merge targets × annotated (co-occurrence) ----------

if run_merge and not df_target.empty:
    st.subheader("Merge & co-occurrence")
    dfs_to_merge = [x for x in [df_target, df_annot] if not x.empty and col_key in x.columns]
    if len(dfs_to_merge) >= 1:
        merged = merge_by_key(dfs_to_merge, col_key)
        _set_state_df_if_changed("merged_df", merged)


        # Show columns clearly (no pyarrow)
        st.markdown("**Columns in merged table:**")
        st.code("\n".join(map(str, merged.columns)))

        # Preview table (dropdown)
        with st.expander("Preview merged table (choose how many rows)"):
            nrows = st.selectbox("Rows to show", [5, 10, 20, 50, 100], index=1, key="merged_preview_rows")
            show_df(merged.head(nrows), caption=f"Merged preview (first {nrows} rows)")
            st.caption(f"{merged.shape[0]} rows × {merged.shape[1]} columns")

        # Mark intersection (if two frames given)
        if len(dfs_to_merge) >= 2:
            inter = merged.loc[(merged.get("in_df1",0)==1) & (merged.get("in_df2",0)==1)]
            st.info(f"Co-occurring rows: {len(inter)}")

            # Exports with semicolon separator
            exports["Co-occurring_compounds.csv"] = inter.to_csv(sep=";", index=False).encode()
            exports["Merged_all.csv"] = merged.to_csv(sep=";", index=False).encode()


# ---------- Global merged preview (if available) ----------
merged_df = st.session_state.get("merged_df", pd.DataFrame())
if not merged_df.empty:
    st.subheader("Merged preview (InchiKey-based)")
    with st.expander("Preview merged table (choose how many rows)"):
        nrows = st.selectbox(
            "Rows to show", [5, 10, 20, 50, 100],
            index=1,
            key="merged_preview_rows_global"
        )
        # --- Inline downloads for the merged (InchiKey-based) table ---
        csv_blob = merged_df.to_csv(sep=';', index=False).encode('utf-8')

        # (optional) Excel download
        from io import BytesIO
        _buf = BytesIO()
        try:
            with pd.ExcelWriter(_buf, engine="xlsxwriter") as _xw:
                merged_df.to_excel(_xw, index=False, sheet_name="merged")
            xlsx_data = _buf.getvalue()
            st.download_button(
                "Download merged.xlsx",
                data=xlsx_data,
                file_name="Merged_all.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                key="dl_merged_xlsx_global",
            )
        except Exception as _e:
            st.caption(f"Excel export unavailable ({_e}). Install 'xlsxwriter' to enable.")

        # (optional) also expose in the bottom 'Downloads' section
        exports["Merged_all.csv"] = csv_blob
        if 'xlsx_data' in locals():
            exports["Merged_all.xlsx"] = xlsx_data
        
        
        show_df(merged_df.head(nrows), caption=f"Merged preview (first {nrows} rows)")
        st.caption(f"{merged_df.shape[0]} rows × {merged_df.shape[1]} columns")
else:
    st.info("Merged table not available yet — import at least two files that contain SMILES so we can compute InchiKey.")

st.markdown(f"Combined compounds sum up to: {merged_df.shape[0]}")

# === Co-Occurrence preview (rows present in ALL input tables) ===
import re

co_df = pd.DataFrame()

st.subheader("Co-Occurrence preview")

if merged_df.empty:
    st.info("No merged data yet.")
else:
    in_cols = [c for c in merged_df.columns if re.match(r"^in_df\d+$", c)]
    if len(in_cols) < 2:
        st.caption("Co-occurrence needs ≥2 input tables (found: " + str(len(in_cols)) + ").")
    else:
        co_mask = (merged_df[in_cols] == 1).all(axis=1)
        co_df = merged_df.loc[co_mask].copy()

        with st.expander("Preview co-occurring rows (present in ALL input tables)"):
            nrows = st.selectbox(
                "Rows to show",
                [0, 1, 10, 20, 50, 100],
                index=1,
                key="cooc_preview_rows"  # keep key unique
            )
            show_df(co_df.head(nrows), caption=f"Co-occurring (first {nrows} rows)")
            st.caption(f"{len(co_df)} co-occurring rows out of {len(merged_df)} total. "
                       f"Required columns: {', '.join(in_cols)}")

            # Downloads (CSV ; and optional XLSX)
            co_csv = co_df.to_csv(sep=';', index=False).encode('utf-8')
            st.download_button(
                "Download co-occurring (CSV ;)",
                data=co_csv,
                file_name="Co_occurring_compounds.csv",
                mime="text/csv",
                key="dl_cooc_csv"
            )

            try:
                from io import BytesIO
                _co_buf = BytesIO()
                with pd.ExcelWriter(_co_buf, engine="xlsxwriter") as _xw:
                    co_df.to_excel(_xw, index=False, sheet_name="co_occurrence")
                co_xlsx = _co_buf.getvalue()
                st.download_button(
                    "Download co-occurring.xlsx",
                    data=co_xlsx,
                    file_name="Co_occurring_compounds.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    key="dl_cooc_xlsx"
                )
            except Exception as _e:
                st.caption(f"Excel export unavailable ({_e}). Install 'xlsxwriter' to enable.")

            # Also expose in bottom 'Downloads' section
            exports["Co_occurring_compounds.csv"] = co_csv
            if 'co_xlsx' in locals():
                exports["Co_occurring_compounds.xlsx"] = co_xlsx

st.markdown(f"Co-occurrent compounds sum up to: {co_df.shape[0]}")

if st.sidebar.button("🗑️ Clear all downloads"):
    st.session_state["exports"].clear()
    st.sidebar.success("Exports cleared.")

# ---------- Check Information from the Co-occurrent Compounds ----------
st.markdown("---")
st.subheader("Check Information from the Co-occurrent Compounds")
st.markdown("""
If the added (or any added) csv file has any caracteristics that is worth mentions, we should select the right Column header to create a wordcloud"
""")

# ---------- Export format to NMRfilter ----------
st.markdown("---")
st.subheader("Export format to NMRfilter for NMR-based dereplication")
st.markdown("""
NMRfilter - check out: 
""")

# ---------- Export format to custom database using MZMine ----------
st.markdown("---")
st.subheader("Export format to custom database using MZMine (MS1 comparison only)")
st.markdown("""
MZMine - check out: 
""")

# ---------- MassQL builders ----------
st.markdown("---")

def _pick_mass_col(df: pd.DataFrame, preferred: str) -> tuple[str, str|None]:
    candidates = [preferred, "MolWeight", "ExactMass", "MonoisotopicMass", "MonoIsotopicMass", "Mass"]
    for c in candidates:
        if c in df.columns:
            return c, (None if c == preferred else f"Using '{c}' for mass (fallback).")
    return "", "No mass-like column found. Add 'MolWeight' or map your mass column."

def _ensure_name_col(df: pd.DataFrame, wanted: str, smiles_col: str) -> tuple[str, str|None]:
    if wanted in df.columns:
        return wanted, None
    for c in ["Name","Compound","Compound_Name","CompoundName","InchiKey",smiles_col]:
        if c in df.columns:
            auto = "__AutoName"
            s = df[c].astype(str).str.strip()
            s = s.where(s != "", None)
            df[auto] = s.fillna(pd.Series([f"entry_{i+1}" for i in range(len(df))], index=df.index))
            return auto, f"No '{wanted}' found; using '{c}' (or generated 'entry_i')."
    auto = "__AutoName"
    df[auto] = [f"entry_{i+1}" for i in range(len(df))]
    return auto, f"No '{wanted}' found; generated sequential names."

if run_massql:
    st.subheader("MassQL builders")

    # choose the source: merged -> target -> annotated
    source_df = pd.DataFrame()
    merged_for_massql = st.session_state.get("merged_df", pd.DataFrame())
    if use_merged_for_massql and not merged_for_massql.empty:
        source_df = merged_for_massql.copy()
        st.caption("Using **merged** table (InchiKey merge) for MassQL input.")
    elif not df_target.empty:
        source_df = df_target.copy()
        st.caption("Using **targets** table for MassQL input.")
    elif not df_annot.empty:
        source_df = df_annot.copy()
        st.caption("Using **annotated** table for MassQL input.")



    if source_df.empty:
        st.warning("No input table available yet. Upload a CSV/SDF first.")
    else:
        # st.code("\n".join(map(str, source_df.columns)))  # uncomment if you also want to echo columns here

        # Pick/ensure columns
        use_mass_col, mass_note = _pick_mass_col(source_df, mass_col)
        use_name_col, name_note = _ensure_name_col(source_df, name_col, col_smiles)
        if mass_note: st.info(mass_note)
        if name_note: st.info(name_note)

        if not use_mass_col:
            st.error("Cannot build MS1: no mass column available.")
        else:
            if st.button("Build MS1 & (optional) MS2 queries"):
                with st.spinner("Building MS1 queries..."):
                    q_ms1 = build_ms1_queries(
                        source_df, use_name_col, use_mass_col, adducts,
                        ppm_ms1, inten_ms1, 5, separate_adducts
                    )

                q_ms2 = {}
                has_ms2 = (frag_col in source_df.columns) and source_df[frag_col].notna().any()
                if has_ms2:
                    with st.spinner("Building MS2 queries..."):
                        q_ms2 = build_ms2_queries(
                            source_df, use_name_col, use_mass_col, frag_col,
                            adducts, ppm_prec, ppm_prod, inten_ms2,
                            5, card_min, card_max, clamp_card, True
                        )
                else:
                    st.info("No 'Fragments' column found or it is empty — skipping MS2 (MS1 only).")

                ms1_buf = io.StringIO(); i = 1
                for k, v in q_ms1.items():
                    if isinstance(v, dict):
                        for ad, q in v.items():
                            ms1_buf.write(f"### {i}. {k} {ad} ###\n{q}\n\n"); i += 1
                    else:
                        ms1_buf.write(f"### {i}. {k} ###\n{v}\n\n"); i += 1
                exports["MS1_queries_by_compound.txt"] = ms1_buf.getvalue().encode()

                if q_ms2:
                    ms2_buf = io.StringIO(); i = 1
                    for k, q in q_ms2.items():
                        ms2_buf.write(f"### {i}. {k} ###\n{q}\n\n"); i += 1
                    exports["MS2_queries_by_compound.txt"] = ms2_buf.getvalue().encode()

                st.success("MassQL file(s) ready. See Downloads below.")


# ---------- Similarity + network ----------
st.markdown("---")
if run_similarity and not df_target.empty:
    st.subheader("Similarity network (Morgan FP)")

    # keep controls visible (outside the dropdown)
    sim_thresh = st.slider("Edge threshold (similarity)", 0.0, 1.0, 0.85, 0.01, key="sim_thresh_main")

    # compute similarities and graph
    smiles_list = df_target[col_smiles].astype(str).tolist()
    inchi_list  = df_target.get("InchiKey", pd.Series([f"mol_{i+1}" for i in range(len(smiles_list))])).astype(str).tolist()
    fps = fingerprints(smiles_list, radius=2, nBits=2048)
    sim_df = pairwise_similarity(inchi_list, fps, metric=sim_metric)

    # edges and graph
    edges = sim_df[sim_df["CORRELATION"] >= float(sim_thresh)].copy()
    G = nx.from_pandas_edgelist(edges, "SOURCE","TARGET")
    nodes_all = set(inchi_list)
    G.add_nodes_from(nodes_all - set(G.nodes()))

    # layout + traces
    pos = nx.spring_layout(G, seed=42, k=0.8)
    x, y, text, deg = [], [], [], []
    for n in G.nodes():
        xx, yy = pos[n]
        x.append(xx); y.append(yy)
        text.append(n)
        deg.append(G.degree[n])
    node_trace = go.Scatter(x=x, y=y, text=text, mode='markers', marker=dict(size=[6+2*d for d in deg]))
    xe, ye = [], []
    for u, v in G.edges():
        xe += [pos[u][0], pos[v][0], None]
        ye += [pos[u][1], pos[v][1], None]
    edge_trace = go.Scatter(x=xe, y=ye, mode='lines', hoverinfo='none')
    fig_net = go.Figure(data=[edge_trace, node_trace])
    fig_net.update_layout(title="Similarity network", showlegend=False, margin=dict(l=0, r=0, t=40, b=0))

    # data for downloads
    edges_csv_blob = edges.to_csv(sep=";", index=False).encode("utf-8")
    net_html_blob = pio.to_html(fig_net, include_plotlyjs='cdn', full_html=True).encode("utf-8")

    G2 = nx.from_pandas_edgelist(edges, "SOURCE", "TARGET")
    iso = list(nodes_all - set(G2.nodes()))
    df_iso = df_target[df_target["InchiKey"].isin(iso)] if "InchiKey" in df_target.columns else pd.DataFrame({"InchiKey": iso})
    iso_csv_blob = df_iso.to_csv(index=False).encode("utf-8")

    # put the heavy preview inside a dropdown
    with st.expander("Show network and sample rows", expanded=False):
        show_df(sim_df.head(10), caption="Similarity rows (sample)")
        st.caption(f"Edges ≥ {sim_thresh:.2f}: {len(edges)}  |  Isolated nodes: {len(iso)}")
        st.plotly_chart(fig_net, use_container_width=True)

    # keep download buttons outside the dropdown (always visible)
    c1, c2, c3 = st.columns(3)
    with c1:
        st.download_button(
            f"Download edges CSV (≥ {sim_thresh:.2f})",
            data=edges_csv_blob,
            file_name=f"DB_compounds_Similarity_{sim_thresh:.2f}.csv",
            key="dl_edges_csv_main"
        )
    with c2:
        st.download_button(
            "Download network HTML",
            data=net_html_blob,
            file_name="similarity_network.html",
            key="dl_network_html_main"
        )
    with c3:
        st.download_button(
            "Download isolated_nodes.csv",
            data=iso_csv_blob,
            file_name="isolated_nodes.csv",
            key="dl_iso_csv_main"
        )

    # also expose in the bottom 'Downloads' section
    exports[f"DB_compounds_Similarity_{sim_thresh:.2f}.csv"] = edges_csv_blob
    exports["similarity_network.html"] = net_html_blob
    exports["isolated_nodes.csv"] = iso_csv_blob

# ---------- Check Information from Similar Compounds ----------
st.markdown("---")
st.subheader("Check Information from Similar Compounds")
st.markdown("""
If the added (or any added) csv file has any caracteristics that is worth mentions, we should select the right Column header to create a wordcloud"
""")


# ---------- Descriptors + PCA / Dendrogram / t-SNE ----------
st.markdown("---")

def sanitize_matrix_for_modeling(df_vals: pd.DataFrame) -> pd.DataFrame:
    """Keep numeric, finite; drop constant-variance columns."""
    X = df_vals.select_dtypes(include=[np.number]).copy()
    if X.empty:
        return X
    # keep only finite values; impute the rest to 0
    X = X.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    # drop constant columns
    variances = X.var(axis=0)
    X = X.loc[:, variances > 0.0]
    return X

if run_descriptors and not df_target.empty and col_smiles in df_target.columns:
    st.subheader("Descriptors → PCA / Dendrogram / t-SNE")

    smiles = df_target[col_smiles].astype(str).tolist()
    mols = [mol_from_smiles(s) for s in smiles]

    with st.spinner("Calculating Mordred descriptors (this may take a while)…"):
        try:
            calc = Calculator(descriptors, ignore_3D=not use_3d)
            df_desc = calc.pandas(mols)
        except Exception as e:
            st.error(f"Mordred failed: {e}")
            df_desc = pd.DataFrame()

    if df_desc.empty:
        st.info("No descriptors available (empty table). Skipping PCA / Dendrogram / t-SNE.")
    else:
        # 1) Sanitize → numeric, finite, non-constant
        df_vals_raw = df_desc.select_dtypes(include=[np.number])
        df_vals = sanitize_matrix_for_modeling(df_vals_raw)

        # Export cleaned descriptors (useful for debugging)
        try:
            exports["descriptors_clean.csv"] = df_vals.to_csv(index=False).encode("utf-8")
        except Exception:
            pass

        # 2) (Optional) cap features to top-k by variance BEFORE scaling
        k = min(1000, df_vals.shape[1])  # expose in sidebar later if you want
        if k < df_vals.shape[1]:
            topk = df_vals.var().sort_values(ascending=False).index[:k]
            df_vals = df_vals[topk]

        # 3) Re-evaluate shapes AFTER capping
        n_samples, n_features = df_vals.shape
        st.caption(f"Descriptors after cleaning: {n_samples} samples × {n_features} features")

        if n_samples < 2 or n_features < 1:
            st.info("Not enough data for PCA/Dendrogram/t-SNE (need ≥2 samples and ≥1 feature).")
        else:
            # 4) Standardize ON the final df_vals (use everywhere below)
            try:
                X = StandardScaler(with_mean=True, with_std=True).fit_transform(df_vals.values)
            except Exception as e:
                st.warning(f"Standardization failed: {e}")
                X = df_vals.values  # fallback

            # Optional IDs for hover labels (aligns with df_target order)
            ids = df_target.get("InchiKey", pd.Series([f"mol_{i+1}" for i in range(n_samples)])).astype(str).tolist()

            # ----------------- PCA (2D) -----------------
            try:
                ncomp = int(min(2, n_features, n_samples))
                if ncomp >= 1:
                    # Fit PCA
                    pca = PCA(n_components=ncomp, svd_solver="auto")
                    scores = pca.fit_transform(X)

                    # ----- Optional color by in_df* flags from merged_df -----
                    pca_color = None
                    merged_for_colors = st.session_state.get("merged_df", pd.DataFrame())
                    if not merged_for_colors.empty:
                        import re as _re
                        in_cols_all = [c for c in merged_for_colors.columns if _re.match(r"^in_df\d+$", c)]

                        # Let user include constant columns so they can pick in_df1 if it's all 1s
                        show_constant = st.checkbox(
                            "Show all in_df… columns (even if constant on PCA set)",
                            value=False,
                            key="pca_show_constant_cols",
                            help="If off, only columns that vary (both 0 and 1) are listed."
                        )

                        # Join key: prefer InchiKey, else col_key if present in both
                        join_key = "InchiKey" if "InchiKey" in merged_for_colors.columns else (
                            col_key if (col_key in merged_for_colors.columns and col_key in df_target.columns) else None
                        )

                        if join_key is not None and in_cols_all:
                            # Align flags to PCA order
                            left  = pd.DataFrame({join_key: df_target.get(join_key, pd.Series(ids)).values})
                            right = merged_for_colors[[join_key] + in_cols_all].drop_duplicates(subset=[join_key])
                            merged_colors = left.merge(right, on=join_key, how="left")

                            # Which columns vary on this subset?
                            varying_cols = []
                            for c in in_cols_all:
                                vals = merged_colors[c].fillna(0).astype(int)
                                if vals.nunique() > 1:
                                    varying_cols.append(c)

                            selectable = in_cols_all if show_constant else varying_cols
                            if selectable:
                                default_idx = selectable.index("in_df2") if "in_df2" in selectable else 0
                                color_col = st.selectbox("Color PCA by", selectable, index=default_idx, key="pca_color_in_df")
                                pca_color = merged_colors[color_col].fillna(0).astype(int).astype(str)  # "0"/"1"
                                if color_col == "in_df1" and not show_constant and "in_df1" not in varying_cols:
                                    st.caption("Note: `in_df1` is constant (=1) on the PCA set; enable the checkbox above to select it.")
                            else:
                                st.info(
                                    "No `in_df…` column varies across the PCA points. "
                                    "This usually happens because PCA uses only `df_target` rows and `in_df1` flags membership in `df_target` (all 1). "
                                    "Turn on the checkbox to select constant columns, or pick another column like `in_df2`."
                                )
                        else:
                            st.caption("Could not align colors for PCA (no common key or no `in_df…` columns).")

                    # ----- Create PCA scatter (with/without color) -----
                    if pca_color is not None:
                        fig_pca = px.scatter(
                            x=scores[:, 0],
                            y=(scores[:, 1] if ncomp > 1 else np.zeros_like(scores[:, 0])),
                            title=f"PCA ({ncomp} component{'s' if ncomp>1 else ''})",
                            labels={"x": "PC1", "y": ("PC2" if ncomp > 1 else "PC2 (zero)")},
                            hover_name=ids,
                            color=pca_color
                        )
                    else:
                        fig_pca = px.scatter(
                            x=scores[:, 0],
                            y=(scores[:, 1] if ncomp > 1 else np.zeros_like(scores[:, 0])),
                            title=f"PCA ({ncomp} component{'s' if ncomp>1 else ''})",
                            labels={"x": "PC1", "y": ("PC2" if ncomp > 1 else "PC2 (zero)")},
                            hover_name=ids
                        )

                    st.plotly_chart(fig_pca, use_container_width=True, #key="pca_chart"
                    )

                    # Exports
                    exports["pca_scores.csv"] = (
                        pd.DataFrame(scores[:, :max(1, ncomp)], columns=["PC1","PC2"][:ncomp]).to_csv(index=False).encode()
                    )
                    exports["pca.html"] = pio.to_html(fig_pca, include_plotlyjs='cdn', full_html=True).encode()
                else:
                    st.info("PCA skipped: fewer than 1 usable component.")
            except Exception as e:
                st.warning(f"PCA failed: {e}")


            # ----------------- Dendrogram: view or save (HTML/PNG) -----------------
            # === Dendrogram controls (compute and store) ===
            with st.expander("Dendrogram controls (compute clusters)"):
                linkage_method = st.selectbox(
                    "Linkage", ["ward","average","complete","single","weighted","centroid","median"], index=0
                )
                dist_metric = st.selectbox(
                    "Distance metric", ["euclidean","cityblock","cosine","correlation","chebyshev","canberra","braycurtis","minkowski"], index=0
                )
                # Ward requires euclidean
                if linkage_method == "ward" and dist_metric != "euclidean":
                    st.info("Linkage 'ward' requires Euclidean; overriding metric to 'euclidean'.")
                    dist_metric = "euclidean"

                # compute condensed distance + linkage
                from scipy.spatial.distance import pdist
                from scipy.cluster.hierarchy import linkage as _linkage, fcluster as _fcluster

                D = pdist(X, metric=dist_metric)              # X is your standardized matrix
                Z = _linkage(D, method=linkage_method)

                # slider to define color threshold
                thr = st.slider("Cut (color threshold)", float(D.min()) if len(D) else 0.0,
                                float(D.max()) if len(D) else 1.0,
                                value=float(np.percentile(D, 75)) if len(D) else 0.5, step=0.001)

                st.session_state["_dendro_Z"] = Z
                st.session_state["_dendro_thr"] = thr

                # flat clusters consistent with the dendrogram cut
                labels_ = _fcluster(Z, t=float(thr), criterion="distance")
                st.session_state["_cluster_labels"] = labels_
                st.caption(f"Formed {len(np.unique(labels_))} clusters at cut = {thr:.3g}.")

            
            with st.expander("Dendrogram (view / save)"):
                Z_use = st.session_state.get("_dendro_Z", None)
                thr_use = st.session_state.get("_dendro_thr", None)

                if Z_use is None or thr_use is None:
                    st.info("Adjust the dendrogram slider once to define the cut (no saved dendrogram state yet).")
                else:
                    # Labels (fallback if InchiKey missing)
                    ids_re = df_target.get(
                        "InchiKey",
                        pd.Series([f"mol_{i+1}" for i in range(n_samples)])
                    ).astype(str).tolist()

                    # Controls
                    c1, c2, c3 = st.columns([1,1,1])
                    with c1:
                        fig_w = st.number_input("Width (inches)", 6.0, 20.0, 10.0, 0.5)
                    with c2:
                        fig_h = st.number_input("Height (inches)", 3.0, 20.0, 6.0, 0.5)
                    with c3:
                        dpi = st.number_input("DPI (PNG only)", 100, 600, 200, 25)

                    # Build Matplotlib dendrogram (shown in app)
                    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=100)
                    dn = dendrogram(
                        Z_use,
                        color_threshold=float(thr_use),
                        labels=ids_re,
                        leaf_rotation=90,
                        leaf_font_size=8,
                        above_threshold_color="grey"
                    )
                    ax.set_title(f"Dendrogram (cut @ {float(thr_use):.3g})")
                    ax.set_ylabel("Distance")
                    ax.set_xlabel("Samples")
                    st.pyplot(fig, clear_figure=False)

                    # --- Downloads ---
                    # 1) Interactive HTML (mpld3)
                    # inside the dendrogram "Downloads" block, after you draw the figure:
                    if HAS_MPLD3:
                        html_str = mpld3.fig_to_html(fig)
                        st.download_button("⬇️ Download dendrogram (HTML, interactive)",
                                           data=html_str, file_name="dendrogram.html", mime="text/html")
                        exports["dendrogram.html"] = html_str.encode("utf-8")
                    else:
                        st.caption("Install `mpld3` to enable interactive HTML export: `pip install mpld3`.")
                    
                    # PNG (static)
                    png_buf = io.BytesIO()
                    fig.savefig(png_buf, format="png", dpi=int(dpi), bbox_inches="tight")
                    png_buf.seek(0)
                    st.download_button("⬇️ Download dendrogram (PNG)",
                                       data=png_buf, file_name="dendrogram.png", mime="image/png")
                    exports["dendrogram.png"] = png_buf.getvalue()


                    # Optional: keep latest export synchronized
                    exports["dendrogram.html"] = html_str.encode("utf-8")
                    #exports["dendrogram.png"] = png_buf.getvalue()


            # ----------------- t-SNE -----------------
            try:
                if n_samples < 3:
                    st.info("t-SNE skipped: need ≥3 samples.")
                else:
                    # keep perplexity valid and reasonable
                    # guideline: < (n_samples - 1); often ~ (n_samples-1)/3
                    max_perp = max(5, min(int(tsne_perp), n_samples - 1, max(5, (n_samples - 1)//3)))
                    tsne = TSNE(
                        n_components=2,
                        learning_rate='auto',
                        random_state=int(tsne_seed),
                        init='random',
                        perplexity=max_perp
                    )
                    ts = tsne.fit_transform(X)  # SAME standardized X
                    
                    # If we have clusters from the dendrogram, use them as categorical colors
                    tsne_color = st.session_state.get("_cluster_labels", None)
                    if tsne_color is not None and len(tsne_color) == n_samples:
                        fig_ts = px.scatter(
                            x=ts[:, 0], y=ts[:, 1],
                            title=f"t-SNE (perplexity={max_perp})",
                            labels={'x': 'tSNE-1', 'y': 'tSNE-2'},
                            hover_name=ids,
                            color=pd.Series(tsne_color, dtype=str)
                        )
                    else:
                        fig_ts = px.scatter(
                            x=ts[:, 0], y=ts[:, 1],
                            title=f"t-SNE (perplexity={max_perp})",
                            labels={'x': 'tSNE-1', 'y': 'tSNE-2'},
                            hover_name=ids
                        )

                    st.plotly_chart(fig_ts, use_container_width=True, #key="tsne_chart"
                    )

                    
                    exports["tsne_scores.csv"] = (
                        pd.DataFrame(ts, columns=["tSNE1","tSNE2"]).to_csv(index=False).encode()
                    )
                    exports["tsne.html"] = pio.to_html(fig_ts, include_plotlyjs='cdn', full_html=True).encode()
            except Exception as e:
                st.warning(f"t-SNE failed: {e}")



# ---------- Exports section ----------
if exports:
    st.subheader("Downloads")
    cols = st.columns(3)
    i = 0
    for fname, blob in exports.items():
        with cols[i%3]:
            st.download_button(f"Download {fname}", data=blob, file_name=fname)
        i += 1

st.markdown(
    """
    ---
    **Notes**
    - CSV reader accepts `;`, `TAB`, or `,` (auto). SDF is supported (converted to a DataFrame with SMILES).
    - Map your SMILES and key columns in the sidebar.
    - MassQL files include MS1/ MS2 builders with adducts and ppm/intensity/cardinality controls.
    - Similarity uses Morgan FP (radius=2, 2048 bits) with Dice or Tanimoto; network HTML export included.
    - Descriptors can be heavy; enable only if needed. All plots are downloadable (HTML/CSV where applicable).
    """
)

# --------------------------
# requirements.txt suggestion:
# streamlit
# rdkit-pypi
# pandas
# numpy
# plotly
# networkx
# mordred
# scikit-learn
# scipy
# matplotlib








