# DBsimilarity (Streamlit)

**Structure Mining • MassQL Builders • Similarity Networks • Descriptors (PCA/Dendrogram/t-SNE)**

[![Streamlit App](https://img.shields.io/badge/Streamlit-Launch%20App-red?logo=streamlit)](https://share.streamlit.io/your-user/DBsimilarity_st)
[![Python](https://img.shields.io/badge/Python-3.11-blue)](#) [![RDKit](https://img.shields.io/badge/RDKit-2024.03.5-green)](#) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](#license)

DBsimilarity is a no-loss Streamlit port of our chemoinformatics toolkit for natural products. It ingests target structures and annotated tables, cleans/enriches them, merges by keys (InChIKey), builds **MassQL** queries (MS1/MS2) with adduct/ppm/cardinality control, computes **Morgan fingerprint** similarities and **interactive networks**, and optionally runs **Mordred descriptors → PCA, hierarchical clustering (dendrogram), and t-SNE**. All artifacts are downloadable (CSV/HTML/PNG/XLSX).

---

## Features

* **Data intake (multi-upload)**: CSV/TXT/TSV and SDF (auto SMILES extraction when missing).
* **Cleaning & enrichment**: drop empty/duplicate SMILES; add **MolFormula, InChI, InChIKey, ExactMass** and common **adduct m/z**.
* **Merge (co-occurrence)**: tolerant key preprocessing; automatic **InChIKey** merge; `in_df1/in_df2/...` flags; **Co-Occurrence preview**.
* **MassQL builders**:

  * **MS1** with configurable **adducts, ppm, intensity**, single or per-adduct queries.
  * **MS2** with **precursor ppm**, **fragment list parsing**, **product ppm**, and **CARDINALITY** range.
* **Similarity networks**:

  * Morgan fingerprints (r=2, 2048 bits), **Dice/Tanimoto** similarity, thresholded edges, isolated-node report.
  * Interactive Plotly network + **HTML export**.
* **Descriptors & ML** *(optional)*:

  * **Mordred 2D/3D**, PCA (Plotly), **Dendrogram** (Matplotlib + interactive HTML via mpld3), **Agglomerative labels**, **t-SNE** (Plotly).
* **One-click exports**: CSV, XLSX, HTML, and PNG for dendrogram.
* **No pyarrow hard-dep** for tables (HTML rendering fallback).

---

## Try it online

> Replace this link after you deploy on Streamlit Cloud.

**App:** [https://share.streamlit.io/**your-user**/**DBsimilarity_st](https://share.streamlit.io/**your-user**/**DBsimilarity_st)**
**Tutorial/Docs:** (this repo)
**Demo data:** see `sample_data/` (optional)

---

## Installation (local)

```bash
git clone https://github.com/your-user/DBsimilarity_st.git
cd DBsimilarity_st
# (Recommended) create a fresh env with Python 3.11
# conda create -n dbsim python=3.11 -y && conda activate dbsim
pip install -r requirements.txt
streamlit run app.py
```

**Requirements** (pinned for Streamlit Cloud):

```
streamlit>=1.38
pandas>=2.2
numpy>=1.26
rdkit-pypi==2024.03.5
mordred==1.2.0
plotly>=5.24
networkx>=3.2
scikit-learn>=1.5
scipy>=1.11
matplotlib>=3.8
mpld3>=0.5.9
xlsxwriter>=3.2
pillow>=10.4
```

`runtime.txt`:

```
python-3.11
```

Optional Streamlit config: `.streamlit/config.toml`

```toml
[server]
headless = true
fileWatcherType = "auto"
maxUploadSize = 200
[theme]
base = "dark"
```

---

## Deploy to Streamlit Cloud

1. Push this repo to GitHub with the structure:

```
DBsimilarity_st/
├─ app.py
├─ requirements.txt
├─ runtime.txt
├─ .streamlit/
│  └─ config.toml
└─ static/
   ├─ LAABio.png
   └─ DBsimilarity.png
```

2. Go to **share.streamlit.io** → **New app** → select the repo.
3. **Main file path**: `app.py` • **Python**: 3.11 (from `runtime.txt`).
4. Deploy.

---

## Input formats

### CSV/TXT/TSV (semicolon, tab, or comma)

* Required column for structures: **`SMILES`** (case-insensitive resolver is used).
* Recommended columns (auto-detected if present): `MolWeight`, `ExactMass`, `Name`, `Fragments`, `file_path` (custom merge key).
* App computes: `MolFormula`, `Inchi`, `InchiKey`, `ExactMass`, `mz [adduct]`.

**Minimal example (`;` as separator):**

```csv
SMILES;Name;MolWeight;Fragments;file_path
CCO;Ethanol;46.07;31,45,29;lib1.smi
CC(=O)O;Acetic acid;60.05;43,45;lib1.smi
```

### SDF

* SDF properties are embedded; if no SMILES field is present, the app derives SMILES from the molecule block.

---

## How to use (workflow)

1. **Upload** one or more target tables (CSV/SDF) and optionally annotated tables (CSV/TXT).
2. Map **SMILES** and **Key** (e.g., `file_path`) in the sidebar.
3. *(Optional)* Tick **Merge** to generate the combined table and **Co-Occurrence** preview.
4. *(Optional)* Build **MassQL**:

   * Pick **adducts** and **ppm/intensity**; provide `Mass` (falls back to `MolWeight`/`ExactMass`).
   * Add `Fragments` if available to get **MS2** with **CARDINALITY**.
5. *(Optional)* Compute **Similarity network**:

   * Choose **Dice** or **Tanimoto**, set **edge threshold**, download **edges CSV** and **network HTML**.
6. *(Optional)* Enable **Descriptors**:

   * Calculate Mordred, run **PCA**, **Dendrogram** (choose linkage/metric; HTML/PNG export), and **t-SNE** (perplexity auto-kept valid).
7. Use the **Downloads** section to retrieve all artifacts.

---

## Screenshots

> Replace with your PNGs or GIFs.

* Upload & cleaning
  `![upload](./docs/img/upload.png)`

* Merge & co-occurrence
  `![merge](./docs/img/merge.png)`

* Similarity network (interactive)
  `![network](./docs/img/network.png)`

* PCA / Dendrogram / t-SNE
  `![ml](./docs/img/ml.png)`

---

## Notes & Tips

* **InChI/InChIKey** generation is wrapped in safe try/except to avoid rare backend issues.
* **Table rendering** uses HTML to avoid pyarrow hard dependencies.
* **Ward linkage** requires **Euclidean** distance; the app enforces this automatically.
* **t-SNE**: perplexity is constrained to a sensible range (`< n_samples`) for stability.
* **Big datasets**: use the **Max rows** cap and disable Descriptors if only exploring networks/MassQL.

---

## Citation

> Borges, R. M., Ferreira, G. de A., Campos, M. M., Teixeira, A. M., Costa, F. N., & Chagas, F. O. (2024).
> **Data Base similarity (DBsimilarity) of natural products to aid compound identification on MS and NMR pipelines, similarity networking, and more.** *Phytochemical Analysis*, 35(1), 93–101. [https://doi.org/10.1002/pca.3277](https://doi.org/10.1002/pca.3277)

---

## Contributing

Pull requests are welcome!
Please open an issue describing the change (feature/bug) and attach a minimal test file (CSV/SDF) when relevant.

---

## Troubleshooting

* **App won’t start on Cloud**: verify `runtime.txt` is `python-3.11` and dependencies match `requirements.txt`.
* **RDKit import error**: ensure you’re using **rdkit-pypi==2024.03.5** (works on Streamlit Cloud).
* **Blank plots**: check that your **SMILES** column is correctly mapped and non-empty; duplicates are removed.
* **No MS2 built**: ensure a `Fragments` column exists (numbers separated by `,` or `;` or spaces).

---

## Contact & Support

**Author:** Ricardo M. Borges & LAABio-IPPN-UFRJ
**Contact:** [ricardo_mborges@yahoo.com.br](mailto:ricardo_mborges@yahoo.com.br)

If this tool helps your work, consider supporting development:

[![Donate with PayPal](https://www.paypalobjects.com/en_US/i/btn/btn_donate_SM.gif)](https://www.paypal.com/donate/?business=2FYTFNDV4F2D4&no_recurring=0&item_name=Support+with+%245+→+Send+receipt+to+tlc2chrom.app@gmail.com+with+your+login+email+→+Access+within+24h!&currency_code=USD)

---

## License

MIT License — see [LICENSE](./LICENSE).
