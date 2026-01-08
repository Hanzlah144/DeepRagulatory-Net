# DeepRegulatoryNet ✅

**circRNA–miRNA–mRNA analysis pipeline**

A reproducible pipeline to predict circRNA–miRNA binding sites, build regulatory networks, perform gene set enrichment, construct PPI networks, and analyze drug–gene interactions.

---

## 🔧 Supported Python & Environment

- Python 3.8+ (tested on 3.8–3.11)
- Recommended: create and activate a virtual environment (venv or conda)

Example (venv):
```bash
python -m venv .venv
# Windows
.\.venv\Scripts\activate
# macOS / Linux
source .venv/bin/activate
```

---

## 📦 Dependencies

A proposed `requirements.txt` is provided for convenience. These are the packages used by the pipeline:

- pandas >= 1.5
- numpy >= 1.24
- networkx >= 3.0
- scikit-learn >= 1.1
- catboost >= 1.2
- openpyxl >= 3.0
- requests >= 2.28
- urllib3 >= 1.26
- matplotlib >= 3.5
- seaborn >= 0.12
- tqdm >= 4.64

Install the dependencies with:

```bash
pip install -r requirements.txt
```

---

## 🗂️ Repository layout

- `DeepRegulatoryNet.py` - pipeline entrypoint (CLI)
- `src/` - pipeline implementation modules
- `Model_Files/` - required ML artifacts (model, encoder, scaler)
- `Test_Data/` - example input files
- `output/` (created at runtime) - pipeline outputs (CSVs, Excel, plots)
- `temp/` (created at runtime) - temporary files

---

## 📋 Input file formats & validation

- `--circ`: circRNA IDs (one per line); IDs should start with `hsa_circ_`.
- `--mirna`: miRNA IDs (one per line); IDs should start with `hsa-miR`.
- `--deg`: Differentially expressed genes (one gene symbol per line).

The pipeline validates ID prefixes for `--circ` and `--mirna` and will raise an error on invalid lines.

Example test files are included in `Test_Data/`.

---

## 🚀 Quick start

1. Create and activate a virtual environment.
2. Install dependencies:

```bash
pip install -r requirements.txt
```

3. Run the pipeline with example data:

```bash
python DeepRegulatoryNet.py \
  --circ Test_Data/DEcircRNA.txt \
  --mirna Test_Data/DEmiRNA.txt \
  --deg Test_Data/DEG.txt
```

- Add `--debug` to enable verbose logging and stack traces for troubleshooting.
- Run `python DeepRegulatoryNet.py -h` for detailed help and examples.

Windows absolute-path example:

```powershell
python .\DeepRegulatoryNet.py --circ .\Test_Data\DEcircRNA.txt --mirna .\Test_Data\DEmiRNA.txt --deg .\Test_Data\DEG.txt
```

---

## 📁 Expected outputs

- `output/overlapping_genes.csv` — overlapping mRNA results
- `output/enrichment_results/` — enrichment analysis outputs
- `output/hub_genes.csv` — hub genes used for drug analysis
- `pipeline.log` — detailed run log
- `output/*.xlsx` — comprehensive Excel outputs if generated
- Network files (if constructed)

---

## ✅ Notes, troubleshooting & tips

- Missing model files: Ensure `Model_Files/` contains:
  - `calibrated_catboost_site_type_model.pkl`
  - `label_encoder.pkl`
  - `robust_scaler.pkl`

  The pipeline will exit with `FileNotFoundError` if any are missing.

- If you see unexpected behavior, run with `--debug` to see full tracebacks in `pipeline.log`.

- Outputs are written to `output/`. If `temp/` or `output/` exist they will be removed at pipeline start (to ensure reproducible runs).

- To reproduce results or run a single part of the pipeline, consider importing modules from `src/` and calling functions directly in a Python session or a Jupyter notebook.

---

## 📋 Contributing

- Fork the repository, create a feature branch, add tests, and open a pull request.
- Add any new dependencies to `requirements.txt` and update this README.

---

## 📝 License

This project is provided under the MIT License — see `LICENSE`.

---