# SPEC2SMILES

**Predict high‑resolution EI mass spectra directly from SMILES** – a full‑stack application that couples a trained machine‑learning model (scikit‑learn + RDKit) to a lightweight Flask API and a Svelte single‑page UI.

---

## ✨  Features

| Layer                       | Highlights                                                                                                                                                                                                       |
| --------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Backend**                 | Flask 2.x REST API, Gunicorn production server, on‑the‑fly RDKit featurisation, trained regressor, JSON output with spectrum & peak list.                                                                        |
| **Pre‑processing pipeline** | Parallelised (joblib) mass‑spectrum featuriser, automatic feature schema builder, variance/NaN filtering, scalable StandardScaler wrapper, artefact writer (`feature_preprocessor.pkl`, `feature_mapping.json`). |
| **Model**                   | Any scikit‑learn compatible regressor (default: RandomForest) saved with `joblib`.                                                                                                                               |
| **Frontend**                | Svelte + Plotly.js for interactive spectrum plots, simple form for SMILES input, fully static bundle served by Nginx.                                                                                            |
| **Dev Ops**                 | Multi‑arch Docker build (arm64 & amd64), two‑service `docker‑compose`, automatic proxying `/api/*` → backend.                                                                                                    |

---

## 📁  Repository layout

```
.
├── backend/               # Flask API + ML artefacts
│   ├── app.py             # Entry‑point (creates Flask app)
│   ├── models/            # spectrum_predictor.pkl, feature_preprocessor.pkl, feature_mapping.json
│   ├── templates/         # Jinja templates (fallback HTML test page)
│   └── ...                # src modules (feature_preprocessor.py, model_handler.py, ...)
├── src/                   # Svelte source
│   └── components/        # UI components
├── public/                # Compiled SPA is copied here (public/build)
├── Dockerfile.backend     # Build Python API image
├── Dockerfile.frontend    # Build static SPA + nginx image
├── docker-compose.yml     # One‑shot dev / prod stack
└── README.md
```

---

## ⚙️  Requirements (manual dev install)

| Tool   | Tested version |
| ------ | -------------- |
| Python | 3.10           |
| Node   | 20.x LTS       |
| RDKit  | 2023.09        |
| npm    | 10.x           |

> **Tip:** Use the provided Docker setup if you don’t want to compile RDKit locally.

Install Python deps (backend):

```bash
cd backend
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```

Install Node deps (frontend):

```bash
npm ci    # installs exact lockfile versions
```

---

## 🏃‍♂️  Running locally (no Docker)

Backend – dev server with hot‑reload:

```bash
export FLASK_APP=app.py
export FLASK_ENV=development
cd backend && flask run -p 5050
```

Frontend – rollup dev server:

```bash
npm run dev   # default on http://localhost:5000
```

Edit `src/services/api.js` if you want to hit a remote backend during dev.

---

## 🐳  One‑command deployment (Docker Compose)

```bash
docker compose build
docker compose up -d      # stack available on http://<host>/
```

* `backend` image exposes **5050** inside the compose network.
* `frontend` image (nginx) exposes **80** and proxies `/api/*` to the backend.

Health‑check:

```bash
curl http://<host>/api/health
```

### Update & redeploy

```bash
git pull
# rebuild only the bit you changed
docker compose build backend   # or frontend
docker compose up -d
```

---

## 🔬  Pre‑processing & model training

The Jupyter‑ready script in **backend/sanity.py** walks through:

1. **Raw spectrum → peaks** (binning, normalisation)
2. **Molecule features** (descriptors + fingerprints)
3. **Feature schema analysis** (dynamic, dataset‑wide)
4. **Scaler fitting & feature masks** (saved to `feature_preprocessor.pkl`)
5. **Train / val / test split** (JSONL lines files)
6. **Model fit** (separate notebook / script of your choice)

After training, drop the following artefacts into `backend/models/`:

* `spectrum_predictor.pkl` – the `sklearn` regressor
* `feature_preprocessor.pkl` – feature masks + scalers
* `feature_mapping.json`   – raw vs filtered feature ordering metadata

The runtime API only needs those three files.

---

## 🔑  REST API

All endpoints accept / return **JSON**.

### `POST /api/predict`

```jsonc
{
  "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"
}
```

**Response**

```jsonc
{
  "smiles": "…",
  "molecular_weight": 180.157,
  "exact_mass": 180.0423,
  "spectrum": {       // binned to 1 m/z increments
    "x": [0, 1, 2, …],
    "y": [0.0, 0.0, 0.03, …]
  },
  "peaks": [          // intensity‑sorted list over threshold 0.01
    {"mz": 43.0, "intensity": 0.91},
    {"mz": 77.0, "intensity": 0.72},
    …
  ]
}
```

### `GET /api/health`

```json
{"status":"healthy","model_loaded":true}
```

---

## 🛠  Customisation

* **Change model** – retrain with your own regressor, just keep `.predict(X)` signature and overwrite `models/spectrum_predictor.pkl`.
* **Add more features** – update `molecule_featurizer.py` and re‑run the preprocessing pipeline – the `FeatureSchemaBuilder` will automatically extend the schema.
* **Front‑end tweaks** – Svelte components live in `src/components`; hot‑reload works out of the box (`npm run dev`).

---

## 📜  License

MIT – see `LICENSE` file.

RDKit binaries are licensed under the BSD 3‑Clause; any spectra data you train on may have its own licence – please check before distribution.

---

## 🙏  Acknowledgements

This project started as part of an MSc thesis at QMUL. Thanks to the open‑source authors of **RDKit**, **scikit‑learn**, **Svelte** and **Plotly**.
