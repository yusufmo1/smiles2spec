# SPEC2SMILES

**Predict high‑resolution EI mass spectra directly from SMILES** – a full‑stack application that couples a trained machine‑learning model (scikit‑learn + RDKit) to a lightweight Flask API and a Svelte single‑page UI.

---

## ✨  Features

| Layer                       | Highlights                                                                                                                                                                                           |
| --------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Backend**                 | Flask 2.x REST API & Gunicorn server, on‑the‑fly RDKit featurisation, pretrained scikit‑learn regressor, JSON output with spectrum, peak list, molecular structure SVG and MSP export.               |
| **Pre‑processing pipeline** | Parallelised (joblib) spectrum featuriser, automatic feature schema builder, variance / NaN filtering, StandardScaler wrapper, artefact writer (`feature_preprocessor.pkl`, `feature_mapping.json`). |
| **Model**                   | Ready‑to‑use RandomForest intensity regressor (`spectrum_predictor.pkl`).                                                                                                                            |
| **Frontend**                | Svelte + Plotly.js for interactive spectrum plots, chemical‑structure panel rendered from SVG, bulk SMILES upload (TXT), arrow navigator, MSP download, fully static bundle served by Nginx.         |
| **Dev Ops**                 | Multi‑arch Docker build (arm64 & amd64), two‑service `docker‑compose`, automatic proxying `/api/*` → backend.                                                                                        |

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

## 🔬  Model provenance & feature pipeline

SMILES in → **RDKit featurisation** (physicochemical descriptors + 7 fingerprint families, **≈ 8 000 raw features**) → variance / NaN filtering → log‑scaling & `StandardScaler` → **RandomForest intensity regressor**.

The resulting artefacts (`feature_preprocessor.pkl`, `feature_mapping.json`, `spectrum_predictor.pkl`) live in `backend/models/` – you **don’t need to retrain**.

---

## 🔑  REST API

All endpoints accept / return **JSON**.

### `POST /api/predict`

Request

```jsonc
{
  "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"
}
```

Response

```jsonc
{
  "smiles": "…",
  "chemical_name": "acetylsalicylic acid",
  "molecular_weight": 180.157,
  "exact_mass": 180.0423,
  "spectrum": {
    "x": [0, 1, 2, …],
    "y": [0.0, 0.0, 0.03, …]
  },
  "peaks": [
    {"mz": 43.0, "intensity": 0.91},
    {"mz": 77.0, "intensity": 0.72}
  ],
  "structure_svg": "<svg…/svg>"
}
```

### `POST /api/export_msp`

```jsonc
{
  "smiles": "CCO"
}
```

Returns a downloadable `.msp` file (intensities scaled to 1000).

### `POST /api/smiles_bulk`

Multipart upload of a `.txt` file (one SMILES per line) → returns JSON list for the front‑end navigator.

### `GET /api/health`

```json
{"status":"healthy","model_loaded":true}
```

---

## 📜  License

MIT – see `LICENSE` file.

RDKit binaries are licensed under the BSD 3‑Clause; any spectra data you train on may have its own licence – please check before distribution.

---

## 🙏  Acknowledgements

Built with **RDKit**, **scikit‑learn**, **Svelte** and **Plotly.js**. Many thanks to their respective communities.
