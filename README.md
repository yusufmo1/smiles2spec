# SMILES2SPEC

**Predict high‑resolution EI mass spectra directly from SMILES** – a full‑stack application that couples a trained machine‑learning model (scikit‑learn + RDKit) to a lightweight Flask API and a Svelte single‑page UI.

---

## ✨  Features

| Layer                       | Highlights                                                                                                                                                                                           |
| --------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Backend**                 | Flask 2.x REST API & Gunicorn server, on‑the‑fly RDKit featurisation, pretrained scikit‑learn regressor, JSON output with spectrum, peak list, molecular structure SVG and MSP export.               |
| **Pre‑processing pipeline** | Parallelised (joblib) spectrum featuriser, automatic feature schema builder, variance / NaN filtering, StandardScaler wrapper, artefact writer (`feature_preprocessor.pkl`, `feature_mapping.json`). |
| **Model**                   | Ready‑to‑use RandomForest intensity regressor (`spectrum_predictor.pkl`).                                                                                                                            |
| **Frontend**                | Svelte + Plotly.js for interactive spectrum plots, chemical‑structure panel rendered from SVG, bulk SMILES upload (TXT), arrow navigator, MSP download, fully static bundle served by Nginx.         |

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
└── README.md
```

---

## ⚙️  Requirements

| Tool   | Tested version |
| ------ | -------------- |
| Python | 3.10           |
| Node   | 20.x LTS       |
| RDKit  | 2023.09        |
| npm    | 10.x           |

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

## 🏃‍♂️  Running locally

Backend – dev server with hot‑reload:

```bash
export FLASK_APP=app.py
export FLASK_ENV=development
cd backend && python3 -m backend.app
```

Frontend – rollup dev server:

```bash
npm run dev   # default on http://localhost:5000
```

Edit `src/services/api.js` if you want to hit a remote backend during dev.

---

## 🔬  Model provenance & feature pipeline

SMILES in → **RDKit featurisation** (physicochemical descriptors ≈ 8,000 raw features) → variance / NaN filtering → log‑scaling & `StandardScaler` → **RandomForest intensity regressor**.

The resulting artefacts (`feature_preprocessor.pkl`, `feature_mapping.json`, `spectrum_predictor.pkl`) live in `backend/models/` – you **don't need to retrain**.

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

RDKit binaries are licensed under the BSD 3‑Clause; any spectra data you train on may have its own licence – please check before distribution.

---

## 🙏  Acknowledgements

Built with **RDKit**, **scikit‑learn**, **Svelte** and **Plotly.js**. Many thanks to their respective communities.

## Chat with Spectrum

The Chat with Spectrum feature provides an interactive chat interface integrated directly into the main application interface to help with:

- Understanding mass spectrometry concepts
- Interpreting SMILES notation
- Learning about chemical structures
- Getting explanations about spectral data

### How It Works

The chat interface connects to a language model backend that has been specialized in chemistry and mass spectrometry concepts. Users can ask questions about:

1. General chemistry concepts
2. SMILES notation interpretation
3. Mass spectrometry principles
4. Interpretation of spectral data
5. Chemical compound information

### Using the Chat

1. Find the "CHAT WITH SPECTRUM" panel on the main interface
2. Type your question in the input box at the bottom of the panel
3. Press Enter or click the send button
4. Receive a response from Spectrum directly within the panel

### Future Enhancements

Planned enhancements for the Chat with Spectrum feature include:

- Integration with your current spectrum data for contextualized answers
- File upload capabilities for spectral data interpretation
- Export of chat conversations for documentation
- Multi-language support
