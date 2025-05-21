# SMILES2SPEC

SMILES2SPEC is a full-stack application that predicts high-resolution electron ionization (EI) mass spectra directly from SMILES strings. It combines a trained machine-learning model (scikit-learn and RDKit) with a lightweight Flask backend and an interactive Svelte frontend.

## ✨ Features

### Backend

* **Flask 2.x REST API** with Gunicorn server
* Real-time molecule featurization using **RDKit**
* Pre-trained scikit-learn regression model
* Outputs JSON-formatted spectra, peak lists, molecular structures (SVG/PNG), and MSP file exports

### Pre-processing Pipeline

* Parallelized feature extraction (joblib)
* Automatic feature schema generation
* Variance and NaN filtering
* Log-scaling and standardization (StandardScaler)
* Preprocessing artifacts saved for consistent deployment

### Model

* Ready-to-use **RandomForest** regressor for spectrum intensity prediction

### Frontend

* Interactive spectrum visualization with **Plotly.js**
* Chemical structures rendered directly from SVG/PNG
* Supports bulk upload of SMILES strings (TXT/CSV)
* Easy navigation between molecules
* Direct MSP file downloads
* Fully static frontend deployment via Nginx

### AI Integration

* **Chat interface** powered by a specialized language model for spectral analysis and chemical explanations
* AI-assisted generation of chemically valid SMILES strings from natural language descriptions

## 📁 Repository Layout

```
.
├── backend/                  # Flask API + ML artifacts
│   ├── app.py                # Entry‑point for the Flask application
│   ├── models/               # Pretrained models and preprocessing artifacts
│   ├── llm_integration/      # AI modules for chat and SMILES generation
│   ├── templates/            # Jinja templates for fallback pages
│   └── ...                   # Additional source modules
├── frontend/                 # Svelte application
│   ├── src/                  # Svelte source code
│   │   ├── components/       # UI components
│   │   └── services/         # API interactions and utilities
│   ├── public/               # Static assets
│   └── ...                   # Configuration files
├── Dockerfiles               # Docker deployment configurations
└── README.md
```

## ⚙️ Requirements

| Tool   | Tested version |
| ------ | -------------- |
| Python | 3.10           |
| Node   | 20.x LTS       |
| RDKit  | 2023.09+       |
| npm    | 10.x           |

### Installation

**Backend:**

```bash
cd backend
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```

**Frontend:**

```bash
npm ci  # Installs exact lockfile versions
```

## 🏃‍♂️ Running Locally

**Backend (with hot-reload):**

```bash
export FLASK_APP=app.py
export FLASK_ENV=development
cd backend && python3 -m backend.app
```

**Frontend:**

```bash
npm run dev  # Available at http://localhost:5000
```

For remote backend during development, adjust `src/services/api.js`.

## 🐳 Docker Deployment

Dockerfiles are provided for easy deployment:

```bash
# Build
 docker build -f backend/Dockerfile -t smiles2spec-backend .
 docker build -f frontend/Dockerfile -t smiles2spec-frontend .

# Run
 docker run -p 5050:5050 smiles2spec-backend
 docker run -p 80:80 smiles2spec-frontend
```

## 🔬 Model and Feature Pipeline

* Input SMILES → RDKit featurization (approx. 8,000 features) → variance/NaN filtering → scaling and standardization → RandomForest regression model
* Model and preprocessing artifacts are stored in `backend/models/`.
* Production deployment supports downloading models from Azure Blob Storage.

## 🔑 REST API Endpoints

* **`POST /api/predict`**: Predict spectra from SMILES strings
* **`POST /api/export_msp`**: Generate downloadable MSP file
* **`POST /api/export_msp_batch`**: Batch MSP file export
* **`POST /api/smiles_bulk`**: Bulk SMILES upload via TXT or CSV
* **`POST /api/chat`**: AI-powered chat interactions
* **`POST /api/generate_smiles`**: Generate SMILES from descriptions
* **`GET /api/structure` and `GET /api/structure/png`**: Retrieve molecular structures
* **`GET /api/health`**: Check API health status

## 🧠 Chat with Spectrum

An integrated chat assistant to support chemical analyses:

* Detailed explanations of SMILES notation and mass spectrometry
* Spectrum interpretation and compound analysis
* Fully integrated within the main user interface
* Context-aware chemical insights

## 🪄 AI SMILES Generation

Generate chemically valid SMILES using AI:

* Quickly produce multiple drug-like molecules
* Guide molecule generation with specific descriptions
* Ideal for batch analysis and predictive workflows

## 🖼️ User Interface

Modern and intuitive UI designed for chemical analysis:

* Interactive Plotly.js visualizations
* Chemical structure rendering and fragment ion analysis
* Detailed peak data and fragment identifications
* Bulk analysis capabilities
* Flexible export options (MSP, PNG)
* Integrated console output for detailed analyses

## 📜 License

Licensed under MIT (see `LICENSE` file).
RDKit binaries are licensed under BSD 3-Clause. Ensure compliance with licensing for any spectral data used in training.

## 🙏 Acknowledgements

Built with RDKit, scikit-learn, Svelte, and Plotly.js. Thanks to their respective communities.
