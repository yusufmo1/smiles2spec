# Mass Spectrometry Prediction API

This is a Flask-based API for predicting mass spectrometry spectra from SMILES strings.

## Setup

### Prerequisites
- Python 3.8+
- RDKit
- Flask
- Flask-CORS
- NumPy
- Other dependencies listed in requirements.txt

### Installation

1. Create a virtual environment:
```bash
python -m venv msvenv
source msvenv/bin/activate  # On Windows: msvenv\Scripts\activate
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Running the API

Start the API server:
```bash
python app.py
```

The server will run on `http://0.0.0.0:8080` by default (configurable in `config.py`).

## API Endpoints

### Predict Spectrum
- **URL**: `/api/predict`
- **Method**: `POST`
- **Content-Type**: `application/json`
- **Request Body**:
  ```json
  {
    "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"
  }
  ```
- **Response**:
  ```json
  {
    "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "molecular_weight": 180.159,
    "exact_mass": 180.042259,
    "spectrum": {
      "x": [0, 1, 2, ..., 499],
      "y": [0, 0.01, 0.55, ..., 0]
    },
    "peaks": [
      {"mz": 120.0, "intensity": 0.85},
      {"mz": 138.0, "intensity": 0.62},
      ...
    ]
  }
  ```
- **Error Response**:
  ```json
  {
    "error": "Invalid SMILES string: X"
  }
  ```

### Health Check
- **URL**: `/api/health`
- **Method**: `GET`
- **Response**:
  ```json
  {
    "status": "healthy",
    "model_loaded": true
  }
  ```

## Configuration

Configuration settings are stored in `config.py`:

- **API_CONFIG**: Contains server settings (host, port, debug mode) and model paths
- **FEATURE_CONFIG**: Controls molecular feature extraction
- **SPECTRAL_CONFIG**: Defines parameters for spectral processing
- **FEATURE_PROCESSING_CONFIG**: Configures feature preprocessing

The default port is `8080` and the host is `0.0.0.0` (accessible from any network interface).

## Web Interface

A simple web interface is available at the root URL (`/`) for testing predictions directly.

## Directory Structure
- `app.py`: Main Flask application and API endpoints
- `prediction_service.py`: Core prediction functionality
- `model_handler.py`: Handles model loading and prediction
- `molecule_featurizer.py`: Extracts molecular features from SMILES
- `spectrum_processor.py`: Processes spectral data
- `feature_preprocessor.py`: Handles feature preprocessing
- `utils.py`: Utility functions
- `config.py`: Configuration settings
- `models/`: Directory for model files
