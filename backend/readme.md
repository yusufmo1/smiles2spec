# SMILES2SPEC Backend

The backend service for SMILES2SPEC - a mass spectrometry prediction system that generates high-resolution electron ionization (EI) mass spectra from SMILES molecular structures.

## Overview

This Flask-based backend provides:
- RESTful API for spectrum prediction
- Production WSGI server with Gunicorn
- Machine learning model serving
- Molecular feature extraction with RDKit
- Real-time chat with AI assistance
- Data export in multiple formats
- Docker containerization with health checks

## Architecture

```
backend/
├── api/                 # API layer - routes, middleware, schemas
├── core/                # Core business logic
├── services/            # High-level service orchestration
├── integrations/        # External service integrations
├── models/              # Trained ML models (not in git)
├── config/              # Configuration management
├── utils/               # Shared utilities
├── main.py              # Application entry point (development)
├── wsgi.py              # WSGI entry point (production)
├── gunicorn.conf.py     # Gunicorn production configuration
├── Dockerfile           # Docker container definition
└── requirements.txt     # Python dependencies including Gunicorn
```

## Quick Start

### Prerequisites

- Python 3.8+
- RDKit (via conda recommended)
- 4GB+ RAM for model loading

### Installation

1. Create and activate a virtual environment:
```bash
cd backend
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Set up environment variables:
```bash
# Create .env file in backend directory
echo "OPENROUTER_API_KEY=your_api_key_here" > .env
```

4. Download model files:
```bash
# Place these files in backend/models/:
# - spectrum_predictor.pkl
# - feature_preprocessor.pkl
# - feature_mapping.json
```

### Running the Server

#### Development Mode
```bash
python main.py
```

#### Production Mode (WSGI Server)
```bash
# Using Gunicorn with optimized configuration
gunicorn --config gunicorn.conf.py backend.wsgi:app

# Alternative: Simple Gunicorn command
gunicorn -w 2 -b 0.0.0.0:5050 --timeout 300 backend.wsgi:app
```

#### Docker Deployment
```bash
# Build and run with Docker Compose
docker-compose up --build

# Or build and run manually
docker build -t smiles2spec-backend .
docker run -p 5050:5050 smiles2spec-backend
```

The API will be available at `http://localhost:5050`

## API Endpoints

### Core Endpoints

- `POST /predict` - Predict mass spectrum from SMILES
- `POST /chat` - Chat with AI about spectra (streaming)
- `GET /health` - Service health check

### Export Endpoints

- `POST /export/msp` - Export single spectrum as MSP
- `POST /export/csv` - Export spectrum data as CSV
- `POST /export/json` - Export spectrum data as JSON

### Utility Endpoints

- `POST /generate_smiles` - Generate SMILES from description
- `POST /upload` - Upload file with multiple SMILES
- `GET /structure` - Get molecular structure as PNG

## Project Structure

### `/api` - API Layer

Handles HTTP requests and responses:
- **routes/** - Endpoint definitions
- **middleware/** - Error handling, validation
- **schemas/** - Request/response validation

### `/core` - Core Business Logic

Contains the ML and chemistry components:
- **ml/** - Model loading and prediction
- **processors/** - Feature extraction and spectrum processing
- **models/** - Data models (Pydantic)

### `/services` - Service Layer

High-level orchestration:
- **prediction_service.py** - Main prediction workflow

### `/integrations` - External Services

Third-party integrations:
- **llm/** - OpenRouter/LLM integration for chat (google/gemma-3-27b-it:free)

### `/config` - Configuration

- **settings.py** - Centralized configuration with environment overrides

### `/utils` - Utilities

Shared functionality:
- **chemistry.py** - SMILES validation, molecular operations
- **formats.py** - Data format conversions (MSP, JCAMP)
- **errors.py** - Custom exception classes
- **logging.py** - Logging configuration

## Configuration

Configuration is managed through:
1. Default values in `config/settings.py`
2. Environment variables in `.env`
3. Runtime environment variables

Key configuration options:

```python
# API Configuration
API_HOST=0.0.0.0
API_PORT=5050
API_DEBUG=true

# LLM Configuration
OPENROUTER_API_KEY=your_key_here

# Model Paths (relative to backend/)
MODEL_PATH=models/spectrum_predictor.pkl
PREPROCESSOR_PATH=models/feature_preprocessor.pkl
FEATURE_MAPPING_PATH=models/feature_mapping.json

# Production WSGI Configuration (gunicorn.conf.py)
GUNICORN_WORKERS=2
GUNICORN_TIMEOUT=300
GUNICORN_LOG_LEVEL=info
```

### Production Configuration

The application includes a production-ready Gunicorn configuration in `gunicorn.conf.py`:

- **Workers**: 2 sync workers (optimized for ML workloads)
- **Timeout**: 300 seconds (extended for ML inference)
- **Preload App**: Enabled for efficient model loading
- **Health Checks**: Integrated Docker health monitoring
- **Logging**: Structured production logging
- **Memory Management**: Optimized for RDKit/ML operations

## Machine Learning Pipeline

### Feature Extraction

The system extracts 2048+ molecular features including:
- RDKit molecular descriptors (200+)
- Multiple fingerprint types (Morgan, MACCS, etc.)
- Electronic features (partial charges)
- Structural counts (atoms, bonds)

### Model Architecture

- **Algorithm**: Random Forest Regressor
- **Input**: 2048+ molecular features
- **Output**: 500 intensity values (0-499 m/z)
- **Training**: ~280,000 NIST spectra

### Prediction Flow

1. SMILES validation and parsing
2. Molecular feature extraction
3. Feature preprocessing (scaling, normalization)
4. Model prediction
5. Spectrum post-processing
6. Response formatting

## Development

### Code Style

- Follow PEP 8 guidelines
- Use type hints for function signatures
- Add docstrings to all public functions
- Keep functions focused and testable

### Testing

Run tests:
```bash
pytest tests/
```

### Debugging

Enable debug logging:
```python
# In your code
from utils import logger
logger.debug("Debug message")
```

### Common Issues

1. **RDKit Import Error**: Install RDKit via conda
2. **Model Not Found**: Ensure model files are in `backend/models/`
3. **Memory Error**: Increase available RAM or reduce batch size
4. **CORS Error**: Check allowed origins in `main.py`

## API Examples

### Predict Spectrum

```bash
curl -X POST http://localhost:5050/predict \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO"}'
```

### Chat with AI

```bash
curl -X POST http://localhost:5050/chat \
  -H "Content-Type: application/json" \
  -d '{
    "messages": [{"role": "user", "content": "What is this molecule?"}],
    "smiles": "CCO",
    "stream": true
  }'
```

### Export MSP

```bash
curl -X POST http://localhost:5050/export/msp \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO"}' \
  --output ethanol.msp
```

### Health Check

```bash
curl http://localhost:5050/health
```

## Performance

### Development Mode
- Prediction time: ~100-500ms per molecule
- Memory usage: ~1-2GB with models loaded
- Single-threaded Flask development server

### Production Mode (Gunicorn)
- Prediction time: ~100-500ms per molecule
- Memory usage: ~2-4GB with 2 workers and preloaded models
- Concurrent requests: 2 sync workers optimized for ML workloads
- Health monitoring: Docker health checks every 30s
- Graceful shutdown: 30s timeout for request completion

## Security

- Input validation on all endpoints
- CORS configured for specific origins
- API key required for LLM features
- No user data persistence

## Contributing

1. Create a feature branch
2. Add tests for new functionality
3. Ensure all tests pass
4. Submit a pull request

## License

This project is part of the SMILES2SPEC application.