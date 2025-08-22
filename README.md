# SMILES2SPEC: AI-Powered Mass Spectrum Prediction System

**MSc AI in Biosciences Dissertation Project**  
**Queen Mary University of London**  
**Author:** Yusuf Mohammed  
**Supervisor:** Mohammed Elbadawi  

<div align="center">

[![MSc Dissertation](https://img.shields.io/badge/MSc%20Dissertation-Queen%20Mary%20University%20of%20London-003E74?style=for-the-badge&logo=graduation-cap)](https://www.qmul.ac.uk/)
[![Status](https://img.shields.io/badge/Status-Production%20Ready-success?style=for-the-badge)](https://github.com/yusufmo1/smiles2spec)
[![Academic Year](https://img.shields.io/badge/Academic%20Year-2024--2025-blue?style=for-the-badge)](https://github.com/yusufmo1)

### Core Technologies

[![Python](https://img.shields.io/badge/Python-3.10--3.11-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![Flask](https://img.shields.io/badge/Flask-3.1.1-000000?style=flat-square&logo=flask&logoColor=white)](https://flask.palletsprojects.com/)
[![SvelteKit](https://img.shields.io/badge/SvelteKit-5.0-FF3E00?style=flat-square&logo=svelte&logoColor=white)](https://kit.svelte.dev/)
[![TypeScript](https://img.shields.io/badge/TypeScript-5.6.3-3178C6?style=flat-square&logo=typescript&logoColor=white)](https://www.typescriptlang.org/)
[![RDKit](https://img.shields.io/badge/RDKit-2025.3.2-07598D?style=flat-square&logo=molecule&logoColor=white)](https://www.rdkit.org/)
[![Plotly](https://img.shields.io/badge/Plotly.js-3.0.1-3F4F75?style=flat-square&logo=plotly&logoColor=white)](https://plotly.com/javascript/)
[![Gunicorn](https://img.shields.io/badge/Gunicorn-21.2.0-499848?style=flat-square&logo=gunicorn&logoColor=white)](https://gunicorn.org/)
[![Docker](https://img.shields.io/badge/Docker-Ready-2496ED?style=flat-square&logo=docker&logoColor=white)](https://www.docker.com/)

### Machine Learning & Performance

[![Machine Learning](https://img.shields.io/badge/ML-Random%20Forest-green?style=flat-square&logo=scikit-learn)](https://scikit-learn.org/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.7.0-EE4C2C?style=flat-square&logo=pytorch&logoColor=white)](https://pytorch.org/)
[![scikit-learn](https://img.shields.io/badge/scikit--learn-1.6.1-F7931E?style=flat-square&logo=scikit-learn&logoColor=white)](https://scikit-learn.org/)
[![Performance](https://img.shields.io/badge/Cosine%20Similarity-0.8138-brightgreen?style=flat-square)](https://github.com/yusufmo1)
[![Samples](https://img.shields.io/badge/Training%20Samples-2,720-blue?style=flat-square)](https://github.com/yusufmo1)
[![Features](https://img.shields.io/badge/Molecular%20Features-2,048+-orange?style=flat-square)](https://github.com/yusufmo1)
[![API Endpoints](https://img.shields.io/badge/API%20Endpoints-20+-purple?style=flat-square)](https://github.com/yusufmo1)
[![Accuracy](https://img.shields.io/badge/R²-0.85+-success?style=flat-square)](https://github.com/yusufmo1)

</div>

---

## Quick Navigation

[Executive Summary](#executive-summary) • 
[System Features](#system-features) • 
[Architecture](#system-architecture) • 
[Installation](#installation-and-setup) • 
[API Documentation](#api-documentation) • 
[Configuration](#configuration) • 
[Technical Specifications](#technical-specifications) • 
[Citation](#citation)

---

## Executive Summary

SMILES2SPEC is a production-ready full-stack application developed as part of an MSc dissertation exploring the integration of artificial intelligence into electronic laboratory notebooks. The system predicts high-resolution electron ionization (EI) mass spectra directly from SMILES molecular notation, combining advanced machine learning models with a robust web-based architecture to provide real-time spectral prediction capabilities for biomedical research applications.

## System Features

### Backend 

* **Flask 2.x REST API** with modular, service-oriented architecture
* **Domain-driven design** with separated concerns:
  - **API Layer**: Routes, middleware, validation schemas
  - **Service Layer**: Business logic and orchestration
  - **Core Layer**: Domain models, ML processors, feature extraction
  - **Integration Layer**: External service abstractions (LLM, databases)
  - **Utils Layer**: Shared utilities and error handling
* **Advanced Molecular Featurization**:
  - 200+ RDKit molecular descriptors
  - 8 different fingerprint types (Morgan, MACCS, RDKit, Avalon, etc.)
  - Electronic properties and structural analysis
  - Parallelised feature extraction for performance
* **Pre-trained Random Forest Model** with 2048+ features
* **AI-powered chat assistant** with spectrum context and chemical knowledge (google/gemma-3-27b-it:free)
* **SMILES generation** from natural language descriptions via OpenRouter LLM integration
* **Comprehensive Export Capabilities**:
  - JSON-formatted spectra with peak lists
  - MSP (Mass Spectral Library) format export
  - Molecular structure images (PNG/SVG)
  - Bulk processing for multiple compounds
* **Production-Ready Features**:
  - **Production WSGI Server**: Gunicorn with optimised worker configuration
  - Pydantic schema validation for all inputs/outputs
  - Global error handling with structured responses
  - Health monitoring and service readiness checks
  - CORS configuration for frontend integration
  - Environment-based configuration management
  - Docker containerization with health checks

### Pre-processing Pipeline

* Parallelised feature extraction (joblib)
* Automatic feature schema generation
* Variance and NaN filtering
* Log-scaling and standardisation
* Support for 2048+ molecular descriptors and fingerprints

### Frontend

* **Svelte/SvelteKit** with TypeScript and modular architecture
* **Performance-Optimised State Management**:
  - Modular store architecture with lazy loading
  - Plot management code loaded on-demand
  - Tree-shaking optimisation for smaller bundles
  - Separate chunks for heavy functionality
* **User-Centred Route Structure**:
  - **Landing Page (`/`)**: Welcome interface with navigation to primary tools
  - **Spectral Simulation (`/spectral-simulation`)**: Full prediction interface with interactive visualisation
  - **How It Works (`/how-it-works`)**: Educational content about mass spectrum prediction
  - **About (`/about`)**: Developer information and project details
  - **Chat with Spectrum (`/chat-with-spectrum`)**: AI-powered spectral analysis assistance
* **Real-time spectrum visualisation** with Plotly.js
* **Interactive molecular structure display**
* **Bulk SMILES processing** with progress tracking
* **Export capabilities** (JSON, MSP, CSV) with format validation
* **Responsive design** with modern UI components and accessibility features

## System Architecture

### Backend Architecture 

```
backend/
├── main.py                    # Application entry point
├── config/                    # Centralized configuration
│   └── settings.py           # Environment-based config
├── api/                      # API layer
│   ├── routes/              # Domain-specific route modules
│   │   ├── prediction.py    # Spectrum prediction endpoints
│   │   ├── export.py        # File export endpoints
│   │   ├── chat.py          # AI chat endpoints
│   │   ├── health.py        # Health check endpoints
│   │   └── upload.py        # File upload endpoints
│   ├── schemas/             # Pydantic validation schemas
│   │   ├── prediction.py    # Prediction request/response schemas
│   │   └── chat.py          # Chat message schemas
│   └── middleware/          # Request/response middleware
│       ├── validation.py    # Input validation decorators
│       └── error_handler.py # Global error handling
├── services/                # Business logic layer
│   └── prediction_service.py # Spectrum prediction orchestration
├── core/                    # Core domain layer
│   ├── models/             # Domain models
│   │   ├── molecule.py     # Molecular data structures
│   │   └── spectrum.py     # Spectrum data structures
│   ├── processors/         # Core processing logic
│   │   ├── feature_processor.py    # Molecular feature extraction
│   │   ├── spectrum_processor.py   # Spectrum processing
│   │   └── feature_preprocessor.py # Feature preprocessing
│   └── ml/                 # Machine learning components
│       └── model_handler.py # ML model operations
├── integrations/           # External service integrations
│   └── llm/               # LLM integration
│       ├── client.py      # LLM API client
│       └── services/      # LLM service implementations
│           ├── chat_service.py    # Chat with spectrum context
│           └── smiles_service.py  # SMILES generation
├── utils/                 # Shared utilities
│   ├── errors.py         # Exception hierarchy
│   ├── logging.py        # Centralized logging
│   ├── chemistry.py      # Chemical utilities
│   ├── data.py           # Data conversion utilities
│   └── formats.py        # File format utilities
└── models/                # ML model files (not in git)
    ├── spectrum_predictor.pkl     # Trained Random Forest model
    ├── feature_preprocessor.pkl   # Feature scaling pipeline
    └── feature_mapping.json       # Feature index mapping
```

### Frontend Architecture

```
frontend/
├── src/
│   ├── routes/                    # SvelteKit pages
│   │   ├── +page.svelte          # Landing page with project overview
│   │   ├── +layout.svelte        # Root layout
│   │   ├── spectral-simulation/  # Main prediction interface
│   │   ├── about/                # About page
│   │   ├── chat-with-spectrum/   # AI chat interface
│   │   └── how-it-works/         # Documentation
│   └── lib/                      # Core library
│       ├── components/           # Reusable components
│       │   ├── landing/          # Landing page components
│       │   ├── smiles-input/     # SMILES input system
│       │   ├── panels/           # UI panel system
│       │   │   ├── simulation/   # Spectrum visualisation
│       │   │   ├── info/         # Information panels
│       │   │   └── about/        # About panels
│       │   ├── chat/             # Chat components
│       │   └── icons/            # SVG icon library
│       ├── services/             # API and business logic
│       │   ├── api.ts           # Backend API client
│       │   ├── plotlyService.ts # Plotting utilities
│       │   └── chatService.ts   # Chat functionality
│       ├── stores/              # Modular state management
│       │   ├── index.ts         # Main store exports with lazy loading
│       │   ├── appState.ts      # Global application state
│       │   ├── panelStore.ts    # Panel definitions and management
│       │   ├── carouselStore.ts # Panel navigation and carousel
│       │   ├── pageStore.ts     # Page routing and navigation
│       │   ├── plotEffects.ts   # Lazy-loaded plot management
│       │   └── types.ts         # TypeScript type definitions
│       └── styles/              # Global styles
│           ├── theme.css        # CSS variables
│           └── tokens.css       # Design tokens
├── static/                      # Static assets
├── package.json                 # Dependencies
└── vite.config.js              # Build configuration
```

### Key Architectural Principles

1. **Separation of Concerns**: Clear boundaries between API, business logic, and data layers
2. **Dependency Injection**: Services are injected rather than tightly coupled
3. **Error Handling**: Structured exception hierarchy with global error handlers
4. **Validation**: Input validation at API boundaries using Pydantic schemas
5. **Testability**: Modular design enables comprehensive unit and integration testing
6. **Maintainability**: Domain-driven organisation makes code easy to understand and modify
7. **Performance**: Lazy loading and code splitting optimise bundle size and initial load times

### Backend Components

#### Core Services
- **`PredictionService`**: Orchestrates spectrum prediction from SMILES strings
- **`ChatService`**: AI-powered chat functionality with spectrum context integration
- **`SMILESService`**: Natural language to SMILES generation using LLMs

#### Machine Learning Pipeline
- **`ModelHandler`**: Manages trained model loading and inference operations
- **`FeatureProcessor`**: Extracts 2048+ molecular features using RDKit
- **`FeaturePreprocessor`**: Applies scaling, filtering, and standardisation
- **`SpectrumProcessor`**: Converts model outputs to formatted spectra and peak lists

#### Data Processing & Validation
- **Pydantic Schemas**: Comprehensive request/response validation
- **Error Handling**: Custom exception hierarchy with structured error responses
- **File Processing**: Support for CSV/TXT bulk uploads and MSP exports
- **Image Generation**: Molecular structure visualisation (PNG/SVG)

#### API Architecture
- **Route Blueprints**: Domain-separated endpoint organisation
- **Middleware**: Global validation, error handling, and CORS configuration
- **Health Monitoring**: Service readiness and model status checking

#### Frontend Performance Architecture
- **Modular Store System**: Separated stores for different concerns (app state, panels, navigation)
- **Lazy Loading**: Plot management functionality loaded on-demand
- **Code Splitting**: Heavy functionality in separate chunks for optimal bundle size
- **Tree Shaking**: Unused store logic excluded from production bundles
- **Memory Management**: Automatic cleanup of plot resources and event listeners

## Installation and Setup

### Prerequisites

* Python 3.10-3.11 (Python 3.12+ has compatibility issues with RDKit)
* Node.js 18+ with npm
* RDKit molecular toolkit (installation via conda-forge required)
* Docker and Docker Compose (optional for containerised deployment)

### Backend Setup

#### Development Mode
```bash
cd backend
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install -r requirements.txt

# Optional: Set up environment variables
echo "OPENROUTER_API_KEY=your_api_key_here" > .env

# Start development server
python main.py
```

#### Production Mode (WSGI Server)
```bash
cd backend
# Install dependencies including Gunicorn
pip install -r requirements.txt

# Run with Gunicorn (production WSGI server)
gunicorn --config gunicorn.conf.py backend.wsgi:app

# Alternative: Run from parent directory
cd ..
gunicorn --config backend/gunicorn.conf.py backend.wsgi:app
```

#### Docker Deployment
```bash
# Full stack deployment with production WSGI server
docker-compose up --build
# Backend: http://localhost:5050 (Production WSGI with Gunicorn)
# Frontend: http://localhost:3001
```

The API will be available at `http://localhost:5050`

### Frontend Setup

```bash
cd frontend
npm install
npm run dev
```

The frontend will be available at `http://localhost:5173`

## API Documentation

### Spectrum Prediction
- `POST /predict` - Predict mass spectrum from SMILES string
- `GET /structure` - Get molecular structure as base64 PNG
- `GET /structure/png` - Alternative endpoint for molecular structure

### File Operations & Export
- `POST /smiles_bulk` - Upload CSV/TXT files with bulk SMILES data
- `POST /export_msp` - Export single spectrum in MSP format
- `POST /export_msp_batch` - Export multiple spectra as combined MSP file

### AI-Powered Features
- `POST /chat` - Interactive chat with AI assistant (includes spectrum context)
- `POST /generate_smiles` - Generate SMILES from natural language descriptions

### System Health & Monitoring
- `GET /health` - Comprehensive service health check with model status

### Example Usage

```bash
# Predict spectrum with detailed response
curl -X POST http://localhost:5050/predict \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"}'

# Response includes:
# - Chemical name and molecular properties
# - High-resolution spectrum data (x/y arrays)
# - Peak list with m/z and intensity values
# - Base64-encoded molecular structure PNG
# - Metadata and processing information

# Get molecular structure image
curl "http://localhost:5050/structure?smiles=CCO"

# Upload bulk SMILES file (CSV or TXT)
curl -X POST http://localhost:5050/smiles_bulk \
  -F "file=@molecules.csv"

# Export spectrum in MSP format (mass spectral library format)
curl -X POST http://localhost:5050/export_msp \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO"}' \
  --output ethanol_spectrum.msp

# Chat with AI about spectrum interpretation
curl -X POST http://localhost:5050/chat \
  -H "Content-Type: application/json" \
  -d '{
    "messages": [{"role": "user", "content": "Explain the fragmentation pattern"}],
    "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "stream": false
  }'

# Generate SMILES from description
curl -X POST http://localhost:5050/generate_smiles \
  -H "Content-Type: application/json" \
  -d '{"description": "anti-inflammatory drug with benzene ring"}'
```


## Configuration

### Environment Variables
- `API_HOST` - Server host (default: 0.0.0.0)
- `API_PORT` - Server port (default: 5050)
- `API_DEBUG` - Debug mode (default: true)
- `OPENROUTER_API_KEY` - OpenRouter API key for LLM features (optional)
- `LOG_LEVEL` - Logging level (default: INFO)

### Backend Configuration

#### API Settings (`config/settings.py`)
- **host/port**: Server binding configuration
- **debug**: Development mode settings
- **model_path**: Path to trained spectrum predictor model
- **preprocessor_path**: Feature preprocessing pipeline location
- **feature_mapping_path**: Molecular feature mapping definitions

#### Production WSGI Configuration (`gunicorn.conf.py`)
- **Workers**: Optimised for ML workload (2 workers default)
- **Worker Class**: Sync workers for RDKit/ML compatibility
- **Timeouts**: Extended timeout (300s) for ML inference operations
- **Memory Management**: Preload app for efficient model loading
- **Logging**: Structured production logging with access logs
- **Health Checks**: Configurable health monitoring endpoints

#### Feature Extraction Configuration
- **Molecular Descriptors**: 200+ RDKit descriptors (MW, LogP, TPSA, etc.)
- **Fingerprints**: Multiple fingerprint types with configurable parameters:
  - Morgan fingerprints (radii 1-3, size 1024)
  - Morgan feature fingerprints (radius 2, size 1024)
  - MACCS keys (166 structural keys)
  - Topological fingerprints (size 1024)
  - RDKit fingerprints (size 2048)
  - Avalon fingerprints (size 1024)
  - Pattern fingerprints (size 1024)
  - Layered fingerprints (size 2048)
- **Count Features**: Bond counts, atom counts, ring analysis
- **Electronic Properties**: HOMO/LUMO estimation, dipole moments

#### Spectral Processing
- **m/z Range**: Configurable mass-to-charge ratio ranges
- **Resolution**: Spectral resolution and binning parameters
- **Peak Detection**: Intensity thresholds and peak finding algorithms

#### LLM Integration
- **Provider**: OpenRouter gateway to multiple LLM providers
- **Model**: google/gemma-3-27b-it:free (primary) for spectrum interpretation and SMILES generation
- **Temperature**: Response creativity (0.7 default)
- **Max Tokens**: Response length limits (1000 default)

## Technical Specifications

### Machine Learning Implementation

The system employs a sophisticated machine learning pipeline developed through extensive experimentation documented in the SMILES2SPEC Foundry research pipeline:

* **Algorithm**: Random Forest Regression optimised for spectral prediction
* **Features**: 2048+ molecular descriptors and fingerprints including:
  - RDKit molecular descriptors (molecular weight, LogP, TPSA, etc.)
  - Morgan fingerprints (multiple radii for different structural patterns)
  - MACCS keys (166 structural keys for pharmacophore analysis)
  - Topological, Avalon, Pattern, and Layered fingerprints
  - Bond/atom counts and electronic properties
* **Training Data**: 2,720 curated electron ionisation mass spectrometry samples
* **Performance**: Cosine similarity of 0.8138 achieved through Bayesian optimisation
* **Preprocessing Pipeline**:
  - Variance and NaN filtering for feature selection
  - Log-scaling and standardisation for optimal model performance
  - Automatic feature schema generation and validation
  - Parallelised processing for high-throughput applications

### Technical Capabilities
* **Request Validation**: Comprehensive Pydantic schemas for all endpoints
* **Error Handling**: Structured exception hierarchy with meaningful error messages
* **File Processing**: Support for CSV/TXT bulk uploads with intelligent parsing
* **Export Formats**: MSP (mass spectral library), JSON, and image formats
* **AI Integration**: OpenRouter/google/gemma-3-27b-it:free for spectrum interpretation and SMILES generation
* **Performance**: Optimised for both single predictions and bulk processing
* **Monitoring**: Health checks with model status and service readiness indicators
* **Security**: Input sanitisation, file upload validation, and environment-based secrets


## Academic Context

This system was developed as one of two complementary implementations demonstrating AI integration into electronic laboratory notebooks. Together with the GUARDIAN pharmaceutical compliance system, it forms the technical foundation of the dissertation "Integrating AI into Electronic Lab Notebooks" submitted for the MSc AI in Biosciences programme at Queen Mary University of London.

## Citation

For academic use of this work, please cite:

```bibtex
@mastersthesis{mohammed2025smiles2spec,
  title = {SMILES2SPEC: AI-Powered Mass Spectrum Prediction System for Electronic Laboratory Notebooks},
  author = {Mohammed, Yusuf},
  year = {2025},
  school = {Queen Mary University of London},
  department = {MSc AI in Biosciences},
  supervisor = {Elbadawi, Mohammed},
  note = {MSc Dissertation Project: Integrating AI into Electronic Lab Notebooks}
}
```

---

<div align="center">

[![Author](https://img.shields.io/badge/Author-Yusuf%20Mohammed-blue?style=flat-square&logo=github)](https://github.com/yusufmo1)
[![Supervisor](https://img.shields.io/badge/Supervisor-Dr%20Mohammed%20Elbadawi-green?style=flat-square&logo=github)](https://github.com/Dr-M-ELBA)
[![Institution](https://img.shields.io/badge/Institution-QMUL-003E74?style=flat-square)](https://www.qmul.ac.uk/)
[![Programme](https://img.shields.io/badge/Programme-MSc%20AI%20in%20Biosciences-purple?style=flat-square)](https://www.qmul.ac.uk/)

</div>

*Developed as part of MSc AI in Biosciences dissertation at Queen Mary University of London (2025)*
