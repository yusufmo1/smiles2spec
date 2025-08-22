# Services Layer

The services layer provides high-level orchestration of business logic, coordinating between the API layer and core components.

## Overview

Services act as the main entry points for complex operations, handling:
- Component initialization
- Workflow orchestration
- Error handling and recovery
- State management

## PredictionService

The main service for spectrum prediction operations.

### Usage

```python
from services.prediction_service import PredictionService

# Initialize service (happens automatically on first use)
service = PredictionService()
service.initialize()

# Single prediction
result = service.predict_spectrum("CCO")

# Batch prediction
results = service.predict_batch(["CCO", "CC(=O)O", "c1ccccc1"])

# Check service status
if service.is_ready():
    print("Service ready for predictions")
```

### Architecture

```
PredictionService
    ├── ModelHandler (ML inference)
    ├── FeatureProcessor (feature extraction)
    └── SpectrumProcessor (data formatting)
```

### Workflow

1. **Initialization**
   - Load ML models from disk
   - Initialize processors
   - Validate configuration

2. **Prediction Flow**
   ```
   SMILES Input
        ↓
   Validation & Cleaning
        ↓
   Feature Extraction
        ↓
   Model Prediction
        ↓
   Spectrum Generation
        ↓
   Result Formatting
   ```

3. **Error Handling**
   - Invalid SMILES → Clear error message
   - Model not loaded → Graceful degradation
   - Feature extraction failure → Detailed logging

### Key Methods

#### `initialize() -> bool`
Loads all required models and processors. Returns success status.

#### `predict_spectrum(smiles: str) -> PredictionResult`
Main prediction method. Returns complete prediction data including:
- Predicted spectrum
- Molecular properties
- Structure visualization
- Peak list

#### `predict_batch(smiles_list: List[str]) -> List[PredictionResult]`
Batch processing for multiple molecules. Continues on individual failures.

#### `is_ready() -> bool`
Check if service is initialized and ready for predictions.

## Service Patterns

### Initialization Pattern

Services use lazy initialization:
```python
def predict_spectrum(self, smiles: str):
    if not self._initialized:
        if not self.initialize():
            raise ModelError("Service initialization failed")
    # ... continue with prediction
```

### Error Handling Pattern

Services provide user-friendly error messages:
```python
try:
    result = self._process(data)
except SpecificError as e:
    logger.error(f"Specific error: {e}")
    raise ServiceError(f"User-friendly message: {e}")
except Exception as e:
    logger.error(f"Unexpected error: {e}")
    raise ServiceError("An unexpected error occurred")
```

### Resource Management

Services manage resources efficiently:
- Models loaded once and reused
- Processors initialized on startup
- Memory-efficient batch processing

## Extending Services

### Adding New Services

1. Create new service class in `/services/`
2. Follow the initialization pattern
3. Implement error handling
4. Add logging for debugging

Example structure:
```python
class NewService:
    def __init__(self):
        self._component = None
        self._initialized = False
    
    def initialize(self) -> bool:
        try:
            self._component = Component()
            self._initialized = True
            return True
        except Exception as e:
            logger.error(f"Initialization failed: {e}")
            return False
    
    def process(self, data):
        if not self._initialized:
            self.initialize()
        # Implementation
```

### Service Dependencies

Services can depend on other services:
```python
class ComplexService:
    def __init__(self):
        self.prediction_service = PredictionService()
        self.other_service = OtherService()
```

## Configuration

Services use settings from `config/settings.py`:
- Model paths
- Processing parameters
- Feature flags

## Testing Services

Services are designed for easy testing:
```python
def test_prediction_service():
    service = PredictionService()
    assert service.initialize()
    
    result = service.predict_spectrum("CCO")
    assert result.smiles == "CCO"
    assert len(result.spectrum.x) > 0
```

## Performance

- Services cache expensive operations
- Batch methods optimize throughput
- Async support planned for future

## Best Practices

1. **Single Responsibility**: Each service handles one domain
2. **Fail Fast**: Validate inputs early
3. **Log Everything**: Detailed logging for debugging
4. **User-Friendly Errors**: Transform technical errors
5. **Resource Cleanup**: Properly manage resources