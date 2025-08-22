# API Layer

The API layer handles all HTTP communication between the frontend and backend services. It provides RESTful endpoints for spectrum prediction, chat functionality, and data export.

## Structure

```
api/
├── routes/          # Endpoint handlers
├── middleware/      # Request/response processing
└── schemas/         # Data validation models
```

## Routes

### `/routes/prediction.py`

Handles spectrum prediction endpoints:
- `POST /predict` - Main prediction endpoint
- `GET /structure` - Molecular structure visualization

### `/routes/chat.py`

AI-powered chat functionality:
- `POST /chat` - Chat with AI about spectra (supports streaming)
- `POST /generate_smiles` - Generate SMILES from text descriptions

### `/routes/export.py`

Data export endpoints:
- `POST /export/msp` - Export single spectrum as MSP file
- `POST /export/csv` - Export spectrum data as CSV
- `POST /export/json` - Export spectrum data as JSON

### `/routes/upload.py`

File upload handling:
- `POST /upload` - Process CSV/TXT files with SMILES

### `/routes/health.py`

System monitoring:
- `GET /health` - Service health check and status

## Middleware

### `/middleware/error_handler.py`

Global error handling that ensures consistent error responses:
- Catches all exceptions and formats them properly
- Logs errors for debugging
- Returns appropriate HTTP status codes

### `/middleware/validation.py`

Request validation decorators:
- `@validate_json(schema)` - Validates JSON body against Pydantic schema
- `@validate_query_params(schema)` - Validates query parameters
- Automatic error messages for invalid requests

## Schemas

### `/schemas/prediction.py`

Data models for prediction endpoints:
- `PredictRequest` - SMILES input validation
- `ExportMSPRequest` - Single MSP export
- `ExportMSPBatchRequest` - Batch MSP export

### `/schemas/chat.py`

Chat-related models:
- `ChatMessage` - Individual message structure
- `ChatRequest` - Full chat request with context
- `GenerateSMILESRequest` - SMILES generation parameters

## Usage Examples

### Basic Prediction Request

```python
# Request
POST /predict
{
    "smiles": "CCO"
}

# Response
{
    "smiles": "CCO",
    "chemical_name": "Ethanol",
    "molecular_weight": 46.07,
    "spectrum": {
        "x": [0, 1, 2, ..., 499],
        "y": [0.0, 0.1, 0.5, ..., 0.0]
    },
    "peaks": [
        {"mz": 45, "intensity": 0.999},
        {"mz": 46, "intensity": 0.30}
    ]
}
```

### Streaming Chat Request

```python
# Request
POST /chat
{
    "messages": [
        {"role": "user", "content": "Explain this spectrum"}
    ],
    "smiles": "CCO",
    "stream": true
}

# Response (Server-Sent Events)
data: {"chunk": "The spectrum of ethanol shows..."}
data: {"chunk": " characteristic peaks at..."}
data: [DONE]
```

## Error Handling

All endpoints return consistent error responses:

```json
{
    "error": "Error message",
    "details": "Additional context (optional)"
}
```

Common HTTP status codes:
- `400` - Bad Request (invalid input)
- `404` - Not Found
- `500` - Internal Server Error

## Adding New Endpoints

1. Create route handler in appropriate file or new file in `/routes/`
2. Define Pydantic schema in `/schemas/` if needed
3. Use validation decorators for input validation
4. Register blueprint in `main.py`

Example:

```python
# routes/new_feature.py
from flask import Blueprint, jsonify
from ..schemas.new_feature import NewFeatureRequest
from ..middleware.validation import validate_json

new_feature_bp = Blueprint('new_feature', __name__)

@new_feature_bp.route('/new_feature', methods=['POST'])
@validate_json(NewFeatureRequest)
def new_feature(validated_data: NewFeatureRequest):
    """New feature endpoint."""
    # Implementation here
    return jsonify({"result": "success"})
```

## Best Practices

1. **Always validate input** using Pydantic schemas
2. **Use blueprints** for route organization
3. **Return consistent responses** across all endpoints
4. **Log errors** but don't expose internal details to clients
5. **Document endpoints** with docstrings
6. **Handle streaming properly** for real-time features