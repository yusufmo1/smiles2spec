# Backend Utilities

The utils module provides shared functionality across the SMILES2SPEC backend, including chemistry operations, data processing, error handling, format conversion, and logging utilities.

## Module Overview

```
utils/
├── __init__.py         # Module initialization and exports
├── chemistry.py        # SMILES validation and molecular operations
├── data.py            # Data transformation and validation utilities  
├── errors.py          # Custom exception classes and error handling
├── formats.py         # File format conversion (MSP, JCAMP, etc.)
└── logging.py         # Centralized logging configuration
```

## Chemistry Utilities (`chemistry.py`)

### SMILES Validation and Processing

Core molecular chemistry operations using RDKit:

```python
from utils.chemistry import (
    validate_smiles,
    normalize_smiles, 
    get_molecule_info,
    calculate_molecular_properties,
    generate_structure_image
)

# Validate SMILES string
is_valid = validate_smiles("CCO")  # Returns: True
is_valid = validate_smiles("invalid")  # Returns: False

# Normalize SMILES (canonicalization)
canonical = normalize_smiles("OCC")  # Returns: "CCO"

# Get basic molecule information
info = get_molecule_info("CCO")
# Returns: MoleculeInfo(formula="C2H6O", mw=46.07, ...)
```

### Molecular Property Calculation

```python
# Calculate comprehensive molecular properties
properties = calculate_molecular_properties("CCO")
# Returns: {
#   'molecular_weight': 46.07,
#   'exact_mass': 46.0419,
#   'formula': 'C2H6O',
#   'logp': -0.31,
#   'tpsa': 20.23,
#   'rotatable_bonds': 1,
#   'h_donors': 1,
#   'h_acceptors': 1,
#   'ring_count': 0,
#   'aromatic_rings': 0
# }
```

### Structure Image Generation

```python
# Generate molecular structure as base64 PNG
structure_png = generate_structure_image(
    smiles="CCO",
    size=(300, 300),
    format="PNG"
)
# Returns: "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAA..."

# Generate SVG for scalable graphics
structure_svg = generate_structure_image(
    smiles="CCO", 
    format="SVG",
    highlight_atoms=None,
    show_atom_labels=True
)
```

### Advanced Chemistry Functions

```python
# SMILES standardization pipeline
def standardize_smiles(smiles: str) -> str:
    """
    Standardize SMILES through multiple operations:
    - Remove salts and solvents
    - Normalize charges
    - Canonicalize tautomers
    - Generate canonical SMILES
    """
    
# Molecular similarity calculations
def calculate_similarity(smiles1: str, smiles2: str, method="tanimoto") -> float:
    """Calculate Tanimoto similarity between molecules"""
    
# Substructure matching
def has_substructure(smiles: str, pattern: str) -> bool:
    """Check if molecule contains specific substructure"""
    
# Molecular descriptors batch calculation
def calculate_descriptors_batch(smiles_list: List[str]) -> pd.DataFrame:
    """Calculate RDKit descriptors for multiple molecules"""
```

## Data Utilities (`data.py`)

### Data Transformation

Utilities for data format conversion and validation:

```python
from utils.data import (
    validate_spectrum_data,
    normalize_spectrum,
    convert_peaks_to_spectrum,
    interpolate_spectrum_data,
    validate_bulk_smiles_data
)

# Validate spectrum data structure
is_valid = validate_spectrum_data({
    'x': [0, 1, 2, 3],
    'y': [0.1, 0.5, 0.3, 0.0]
})

# Normalize spectrum to 0-1 range
normalized = normalize_spectrum(spectrum_data, method="max")

# Convert peak list to full spectrum
spectrum = convert_peaks_to_spectrum(
    peaks=[{'mz': 45, 'intensity': 100}, {'mz': 46, 'intensity': 30}],
    mz_range=(0, 500),
    resolution=1
)
```

### Bulk Data Processing

```python
# Process bulk SMILES uploads
def process_bulk_smiles(
    file_content: str, 
    format: str = "auto"
) -> List[Dict]:
    """
    Process bulk SMILES from various formats:
    - CSV files with SMILES column
    - Plain text (one SMILES per line)
    - JSON arrays
    - SDF files (future)
    """

# Validate CSV structure
def validate_csv_structure(content: str) -> Dict:
    """
    Analyze CSV structure and detect:
    - Column headers
    - SMILES column candidates
    - Data quality issues
    - Recommended processing approach
    """

# Data cleaning utilities
def clean_smiles_list(smiles_list: List[str]) -> List[str]:
    """Remove invalid, duplicate, and problematic SMILES"""
```

### Data Validation

```python
# Comprehensive data validators
class DataValidator:
    @staticmethod
    def validate_prediction_request(data: dict) -> ValidationResult:
        """Validate prediction API request data"""
    
    @staticmethod
    def validate_export_request(data: dict) -> ValidationResult:
        """Validate export API request data"""
    
    @staticmethod
    def validate_file_upload(file_data: bytes, filename: str) -> ValidationResult:
        """Validate uploaded file content and format"""

# Spectrum data validation
def validate_spectrum_format(spectrum: dict) -> bool:
    """Ensure spectrum has required fields and valid data"""

# Peak list validation  
def validate_peak_list(peaks: List[dict]) -> bool:
    """Validate peak list structure and ranges"""
```

## Error Handling (`errors.py`)

### Custom Exception Hierarchy

Structured exception classes for different error types:

```python
from utils.errors import (
    SMILES2SpecError,
    SMILESError, 
    ModelError,
    ValidationError,
    ProcessingError,
    FileFormatError
)

# Base exception class
class SMILES2SpecError(Exception):
    """Base exception for all SMILES2SPEC errors"""
    
    def __init__(self, message: str, error_code: str = None, details: dict = None):
        self.message = message
        self.error_code = error_code
        self.details = details or {}
        super().__init__(message)

# Specific error types
class SMILESError(SMILES2SpecError):
    """Errors related to SMILES processing"""
    pass

class ModelError(SMILES2SpecError):
    """Errors related to ML model operations"""
    pass

class ValidationError(SMILES2SpecError):
    """Data validation errors"""
    pass

class ProcessingError(SMILES2SpecError):
    """General processing errors"""
    pass

class FileFormatError(SMILES2SpecError):
    """File format and parsing errors"""
    pass
```

### Error Handling Utilities

```python
# Error context managers
from contextlib import contextmanager

@contextmanager
def handle_chemistry_errors():
    """Context manager for RDKit operations"""
    try:
        yield
    except Exception as e:
        if "RDKit" in str(e):
            raise SMILESError(f"Invalid molecular structure: {e}")
        raise ProcessingError(f"Chemistry operation failed: {e}")

# Usage example
with handle_chemistry_errors():
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SMILESError("Invalid SMILES string")

# Error response formatting
def format_error_response(error: SMILES2SpecError) -> dict:
    """Format error for API response"""
    return {
        'error': error.message,
        'error_code': error.error_code,
        'details': error.details,
        'timestamp': datetime.utcnow().isoformat()
    }
```

### Error Recovery

```python
# Retry decorators for transient errors
from functools import wraps
import time

def retry_on_failure(max_retries=3, delay=1.0, backoff=2.0):
    """Decorator for retrying failed operations"""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            retries = 0
            while retries < max_retries:
                try:
                    return func(*args, **kwargs)
                except (ProcessingError, ModelError) as e:
                    retries += 1
                    if retries >= max_retries:
                        raise e
                    time.sleep(delay * (backoff ** (retries - 1)))
            return None
        return wrapper
    return decorator

# Usage
@retry_on_failure(max_retries=3)
def predict_spectrum_with_retry(smiles: str):
    """Predict spectrum with automatic retry on failure"""
    return model_handler.predict(smiles)
```

## Format Conversion (`formats.py`)

### MSP Format Support

Mass Spectral Library format conversion:

```python
from utils.formats import (
    export_to_msp,
    parse_msp_file,
    export_to_csv,
    export_to_json,
    export_batch_msp
)

# Export single spectrum to MSP
msp_content = export_to_msp(
    smiles="CCO",
    spectrum_data=spectrum,
    metadata={
        'name': 'Ethanol',
        'molecular_weight': 46.07,
        'source': 'SMILES2SPEC_v1.0'
    }
)

# Batch export multiple spectra
batch_msp = export_batch_msp([
    {
        'smiles': 'CCO',
        'spectrum': spectrum1,
        'name': 'Ethanol'
    },
    {
        'smiles': 'CCC', 
        'spectrum': spectrum2,
        'name': 'Propane'
    }
])
```

### CSV Export

```python
# Export spectrum data as CSV
csv_content = export_to_csv(
    spectrum_data=spectrum,
    include_metadata=True,
    decimal_places=4
)

# Export peak table
peak_csv = export_peaks_to_csv(
    peaks=peak_list,
    columns=['mz', 'intensity', 'relative_intensity', 'annotation']
)
```

### JSON Export

```python
# Export complete prediction data as JSON
json_content = export_to_json(
    prediction_result=result,
    include_structure=True,
    pretty_print=True
)

# Export with compression
compressed_json = export_to_json(
    data=large_dataset,
    compress=True,
    compression_level=9
)
```

### JCAMP-DX Support (Future)

```python
# JCAMP-DX format for spectral data exchange
def export_to_jcamp(spectrum_data: dict, metadata: dict) -> str:
    """Export spectrum in JCAMP-DX format"""
    
def parse_jcamp_file(content: str) -> dict:
    """Parse JCAMP-DX file content"""
```

## Logging Configuration (`logging.py`)

### Centralized Logging Setup

```python
from utils.logging import setup_logging, get_logger

# Initialize logging system
setup_logging(
    log_level="INFO",
    log_file="smiles2spec.log",
    enable_console=True,
    enable_file=True,
    max_bytes=10*1024*1024,  # 10MB
    backup_count=5
)

# Get module-specific logger
logger = get_logger(__name__)

# Usage in modules
logger.info("Starting spectrum prediction")
logger.error("Failed to process SMILES: %s", smiles)
logger.debug("Feature extraction completed in %.2fs", elapsed_time)
```

### Logging Configuration

```python
# Structured logging configuration
LOGGING_CONFIG = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'detailed': {
            'format': '%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s'
        },
        'simple': {
            'format': '%(levelname)s - %(message)s'
        },
        'json': {
            'format': '{"timestamp": "%(asctime)s", "level": "%(levelname)s", "module": "%(name)s", "message": "%(message)s"}'
        }
    },
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'level': 'INFO',
            'formatter': 'simple',
            'stream': 'ext://sys.stdout'
        },
        'file': {
            'class': 'logging.handlers.RotatingFileHandler',
            'level': 'DEBUG',
            'formatter': 'detailed',
            'filename': 'logs/smiles2spec.log',
            'maxBytes': 10485760,  # 10MB
            'backupCount': 5
        },
        'error_file': {
            'class': 'logging.handlers.RotatingFileHandler',
            'level': 'ERROR',
            'formatter': 'json',
            'filename': 'logs/errors.log',
            'maxBytes': 5242880,  # 5MB
            'backupCount': 3
        }
    },
    'loggers': {
        'smiles2spec': {
            'level': 'DEBUG',
            'handlers': ['console', 'file', 'error_file'],
            'propagate': False
        },
        'rdkit': {
            'level': 'WARNING',
            'handlers': ['file'],
            'propagate': False
        }
    },
    'root': {
        'level': 'INFO',
        'handlers': ['console']
    }
}
```

### Performance Logging

```python
# Performance monitoring utilities
import time
from functools import wraps

def log_performance(operation_name: str = None):
    """Decorator to log operation performance"""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            operation = operation_name or func.__name__
            start_time = time.time()
            
            logger = get_logger(func.__module__)
            logger.debug(f"Starting {operation}")
            
            try:
                result = func(*args, **kwargs)
                elapsed = time.time() - start_time
                logger.info(f"{operation} completed in {elapsed:.3f}s")
                return result
            except Exception as e:
                elapsed = time.time() - start_time
                logger.error(f"{operation} failed after {elapsed:.3f}s: {e}")
                raise
        return wrapper
    return decorator

# Usage
@log_performance("spectrum_prediction")
def predict_spectrum(smiles: str):
    """Predict mass spectrum for given SMILES"""
    # Implementation
```

### Request Logging

```python
# HTTP request logging for API endpoints
def log_request(request, response, elapsed_time):
    """Log API request details"""
    logger = get_logger('api')
    
    log_data = {
        'method': request.method,
        'path': request.path,
        'status_code': response.status_code,
        'elapsed_time': elapsed_time,
        'user_agent': request.headers.get('User-Agent'),
        'ip_address': request.remote_addr,
        'content_length': response.content_length
    }
    
    if response.status_code >= 400:
        logger.warning("Request failed: %s", log_data)
    else:
        logger.info("Request completed: %s", log_data)
```

## Usage Examples

### Complete Error Handling Pipeline

```python
from utils.errors import SMILESError, ValidationError, format_error_response
from utils.logging import get_logger
from utils.chemistry import validate_smiles, get_molecule_info

logger = get_logger(__name__)

def process_smiles_safely(smiles: str) -> dict:
    """Process SMILES with comprehensive error handling"""
    try:
        # Validate input
        if not smiles or not isinstance(smiles, str):
            raise ValidationError("SMILES must be a non-empty string")
        
        # Chemistry validation
        if not validate_smiles(smiles):
            raise SMILESError(f"Invalid SMILES string: {smiles}")
        
        # Process molecule
        logger.info(f"Processing SMILES: {smiles}")
        molecule_info = get_molecule_info(smiles)
        
        logger.info(f"Successfully processed {molecule_info.formula}")
        return {
            'success': True,
            'data': molecule_info.__dict__
        }
        
    except (SMILESError, ValidationError) as e:
        logger.warning(f"Processing failed for {smiles}: {e}")
        return format_error_response(e)
    
    except Exception as e:
        logger.error(f"Unexpected error processing {smiles}: {e}")
        error = ProcessingError("An unexpected error occurred")
        return format_error_response(error)
```

### Bulk Data Processing

```python
from utils.data import process_bulk_smiles, validate_csv_structure
from utils.formats import export_batch_msp
from utils.logging import get_logger

logger = get_logger(__name__)

def process_bulk_upload(file_content: str, filename: str) -> dict:
    """Process bulk SMILES upload with validation and export"""
    try:
        # Validate file structure
        if filename.endswith('.csv'):
            validation = validate_csv_structure(file_content)
            if not validation['is_valid']:
                raise FileFormatError(validation['error_message'])
        
        # Process SMILES
        logger.info(f"Processing bulk upload: {filename}")
        smiles_data = process_bulk_smiles(file_content)
        
        # Generate predictions for all valid SMILES
        results = []
        for item in smiles_data:
            if item['is_valid']:
                # Process prediction
                result = predict_spectrum(item['smiles'])
                results.append(result)
        
        # Export results
        msp_content = export_batch_msp(results)
        
        logger.info(f"Processed {len(results)} molecules from {filename}")
        return {
            'success': True,
            'processed_count': len(results),
            'total_count': len(smiles_data),
            'export_data': msp_content
        }
        
    except Exception as e:
        logger.error(f"Bulk processing failed: {e}")
        return format_error_response(ProcessingError(str(e)))
```

## Testing Utilities

### Test Helpers

```python
# Test data generators
def generate_test_smiles(count: int = 10) -> List[str]:
    """Generate valid test SMILES strings"""
    return [
        "CCO",           # Ethanol
        "CCC",           # Propane  
        "c1ccccc1",      # Benzene
        "CC(=O)O",       # Acetic acid
        "CC(C)O",        # Isopropanol
        # ... more test molecules
    ][:count]

def generate_test_spectrum(length: int = 500) -> dict:
    """Generate synthetic test spectrum data"""
    import numpy as np
    
    x = list(range(length))
    y = np.random.exponential(0.1, length).tolist()
    
    return {'x': x, 'y': y}

# Mock data for testing
class MockMoleculeInfo:
    def __init__(self, smiles: str):
        self.smiles = smiles
        self.formula = "C2H6O"  # Default formula
        self.molecular_weight = 46.07
        self.exact_mass = 46.0419
```

### Validation Testing

```python
import pytest
from utils.chemistry import validate_smiles
from utils.data import validate_spectrum_data
from utils.errors import SMILESError

def test_smiles_validation():
    """Test SMILES validation functionality"""
    # Valid SMILES
    assert validate_smiles("CCO") == True
    assert validate_smiles("c1ccccc1") == True
    
    # Invalid SMILES
    assert validate_smiles("invalid") == False
    assert validate_smiles("") == False
    assert validate_smiles(None) == False

def test_spectrum_validation():
    """Test spectrum data validation"""
    valid_spectrum = {'x': [0, 1, 2], 'y': [0.1, 0.5, 0.2]}
    assert validate_spectrum_data(valid_spectrum) == True
    
    invalid_spectrum = {'x': [0, 1], 'y': [0.1]}  # Mismatched lengths
    assert validate_spectrum_data(invalid_spectrum) == False

def test_error_handling():
    """Test custom error classes"""
    with pytest.raises(SMILESError):
        raise SMILESError("Test SMILES error")
```

## Performance Considerations

### Optimization Strategies

```python
# Caching for expensive operations
from functools import lru_cache

@lru_cache(maxsize=1000)
def get_molecule_properties_cached(smiles: str) -> dict:
    """Cache molecular property calculations"""
    return calculate_molecular_properties(smiles)

# Batch processing optimization
def process_smiles_batch(smiles_list: List[str], batch_size: int = 100):
    """Process SMILES in batches for memory efficiency"""
    for i in range(0, len(smiles_list), batch_size):
        batch = smiles_list[i:i + batch_size]
        yield process_batch(batch)

# Memory-efficient file processing
def process_large_file_stream(file_path: str):
    """Process large files without loading entire content"""
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip():
                yield process_smiles_line(line.strip())
```

## Best Practices

### Error Handling
1. **Use specific exception types** for different error categories
2. **Include context** in error messages for debugging
3. **Log errors appropriately** with proper log levels
4. **Provide user-friendly messages** while preserving technical details

### Data Processing
1. **Validate inputs early** to catch issues before processing
2. **Use consistent data formats** across the application
3. **Handle edge cases** like empty data or malformed input
4. **Optimize for common use cases** while supporting edge cases

### Performance
1. **Cache expensive computations** using appropriate cache strategies
2. **Process data in batches** for memory efficiency
3. **Use appropriate data structures** for different operations
4. **Profile and optimize** bottlenecks in critical paths

### Logging
1. **Use appropriate log levels** (DEBUG, INFO, WARNING, ERROR)
2. **Include relevant context** in log messages
3. **Avoid logging sensitive data** like API keys
4. **Use structured logging** for better searchability

## Future Enhancements

- **Async Support**: Asynchronous versions of utility functions
- **Database Integration**: Utilities for database operations
- **Caching Layer**: Redis-based caching for distributed systems
- **Monitoring**: Integration with monitoring and alerting systems
- **Configuration Management**: Dynamic configuration updates
- **Plugin System**: Extensible utility framework