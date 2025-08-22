# Core Business Logic

The core module contains the fundamental business logic for mass spectrometry prediction, including machine learning models, molecular processing, and data structures.

## Structure

```
core/
├── ml/              # Machine learning components
├── processors/      # Data processing pipelines
└── models/          # Data models and structures
```

## Machine Learning (`/ml`)

### ModelHandler

The main interface for ML model operations:

```python
from core.ml.model_handler import ModelHandler

# Initialize
handler = ModelHandler()
handler.initialize(model_path, preprocessor_path, feature_mapping_path)

# Predict
features = MolecularFeatures(...)
spectrum = handler.predict(features)
```

Key features:
- Lazy loading of models
- Feature mapping support for consistency
- Preprocessor integration
- Error handling with fallbacks

## Processors (`/processors`)

### FeatureProcessor

Extracts comprehensive molecular features from SMILES:

```python
from core.processors.feature_processor import FeatureProcessor

processor = FeatureProcessor()
features = processor.extract_features("CCO")  # Returns MolecularFeatures
info = processor.extract_molecule_info("CCO")  # Returns MoleculeInfo
```

Features extracted:
- **RDKit Descriptors**: 200+ molecular properties
- **Fingerprints**: Morgan, MACCS, RDKit, Avalon, Pattern, Layered
- **Electronic Features**: Partial charges, PEOE VSA descriptors
- **Structural Counts**: Atoms, bonds by type

### SpectrumProcessor

Handles spectrum data operations:

```python
from core.processors.spectrum_processor import SpectralProcessor

processor = SpectralProcessor()

# Create from prediction
spectrum = processor.create_spectrum_from_prediction(intensities)

# Create from peaks
spectrum = processor.create_spectrum_from_peaks(peak_list)

# Process peaks for ML
processed = processor.process_peaks(peak_data)
```

Operations:
- Binning and normalization
- Peak filtering and sorting
- Format conversion
- Padding for consistent dimensions

### FeaturePreprocessor

Prepares features for ML models:

```python
from core.processors.feature_preprocessor import FeaturePreprocessor

preprocessor = FeaturePreprocessor()
preprocessor.fit(X_train)
X_scaled = preprocessor.transform(X_new)

# Save/load for production
preprocessor.save("preprocessor.pkl")
loaded = FeaturePreprocessor.load("preprocessor.pkl")
```

Preprocessing steps:
- NaN/infinity handling
- Binary feature detection
- Variance-based filtering
- Sign-preserving log transform
- Standard scaling

## Models (`/models`)

### Data Structures

#### MolecularFeatures
```python
@dataclass
class MolecularFeatures:
    smiles: str
    descriptors: Optional[np.ndarray]
    descriptor_names: Optional[List[str]]
    fingerprints: Dict[str, np.ndarray]
    electronic_features: Optional[np.ndarray]
    atom_counts: Dict[str, int]
    bond_counts: Dict[str, int]
```

#### Spectrum
```python
@dataclass
class Spectrum:
    x: List[float]  # m/z values
    y: List[float]  # intensities
    peaks: List[Peak]  # significant peaks
```

#### PredictionResult
```python
@dataclass
class PredictionResult:
    smiles: str
    chemical_name: str
    molecular_weight: float
    exact_mass: float
    spectrum: Spectrum
    structure_png: Optional[str]
```

## Feature Extraction Pipeline

1. **SMILES Validation**: Ensure valid molecular structure
2. **Descriptor Calculation**: Compute RDKit descriptors
3. **Fingerprint Generation**: Multiple fingerprint types
4. **Electronic Features**: Charge distribution properties
5. **Feature Assembly**: Combine into MolecularFeatures

## ML Pipeline Flow

```
SMILES → FeatureProcessor → MolecularFeatures
                                ↓
                        FeaturePreprocessor
                                ↓
                          ModelHandler
                                ↓
                        Raw Predictions
                                ↓
                        SpectrumProcessor
                                ↓
                          Final Spectrum
```

## Configuration

Feature extraction is configured in `config/settings.py`:

```python
features = FeatureConfig(
    extract_descriptors=True,
    fingerprints={
        'morgan': {'enabled': True, 'radii': [1, 2, 3], 'size': 1024},
        'maccs': {'enabled': True},
        # ... more fingerprints
    },
    extract_electronic=True
)
```

## Performance Considerations

- **Feature Caching**: Features are computed once per molecule
- **Batch Processing**: Processors support batch operations
- **Memory Management**: Large arrays are handled efficiently
- **Parallel Processing**: Can leverage multiple cores

## Error Handling

All processors include comprehensive error handling:
- Invalid SMILES → SMILESError
- Model failures → ModelError
- Feature extraction issues → Graceful degradation

## Extending the Core

### Adding New Features

1. Update `FeatureProcessor._extract_*` methods
2. Add to `MolecularFeatures` dataclass
3. Update feature mapping if needed

### Adding New Processors

1. Create new processor class
2. Follow existing patterns (init, process, error handling)
3. Add configuration to settings

### Model Updates

1. Retrain model with new features
2. Update preprocessor
3. Generate new feature mapping
4. Place files in `backend/models/`