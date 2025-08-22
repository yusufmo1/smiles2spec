# SMILES Input System

The SMILES input system provides a comprehensive interface for molecular structure input, supporting text entry, file uploads, AI generation, and batch processing.

## Architecture Overview

```
smiles-input/
├── SmilesInput.svelte        # Main input component
├── index.js                  # Exports
├── components/               # Sub-components
│   ├── ActionButtons.svelte  # Action controls
│   ├── DragDropZone.svelte   # File drop interface
│   ├── GenerateModal.svelte  # AI generation modal
│   ├── SmilesTextArea.svelte # Text input area
│   └── UploadModal.svelte    # File upload modal
├── services/                 # Business logic
│   └── smilesInputService.js # Input processing service
└── stores/                   # State management
    └── smilesInputStore.ts   # Input state store
```

## Core Component

### SmilesInput.svelte

The main input interface that orchestrates all input methods:

```svelte
<script>
  import { SmilesInput } from '$lib/components/smiles-input';
</script>

<SmilesInput
  on:predict={handlePredict}
  on:bulkProcess={handleBulkProcess}
  placeholder="Enter SMILES string or upload file..."
/>
```

**Features:**

- Multi-line text input
- Real-time validation
- Bulk SMILES detection
- File upload integration
- AI-powered generation
- Format detection

**Events:**

- `predict` - Single SMILES prediction request
- `bulkProcess` - Bulk SMILES processing request
- `clear` - Input cleared
- `error` - Validation or processing error

## Sub-Components

### SmilesTextArea.svelte

Core text input with enhanced functionality:

```svelte
<SmilesTextArea
  value={smilesText}
  placeholder="Enter SMILES..."
  on:input={handleInput}
  on:validate={handleValidation}
/>
```

**Features:**

- Auto-resizing textarea
- Real-time SMILES validation
- Syntax highlighting (future)
- Multi-line support
- Keyboard shortcuts

### ActionButtons.svelte

Action controls for input operations:

```svelte
<ActionButtons {canPredict} {canClear} {isBulk} on:predict on:clear on:generate on:upload />
```

**Button States:**

- **Predict** - Enabled when valid SMILES present
- **Clear** - Enabled when input has content
- **Generate** - Always available for AI generation
- **Upload** - Always available for file upload

### DragDropZone.svelte

File drop interface for bulk uploads:

```svelte
<DragDropZone
  accept=".csv,.txt"
  on:filesDropped={handleFiles}
  on:dragEnter={handleDragEnter}
  on:dragLeave={handleDragLeave}
/>
```

**Features:**

- Drag and drop file upload
- Visual drag feedback
- File type validation
- Error handling
- Progressive enhancement

**Supported Formats:**

- CSV files with SMILES column
- Plain text files (one SMILES per line)
- Maximum file size: 10MB

### UploadModal.svelte

File upload modal with advanced options:

```svelte
<UploadModal isOpen={showUpload} on:close={closeModal} on:upload={processFile} />
```

**Features:**

- File format selection
- Column mapping for CSV
- Preview before processing
- Upload progress
- Error reporting

### GenerateModal.svelte

AI-powered SMILES generation interface:

```svelte
<GenerateModal isOpen={showGenerate} on:close={closeModal} on:generate={processGeneration} />
```

**Features:**

- Natural language to SMILES conversion
- Multiple generation options
- Chemical category selection
- Validation of generated SMILES
- Batch generation support

## Services

### smilesInputService.js

Core business logic for input processing:

```javascript
import { smilesInputService } from '$lib/components/smiles-input/services';

// Validate SMILES
const isValid = await smilesInputService.validateSmiles('CCO');

// Process bulk input
const results = await smilesInputService.processBulkInput(fileContent);

// Generate SMILES
const generated = await smilesInputService.generateSmiles('aromatic alcohol');
```

**Methods:**

#### `validateSmiles(smiles: string): Promise<boolean>`

Validates a single SMILES string using RDKit validation.

#### `processBulkInput(content: string): Promise<string[]>`

Processes bulk input from text or CSV, extracting valid SMILES.

#### `detectFormat(content: string): 'single' | 'multiple' | 'csv'`

Automatically detects input format and suggests processing method.

#### `generateSmiles(description: string, count?: number): Promise<string[]>`

Generates SMILES from natural language description using AI.

#### `parseCSV(content: string, smilesColumn?: string): string[]`

Parses CSV content and extracts SMILES from specified column.

## State Management

### smilesInputStore.ts

Manages input state and processing status:

```typescript
import { smilesInputStore } from '$lib/components/smiles-input/stores';

// Current input state
$: inputValue = $smilesInputStore.value;
$: isValid = $smilesInputStore.isValid;
$: isBulk = $smilesInputStore.isBulk;

// Actions
smilesInputStore.setValue('CCO');
smilesInputStore.setValidation(true);
smilesInputStore.setBulkMode(true);
smilesInputStore.clear();
```

**State Properties:**

- `value` - Current input text
- `isValid` - Validation status
- `isBulk` - Bulk processing mode
- `isProcessing` - Processing status
- `error` - Error message
- `suggestions` - Auto-complete suggestions

**Actions:**

- `setValue()` - Update input value
- `setValidation()` - Set validation status
- `setBulkMode()` - Toggle bulk processing
- `setError()` - Set error message
- `clear()` - Reset all state

## Input Processing Flow

### Single SMILES Processing

```
User Input → Validation → Structure Display → Prediction Request
```

1. **Input**: User enters SMILES string
2. **Validation**: Real-time RDKit validation
3. **Display**: Show molecular structure preview
4. **Prediction**: Send to spectrum prediction API

### Bulk Processing Flow

```
File Upload → Format Detection → Parsing → Validation → Batch Processing
```

1. **Upload**: User drops/selects file
2. **Detection**: Automatic format detection
3. **Parsing**: Extract SMILES from file
4. **Validation**: Validate each SMILES
5. **Processing**: Batch prediction workflow

### AI Generation Flow

```
Description → LLM Request → SMILES Generation → Validation → Selection
```

1. **Description**: User provides natural language
2. **Request**: Send to AI service
3. **Generation**: Multiple SMILES candidates
4. **Validation**: Validate generated structures
5. **Selection**: User selects preferred SMILES

## Validation System

### Real-time Validation

```javascript
// Debounced validation
const validateInput = debounce(async (smiles) => {
  const isValid = await smilesInputService.validateSmiles(smiles);
  smilesInputStore.setValidation(isValid);
}, 300);
```

### Validation Types

1. **Syntax Validation**: Basic SMILES syntax checking
2. **Chemical Validation**: RDKit molecular validation
3. **Structure Validation**: 3D conformer generation test
4. **Aromaticity Validation**: Aromatic structure validation

### Error Handling

```javascript
try {
  const result = await processSmiles(input);
} catch (error) {
  if (error.type === 'VALIDATION_ERROR') {
    showValidationError(error.message);
  } else if (error.type === 'NETWORK_ERROR') {
    showNetworkError();
  } else {
    showGenericError();
  }
}
```

## File Processing

### Supported Formats

#### CSV Files

```csv
Name,SMILES,MW
Ethanol,CCO,46.07
Benzene,c1ccccc1,78.11
```

#### Text Files

```
CCO
c1ccccc1
CC(=O)O
```

### CSV Column Detection

Automatic detection of SMILES columns:

- Common column names: `smiles`, `SMILES`, `canonical_smiles`
- Pattern matching: Strings that look like SMILES
- User selection: Manual column selection option

### File Size Limits

- Maximum file size: 10MB
- Maximum SMILES count: 1000 per file
- Timeout: 30 seconds per processing job

## AI Integration

### Natural Language Processing

Users can generate SMILES from descriptions:

```
"aromatic alcohol" → ["c1ccc(O)cc1", "Oc1ccccc1"]
"simple alkane" → ["CCC", "CCCC", "CC"]
"drug-like molecule" → [complex SMILES]
```

### Generation Parameters

- **Count**: Number of SMILES to generate (1-10)
- **Category**: Chemical category hint
- **Complexity**: Simple, moderate, complex
- **Properties**: MW range, LogP range, etc.

## Performance Optimization

### Debouncing

Input validation is debounced to prevent excessive API calls:

```javascript
const debouncedValidation = debounce(validateSmiles, 300);
```

### Caching

Validation results are cached:

```javascript
const validationCache = new Map();

async function validateWithCache(smiles) {
  if (validationCache.has(smiles)) {
    return validationCache.get(smiles);
  }

  const result = await validateSmiles(smiles);
  validationCache.set(smiles, result);
  return result;
}
```

### Lazy Loading

Heavy processing is lazy-loaded:

```javascript
const processLargeFile = async (file) => {
  const { processCSV } = await import('./bulkProcessor.js');
  return processCSV(file);
};
```

## Accessibility

### Keyboard Navigation

- **Tab**: Navigate between input areas
- **Enter**: Trigger prediction (single line)
- **Ctrl+Enter**: Trigger prediction (multi-line)
- **Escape**: Close modals

### Screen Reader Support

```svelte
<textarea
  aria-label="SMILES input"
  aria-describedby="smiles-help"
  aria-invalid={!isValid}
>
```

### Error Announcements

```svelte
<div role="alert" aria-live="polite" class:hidden={!error}>
  {error}
</div>
```

## Usage Examples

### Basic Usage

```svelte
<script>
  import { SmilesInput } from '$lib/components/smiles-input';

  function handlePredict(event) {
    const { smiles } = event.detail;
    // Process single SMILES
  }

  function handleBulkProcess(event) {
    const { smilesList } = event.detail;
    // Process multiple SMILES
  }
</script>

<SmilesInput on:predict={handlePredict} on:bulkProcess={handleBulkProcess} />
```

### Advanced Configuration

```svelte
<SmilesInput
  placeholder="Enter molecular structure..."
  maxFileSize={5 * 1024 * 1024}
  supportedFormats={['.csv', '.txt', '.sdf']}
  enableGeneration={true}
  validationDelay={500}
  on:predict={handlePredict}
  on:error={handleError}
/>
```

### Custom Validation

```svelte
<script>
  import { SmilesInput } from '$lib/components/smiles-input';

  async function customValidator(smiles) {
    // Custom validation logic
    const isValid = await myValidationService(smiles);
    return { isValid, message: 'Custom validation message' };
  }
</script>

<SmilesInput {customValidator} on:predict={handlePredict} />
```

## Testing

### Component Testing

```javascript
import { render, fireEvent } from '@testing-library/svelte';
import SmilesInput from './SmilesInput.svelte';

test('validates SMILES input', async () => {
  const { getByTestId } = render(SmilesInput);
  const input = getByTestId('smiles-textarea');

  await fireEvent.input(input, { target: { value: 'CCO' } });

  // Wait for validation
  await waitFor(() => {
    expect(getByTestId('validation-status')).toHaveClass('valid');
  });
});
```

### Service Testing

```javascript
import { smilesInputService } from './services/smilesInputService';

test('validates SMILES correctly', async () => {
  const valid = await smilesInputService.validateSmiles('CCO');
  expect(valid).toBe(true);

  const invalid = await smilesInputService.validateSmiles('invalid');
  expect(invalid).toBe(false);
});
```

## Best Practices

### Input Handling

1. **Debounce validation** to prevent excessive API calls
2. **Cache results** for repeated validations
3. **Provide immediate feedback** for user actions
4. **Handle errors gracefully** with clear messages

### File Processing

1. **Validate file types** before processing
2. **Show progress** for large files
3. **Limit file sizes** to prevent performance issues
4. **Provide preview** before bulk processing

### State Management

1. **Keep state minimal** and focused
2. **Use reactive statements** for derived data
3. **Clear state** appropriately
4. **Handle async operations** properly

## Future Enhancements

- **SMILES Auto-complete**: Intelligent suggestions
- **Structure Editor**: Visual molecular editor
- **Format Conversion**: Support for more chemical formats
- **Batch Validation**: Parallel validation for bulk input
- **History**: Recent SMILES history
- **Favorites**: Save frequently used structures
