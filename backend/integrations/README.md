# External Integrations

This module contains integrations with external services, currently focused on Large Language Model (LLM) capabilities for enhanced user interactions.

## Structure

```
integrations/
└── llm/                 # LLM integration
    ├── client.py        # OpenRouter API client
    └── services/        # LLM-powered services
        ├── chat_service.py
        └── smiles_service.py
```

## LLM Integration

### Overview

The LLM integration uses OpenRouter as a gateway to access various language models including:
- google/gemma-3-27b-it:free (primary)
- GPT-4
- Other models as configured

### LLMClient

Base client for API communication:

```python
from integrations.llm.client import LLMClient

client = LLMClient()
response = client.chat_completion(
    messages=[
        {"role": "system", "content": "You are a chemistry expert"},
        {"role": "user", "content": "What is benzene?"}
    ],
    model="google/gemma-3-27b-it:free",
    stream=False
)
```

Features:
- Automatic header configuration
- Streaming support (SSE)
- Error handling with retries
- Response parsing

### ChatService

AI-powered chat about mass spectrometry:

```python
from integrations.llm.services.chat_service import ChatService

service = ChatService()

# Basic chat
response = service.generate_response(
    messages=[{"role": "user", "content": "Explain this spectrum"}],
    spectrum_context=prediction_result,  # Optional context
    stream=False
)

# Streaming chat
for chunk in service.generate_response(messages, stream=True):
    print(chunk, end='')
```

Context features:
- Automatic spectrum context injection
- Molecular property awareness
- Scientific accuracy focus

### SMILESService

Generate SMILES from natural language:

```python
from integrations.llm.services.smiles_service import SMILESService

service = SMILESService()

# Generate SMILES
smiles_list = service.generate_smiles(
    description="aromatic ring with hydroxyl group",
    count=3
)
# Returns: ['c1ccc(O)cc1', 'Oc1ccccc1', 'c1cc(O)ccc1']
```

Features:
- Multiple candidate generation
- Validation of generated SMILES
- Fallback mechanisms
- Retry logic

## Configuration

Set up in `.env`:
```bash
OPENROUTER_API_KEY=your_key_here
```

Configuration in `config/settings.py`:
```python
llm = LLMConfig(
    openrouter_api_key=os.getenv('OPENROUTER_API_KEY'),
    site_url="https://smiles2spec.app",
    site_name="SMILES2Spec App"
)
```

## API Usage

### Chat Endpoint Integration

The chat routes use these services:

```python
@chat_bp.route('/chat', methods=['POST'])
def chat(validated_data):
    # Get spectrum context if SMILES provided
    spectrum_context = None
    if validated_data.smiles:
        spectrum_context = prediction_service.predict_spectrum(validated_data.smiles)
    
    # Generate response
    response = chat_service.generate_response(
        messages=validated_data.messages,
        spectrum_context=spectrum_context,
        stream=validated_data.stream
    )
```

### SMILES Generation Endpoint

```python
@chat_bp.route('/generate_smiles', methods=['POST'])
def generate_smiles_route(validated_data):
    smiles_list = smiles_service.generate_smiles(
        description=validated_data.description,
        count=validated_data.count
    )
    return jsonify({"smiles": smiles_list})
```

## Streaming Implementation

For real-time chat responses:

1. **Client sends** streaming request
2. **Server establishes** SSE connection
3. **LLM streams** tokens
4. **Server forwards** as chunks
5. **Client renders** progressively

Example SSE format:
```
data: {"chunk": "The spectrum shows"}
data: {"chunk": " a molecular ion peak"}
data: [DONE]
```

## Error Handling

### Network Errors
- Automatic retry with exponential backoff
- Graceful fallback responses
- User-friendly error messages

### Invalid Responses
- Response validation
- Parsing error recovery
- Logging for debugging

### Rate Limiting
- Respect API limits
- Queue management (future)
- Cost tracking

## Extending Integrations

### Adding New LLM Providers

1. Create new client class:
```python
class NewLLMClient:
    def __init__(self):
        self.api_key = settings.new_llm.api_key
    
    def chat_completion(self, messages, **kwargs):
        # Implementation
```

2. Update service to use new client:
```python
class ChatService:
    def __init__(self, client_type='openrouter'):
        if client_type == 'new_provider':
            self.client = NewLLMClient()
```

### Adding New Services

1. Create service class in `/services/`
2. Implement core methods
3. Add error handling
4. Create route handler

## Best Practices

1. **API Key Security**: Never commit API keys
2. **Error Messages**: Don't expose internal errors
3. **Streaming**: Handle connection drops gracefully
4. **Context Limits**: Monitor token usage
5. **Response Quality**: Validate generated content

## Performance Considerations

- Stream for long responses
- Cache common queries (future)
- Batch similar requests
- Monitor API costs

## Future Enhancements

- Additional LLM providers
- Embedding-based search
- Fine-tuned models for chemistry
- Multi-modal inputs (structure images)
- Conversation memory