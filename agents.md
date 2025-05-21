# `agents.md` - LLM Integration Refactoring Plan

## Overview

This document outlines the implementation plan for refactoring the LLM integration to create separate components for:
1. The SMILES Generator - generates valid molecular structures from text descriptions
2. Chat with Spectra Assistant - interactive chemistry assistant with spectrum visualization

The refactoring will adhere to KISS (Keep It Simple, Stupid) and SRP (Single Responsibility Principle) while adding streaming responses and visual spectrum context.

## Architecture Changes

### Current Issues
- No clear separation between AI components
- No streaming response support
- Spectrum visualizations not integrated with chat
- Using outdated API request patterns
- Code organization needs improvement

### New Directory Structure
```
llm_integration/
├── __init__.py             # Export main functions
├── common/                 # Shared functionality
│   ├── __init__.py
│   ├── config.py           # Configuration and models
│   ├── streaming.py        # Streaming implementation
│   └── api.py              # API request handling
├── spectra_chat/           # Chat assistant module
│   ├── __init__.py
│   ├── service.py          # Chat implementation
│   ├── prompts.py          # System/user prompts
│   └── formatter.py        # Spectrum visualization
└── smiles_generator/       # SMILES generation module
    ├── __init__.py
    ├── service.py          # Generation implementation
    ├── prompts.py          # System prompts
    └── validator.py        # SMILES validation
```

## Backend Implementation

### 1. Common Modules

#### `common/config.py`
```python
"""Shared configuration for LLM integration."""
import os
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# API configurations
OPENROUTER_API_KEY = os.environ.get("OPENROUTER_API_KEY")
OPENROUTER_BASE_URL = "https://openrouter.ai/api/v1"

# Site information
SITE_URL = os.environ.get("SITE_URL", "https://smiles2spec.app")
SITE_NAME = os.environ.get("SITE_NAME", "SMILES2Spec App")

# Model configurations
MODELS = {
    "chat": "anthropic/claude-3-haiku:free",  # For general chat
    "chemistry": "deepseek/deepseek-chat-v3-0324:free",  # For chemistry tasks
}

def get_headers():
    """Get standard headers for API requests."""
    if not OPENROUTER_API_KEY:
        raise ValueError("OPENROUTER_API_KEY is not set")
        
    return {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {OPENROUTER_API_KEY}",
        "HTTP-Referer": SITE_URL,
        "X-Title": SITE_NAME,
    }
```

#### `common/streaming.py`
```python
"""Streaming implementation for LLM responses."""
import json
import requests
from typing import Dict, Any, Generator, Optional, Callable
from .config import OPENROUTER_BASE_URL, get_headers

def stream_response(
    messages: list, 
    model: str,
    callback: Optional[Callable[[str], None]] = None
) -> Generator[str, None, None]:
    """Stream a response from the LLM."""
    url = f"{OPENROUTER_BASE_URL}/chat/completions"
    headers = get_headers()
    
    payload = {
        "model": model,
        "messages": messages,
        "stream": True
    }
    
    buffer = ""
    
    with requests.post(url, headers=headers, json=payload, stream=True) as r:
        if r.status_code != 200:
            error_msg = f"API request failed: {r.status_code} - {r.text}"
            if callback:
                callback(error_msg)
            yield error_msg
            return
            
        for chunk in r.iter_content(chunk_size=1024, decode_unicode=True):
            if not chunk:
                continue
                
            buffer += chunk
            
            while True:
                line_end = buffer.find('\n')
                if line_end == -1:
                    break
                    
                line = buffer[:line_end].strip()
                buffer = buffer[line_end + 1:]
                
                if line.startswith('data: '):
                    data = line[6:]
                    if data == '[DONE]':
                        return
                        
                    try:
                        data_obj = json.loads(data)
                        content = data_obj["choices"][0]["delta"].get("content", "")
                        if content:
                            if callback:
                                callback(content)
                            yield content
                    except json.JSONDecodeError:
                        # Skip invalid JSON (like comments)
                        pass


def complete_response(messages: list, model: str) -> str:
    """Get a complete response (non-streaming)."""
    url = f"{OPENROUTER_BASE_URL}/chat/completions"
    headers = get_headers()
    
    payload = {
        "model": model,
        "messages": messages
    }
    
    response = requests.post(url, headers=headers, json=payload)
    
    if response.status_code != 200:
        raise Exception(f"API request failed: {response.status_code} - {response.text}")
        
    result = response.json()
    
    return result["choices"][0]["message"]["content"]
```

### 2. SMILES Generator Implementation

#### `smiles_generator/prompts.py`
```python
"""Prompts for SMILES generation."""

SYSTEM_PROMPT = """You are a chemistry expert who specializes in generating valid SMILES strings.
Your task is to generate chemically valid SMILES strings based on text descriptions.

Guidelines:
1. Always generate valid SMILES strings that can be processed by RDKit
2. Output ONLY the SMILES strings, one per line, nothing else
3. Ensure each SMILES follows standard chemical notation rules
4. Do not include any explanations, headers, or other text
5. If a description is vague, generate compounds that generally match that category
"""

def get_user_prompt(description: str, count: int) -> str:
    """Generate the user prompt for SMILES generation."""
    return f"Generate {count} chemically valid SMILES strings for: {description}"
```

#### `smiles_generator/validator.py`
```python
"""Validation of generated SMILES strings."""
import re
from rdkit import Chem
from typing import List, Tuple

def validate_smiles(smiles_list: List[str]) -> List[Tuple[str, bool]]:
    """Validate a list of SMILES strings using RDKit."""
    results = []
    
    for smiles in smiles_list:
        # Basic cleanup
        cleaned = smiles.strip()
        
        # Skip empty lines
        if not cleaned:
            continue
            
        # Validate with RDKit
        mol = Chem.MolFromSmiles(cleaned)
        is_valid = mol is not None
        
        results.append((cleaned, is_valid))
        
    return results

def extract_smiles_from_text(text: str) -> List[str]:
    """Extract SMILES strings from the model's response text."""
    # Split by lines and clean up
    lines = [line.strip() for line in text.split('\n')]
    
    # Remove code blocks markers if present
    if '```' in text:
        # Extract content from code blocks
        code_blocks = re.findall(r'```(?:smiles)?\n(.*?)```', text, re.DOTALL)
        if code_blocks:
            lines = []
            for block in code_blocks:
                lines.extend([line.strip() for line in block.split('\n')])
    
    # Remove numbering patterns like "1. " or "1) "
    cleaned_lines = []
    for line in lines:
        # Skip empty lines
        if not line:
            continue
            
        # Remove numbering
        cleaned = re.sub(r'^\d+[\.\)]\s*', '', line)
        # Remove backticks if any
        cleaned = cleaned.replace('`', '')
        
        cleaned_lines.append(cleaned)
    
    return cleaned_lines
```

#### `smiles_generator/service.py`
```python
"""SMILES generation service implementation."""
import logging
from typing import List, Optional
from rdkit import Chem

from ..common.config import MODELS
from ..common.streaming import complete_response
from .prompts import SYSTEM_PROMPT, get_user_prompt
from .validator import extract_smiles_from_text, validate_smiles

logger = logging.getLogger(__name__)

def generate_smiles(
    description: str, 
    count: int = 1, 
    fallback: str = "C"
) -> List[str]:
    """Generate valid SMILES strings from a description."""
    try:
        # Create messages
        messages = [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": get_user_prompt(description, count)}
        ]
        
        # Get model response
        model = MODELS["chemistry"]
        logger.info(f"Requesting SMILES generation for: {description}")
        response_text = complete_response(messages, model)
        
        # Extract and validate SMILES
        smiles_candidates = extract_smiles_from_text(response_text)
        validated_results = validate_smiles(smiles_candidates)
        
        # Filter valid SMILES
        valid_smiles = [s for s, is_valid in validated_results if is_valid]
        
        if valid_smiles:
            # Limit to requested count
            return valid_smiles[:count]
        else:
            logger.warning(f"No valid SMILES generated for: {description}")
            return [fallback] * count
            
    except Exception as e:
        logger.error(f"SMILES generation error: {str(e)}")
        return [fallback] * count
```

### 3. Chat Assistant Implementation

#### `spectra_chat/prompts.py`
```python
"""Prompts for the Spectra Chat assistant."""

SYSTEM_PROMPT_TEMPLATE = """You are Spectra, a helpful AI assistant specializing in mass spectrometry and molecular chemistry.

You can help users understand:
- Mass spectrometry principles and interpretation
- Chemical structures and SMILES notation
- Molecular properties and fragmentation patterns
- Spectrum analysis and interpretation

{spectrum_context}

Keep your responses concise, accurate and helpful. If you don't know something, be honest about it.
"""

def get_system_prompt(spectrum_data=None):
    """Get the system prompt, optionally including spectrum information."""
    if spectrum_data:
        spectrum_context = f"""
Currently, you're analyzing a spectrum for the molecule:
- SMILES: {spectrum_data.get('smiles', 'Unknown')}
- Chemical name: {spectrum_data.get('chemical_name', 'Unknown')}
- Molecular weight: {spectrum_data.get('molecular_weight', 'Unknown')}
- Exact mass: {spectrum_data.get('exact_mass', 'Unknown')}

The spectrum has been predicted with our mass spectrometry model.
        """
    else:
        spectrum_context = "No specific spectrum is being analyzed currently."
        
    return SYSTEM_PROMPT_TEMPLATE.format(spectrum_context=spectrum_context)
```

#### `spectra_chat/formatter.py`
```python
"""Format spectrum data for inclusion in chat."""
import base64
from io import BytesIO
import matplotlib.pyplot as plt
import numpy as np

def spectrum_to_png_base64(spectrum_data):
    """Convert spectrum data to a PNG image encoded as base64."""
    # Extract spectrum data
    if not spectrum_data or 'spectrum' not in spectrum_data:
        return None
        
    try:
        x = spectrum_data['spectrum']['x']
        y = spectrum_data['spectrum']['y']
        
        # Create plot
        plt.figure(figsize=(10, 6))
        plt.bar(x, y, width=0.5, alpha=0.7)
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.title(f"Predicted Spectrum: {spectrum_data.get('chemical_name', spectrum_data.get('smiles', 'Unknown'))}")
        
        # Convert plot to PNG
        buf = BytesIO()
        plt.savefig(buf, format='png', dpi=100)
        plt.close()
        buf.seek(0)
        
        # Encode as base64
        png_base64 = base64.b64encode(buf.read()).decode('utf-8')
        return png_base64
        
    except Exception as e:
        print(f"Error creating spectrum image: {e}")
        return None

def create_spectrum_message(spectrum_data):
    """Create a message object with the spectrum image."""
    png_base64 = spectrum_to_png_base64(spectrum_data)
    
    if not png_base64:
        return None
        
    # Create message with image
    return {
        "role": "user",
        "content": [
            {
                "type": "text",
                "text": "Here is the predicted mass spectrum for the molecule:"
            },
            {
                "type": "image_url",
                "image_url": {
                    "url": f"data:image/png;base64,{png_base64}"
                }
            }
        ]
    }
```

#### `spectra_chat/service.py`
```python
"""Chat service for Spectra assistant."""
import logging
from typing import List, Dict, Any, Generator, Optional

from ..common.config import MODELS
from ..common.streaming import stream_response, complete_response
from .prompts import get_system_prompt
from .formatter import create_spectrum_message

logger = logging.getLogger(__name__)

def generate_chat_response(
    messages: List[Dict[str, str]],
    spectrum_data: Optional[Dict[str, Any]] = None,
    stream: bool = False,
    callback: Optional[callable] = None
) -> Any:
    """Generate a chat response from the Spectra assistant."""
    try:
        # Create system message with spectrum context if available
        system_prompt = get_system_prompt(spectrum_data)
        system_message = {"role": "system", "content": system_prompt}
        
        # Prepare full message list
        full_messages = [system_message]
        
        # Add spectrum image if available
        if spectrum_data:
            spectrum_message = create_spectrum_message(spectrum_data)
            if spectrum_message:
                full_messages.append(spectrum_message)
        
        # Add user messages
        full_messages.extend(messages)
        
        # Use the chat model
        model = MODELS["chat"]
        
        # Generate response
        if stream:
            return stream_response(full_messages, model, callback)
        else:
            return complete_response(full_messages, model)
            
    except Exception as e:
        error_msg = f"Chat response generation error: {str(e)}"
        logger.error(error_msg)
        if stream:
            if callback:
                callback(error_msg)
            return (yield error_msg)
        else:
            return error_msg
```

### 4. Main Entry Point

#### `llm_integration/__init__.py`
```python
"""LLM integration package for the mass spectrometry prediction API."""

# Import main functions for easy access
from .smiles_generator.service import generate_smiles
from .spectra_chat.service import generate_chat_response

__all__ = [
    'generate_smiles',
    'generate_chat_response'
]
```

### 5. API Endpoints Update

#### Update to `app.py` endpoints
```python
@app.route('/api/chat', methods=['POST'])
def chat():
    """API endpoint for chat functionality."""
    try:
        data = request.json
        messages = data.get('messages', [])
        stream = data.get('stream', False)
        smiles = data.get('smiles')
        
        if not messages:
            return jsonify({"error": "No messages provided"}), 400
        
        # Get spectrum data if SMILES provided
        spectrum_data = None
        if smiles:
            spectrum_data = prediction_service.predict_spectrum_from_smiles(smiles)
            
        if stream:
            def generate():
                for chunk in generate_chat_response(messages, spectrum_data, stream=True):
                    yield f"data: {json.dumps({'chunk': chunk})}\n\n"
                yield "data: [DONE]\n\n"
                
            return Response(generate(), mimetype='text/event-stream')
        else:
            message = generate_chat_response(messages, spectrum_data)
            return jsonify({"message": message})
            
    except Exception as e:
        logger.error(f"Chat API error: {str(e)}")
        return jsonify({"error": str(e)}), 500

@app.route('/api/generate_smiles', methods=['POST'])
def generate_smiles_endpoint():
    """API endpoint for generating SMILES strings."""
    try:
        data = request.json or {}
        count = max(1, min(10, data.get('count', 1)))  # Limit between 1-10
        description = data.get('description', 'organic molecule')
        
        smiles_list = generate_smiles(description, count)
        
        return jsonify({"smiles": smiles_list})
    except Exception as e:
        logger.error(f"SMILES generation error: {str(e)}")
        # Return a simple molecule if generation fails
        return jsonify({"smiles": ["C"] * max(1, min(10, request.json.get('count', 1)))}), 500
```

## Frontend Integration

### 1. SMILES Generator Component

Update the SMILES generator component to connect to the new backend:

```javascript
// src/components/SmilesGenerator.js
import { useState } from 'react';
import axios from 'axios';

export function SmilesGenerator() {
  const [description, setDescription] = useState('');
  const [count, setCount] = useState(1);
  const [smiles, setSmiles] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const generateSmiles = async () => {
    setLoading(true);
    setError(null);
    
    try {
      const response = await axios.post('/api/generate_smiles', {
        description,
        count: parseInt(count, 10)
      });
      
      setSmiles(response.data.smiles || []);
    } catch (err) {
      setError('Failed to generate SMILES: ' + (err.response?.data?.error || err.message));
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="smiles-generator">
      <h2>Generate SMILES</h2>
      
      <div className="form-group">
        <label>Description:</label>
        <input 
          type="text" 
          value={description} 
          onChange={(e) => setDescription(e.target.value)}
          placeholder="e.g., aromatic compounds with nitrogen"
        />
      </div>
      
      <div className="form-group">
        <label>Count:</label>
        <input 
          type="number" 
          min="1" 
          max="10" 
          value={count} 
          onChange={(e) => setCount(e.target.value)}
        />
      </div>
      
      <button 
        onClick={generateSmiles} 
        disabled={loading || !description}
      >
        {loading ? 'Generating...' : 'Generate SMILES'}
      </button>
      
      {error && (
        <div className="error-message">{error}</div>
      )}
      
      {smiles.length > 0 && (
        <div className="results">
          <h3>Generated SMILES</h3>
          <ul>
            {smiles.map((s, index) => (
              <li key={index}>
                <code>{s}</code>
                <button 
                  className="use-button"
                  onClick={() => {
                    // Add functionality to use this SMILES
                    window.dispatchEvent(new CustomEvent('use-smiles', { 
                      detail: { smiles: s } 
                    }));
                  }}
                >
                  Use
                </button>
              </li>
            ))}
          </ul>
        </div>
      )}
    </div>
  );
}
```

### 2. Chat with Spectra Component (with Streaming)

Implement the Chat with Spectra component with streaming support:

```javascript
// src/components/SpectraChat.js
import { useState, useEffect, useRef } from 'react';
import axios from 'axios';

export function SpectraChat({ smiles }) {
  const [messages, setMessages] = useState([]);
  const [input, setInput] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);
  const messagesEndRef = useRef(null);
  const eventSourceRef = useRef(null);

  useEffect(() => {
    // Scroll to bottom when messages change
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [messages]);

  // Clean up any open event sources on unmount
  useEffect(() => {
    return () => {
      if (eventSourceRef.current) {
        eventSourceRef.current.close();
      }
    };
  }, []);

  const sendMessage = async () => {
    if (!input.trim() || isLoading) return;
    
    const userMessage = { role: 'user', content: input };
    setMessages(prev => [...prev, userMessage]);
    setInput('');
    setIsLoading(true);
    setError(null);
    
    // Close any existing connection
    if (eventSourceRef.current) {
      eventSourceRef.current.close();
    }
    
    try {
      // Add the new user message to the history
      const chatHistory = [...messages, userMessage];
      
      // Prepare the API request data
      const requestData = {
        messages: chatHistory.map(msg => ({
          role: msg.role,
          content: msg.content
        })),
        stream: true
      };
      
      // Add SMILES if available
      if (smiles) {
        requestData.smiles = smiles;
      }
      
      // Start streaming
      const response = await fetch('/api/chat', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(requestData)
      });
      
      // Create placeholder for assistant response
      const assistantMessage = { role: 'assistant', content: '' };
      setMessages(prev => [...prev, assistantMessage]);
      
      // Handle streaming response
      const reader = response.body.getReader();
      const decoder = new TextDecoder();
      let buffer = '';
      
      while (true) {
        const { done, value } = await reader.read();
        if (done) break;
        
        // Decode the chunk and add to buffer
        buffer += decoder.decode(value, { stream: true });
        
        // Process complete lines from buffer
        while (true) {
          const lineEnd = buffer.indexOf('\n');
          if (lineEnd === -1) break;
          
          const line = buffer.slice(0, lineEnd).trim();
          buffer = buffer.slice(lineEnd + 1);
          
          if (line.startsWith('data: ')) {
            const data = line.slice(6);
            if (data === '[DONE]') break;
            
            try {
              const parsed = JSON.parse(data);
              const content = parsed.chunk;
              
              if (content) {
                // Update the last message with the new content
                setMessages(prev => {
                  const newMessages = [...prev];
                  const lastMsg = newMessages[newMessages.length - 1];
                  newMessages[newMessages.length - 1] = {
                    ...lastMsg,
                    content: lastMsg.content + content
                  };
                  return newMessages;
                });
              }
            } catch (e) {
              // Ignore invalid JSON
            }
          }
        }
      }
    } catch (err) {
      setError('Error: ' + (err.response?.data?.error || err.message));
      
      // Add error message to chat
      setMessages(prev => [
        ...prev.slice(0, -1), // Remove the loading message
        { role: 'assistant', content: 'Sorry, I encountered an error. Please try again.' }
      ]);
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <div className="spectra-chat">
      <div className="chat-messages">
        {messages.map((msg, idx) => (
          <div key={idx} className={`message ${msg.role}`}>
            <div className="message-content">
              {msg.content || (isLoading && idx === messages.length - 1 ? 'Thinking...' : '')}
            </div>
          </div>
        ))}
        <div ref={messagesEndRef} />
      </div>
      
      <div className="chat-input">
        <textarea
          value={input}
          onChange={(e) => setInput(e.target.value)}
          placeholder="Ask about chemistry or spectra..."
          onKeyDown={(e) => {
            if (e.key === 'Enter' && !e.shiftKey) {
              e.preventDefault();
              sendMessage();
            }
          }}
          disabled={isLoading}
        />
        <button 
          onClick={sendMessage}
          disabled={isLoading || !input.trim()}
        >
          {isLoading ? 'Sending...' : 'Send'}
        </button>
      </div>
      
      {error && (
        <div className="chat-error">{error}</div>
      )}
    </div>
  );
}
```

### 3. Main Application Integration

Update the main application to integrate both components:

```javascript
// src/components/AITools.js
import { useState } from 'react';
import { SmilesGenerator } from './SmilesGenerator';
import { SpectraChat } from './SpectraChat';

export function AITools({ currentSmiles }) {
  const [activeTab, setActiveTab] = useState('chat');
  const [selectedSmiles, setSelectedSmiles] = useState(currentSmiles || '');
  
  // Listen for SMILES selection events
  useEffect(() => {
    const handleUseSmiles = (event) => {
      setSelectedSmiles(event.detail.smiles);
    };
    
    window.addEventListener('use-smiles', handleUseSmiles);
    return () => {
      window.removeEventListener('use-smiles', handleUseSmiles);
    };
  }, []);
  
  // Update when current SMILES changes from parent
  useEffect(() => {
    setSelectedSmiles(currentSmiles);
  }, [currentSmiles]);

  return (
    <div className="ai-tools">
      <div className="tabs">
        <button 
          className={activeTab === 'chat' ? 'active' : ''}
          onClick={() => setActiveTab('chat')}
        >
          Chat with Spectra
        </button>
        <button 
          className={activeTab === 'generator' ? 'active' : ''}
          onClick={() => setActiveTab('generator')}
        >
          SMILES Generator
        </button>
      </div>
      
      <div className="tab-content">
        {activeTab === 'chat' ? (
          <SpectraChat smiles={selectedSmiles} />
        ) : (
          <SmilesGenerator />
        )}
      </div>
    </div>
  );
}
```

## Implementation Plan

### Phase 1: Backend Refactoring (Week 1)
1. Create the new directory structure
2. Implement common modules
3. Implement SMILES generator
4. Implement Chat with Spectra
5. Update API endpoints
6. Write unit tests for all components
7. Integration testing of backend services

### Phase 2: Frontend Integration (Week 2)
1. Create/update frontend components
2. Implement streaming functionality in chat interface
3. Add spectrum visualization
4. Connect frontend to backend APIs
5. UI/UX testing and refinement
6. Cross-browser compatibility testing

### Phase 3: Deployment (Week 3)
1. Staging environment deployment
2. Performance testing
3. Final bug fixes
4. Production deployment
5. Monitoring and feedback gathering

## Testing Strategy

### Backend Tests
- Unit tests for each module
- Integration tests for API endpoints
- Performance tests for streaming functionality

### Frontend Tests
- Component tests with React Testing Library
- End-to-end tests with Cypress
- Visual regression testing
- User acceptance testing

## Success Criteria
- Clear separation of SMILES generator and Chat with Spectra
- Streaming responses working in the chat interface
- Spectrum visualization integrated in chat
- Response times within acceptable limits (< 1s for SMILES, < 3s for chat initial response)
- Code follows KISS and SRP principles

## Recommendations
1. Use TypeScript for frontend implementation to improve type safety
2. Add a caching layer for spectrum predictions
3. Consider implementing retry mechanisms for failed API calls
4. Add comprehensive error handling and user-friendly error messages
5. Implement analytics to monitor usage patterns and performance