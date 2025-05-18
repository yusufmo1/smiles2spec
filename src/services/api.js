const API_URL = 'http://localhost:5050';

export async function predictSpectrum(smiles) {
  const response = await fetch(`${API_URL}/api/predict`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
    },
    body: JSON.stringify({ smiles }),
  });
  
  if (!response.ok) {
    throw new Error('Failed to predict spectrum');
  }
  
  return response.json();
}

export async function exportMsp(smiles) {
  const res = await fetch(`${API_URL}/api/export_msp`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ smiles })
  });
  if (!res.ok) throw new Error("MSP export failed");
  return res.blob();                       // text/plain blob
}

export async function exportMspBatch(smilesList) {
  const res = await fetch(`${API_URL}/api/export_msp_batch`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ smiles_list: smilesList })
  });
  if (!res.ok) throw new Error("Batch MSP export failed");
  return res.blob();                       // text/plain blob
}

export async function uploadSmilesFile(file) {
  const form = new FormData();
  form.append("file", file);

  const res = await fetch(`${API_URL}/api/smiles_bulk`, {
    method: "POST",
    body: form
  });
  if (!res.ok) throw new Error("Bulk upload failed");
  return res.json();        // { smiles: [...] }
}

/**
 * Chat with the AI model using OpenRouter API
 * @param {Array} messages - Array of message objects {role, content}
 * @returns {AsyncIterableIterator<string>} - Stream of response chunks
 */
export async function chat(messages) {
  const apiKey = import.meta.env.VITE_OPENROUTER_API_KEY;
  if (!apiKey) {
    throw new Error('OpenRouter API key not found. Please set VITE_OPENROUTER_API_KEY in .env');
  }

  const response = await fetch("https://openrouter.ai/api/v1/chat/completions", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
      "Authorization": `Bearer ${apiKey}`,
      "HTTP-Referer": window.location.origin,
      "X-Title": "Spectra Chat"
    },
    body: JSON.stringify({
      model: "openai/gpt-4o-mini",
      messages,
      stream: true,
      temperature: 0.7
    })
  });

  if (!response.ok) {
    const error = await response.text();
    throw new Error(`Chat API error: ${error}`);
  }

  // Create a stream parser for SSE
  const reader = response.body.getReader();
  const decoder = new TextDecoder("utf-8");
  
  async function* streamMessages() {
    let buffer = '';
    
    while (true) {
      const { done, value } = await reader.read();
      if (done) break;
      
      buffer += decoder.decode(value, { stream: true });
      const lines = buffer.split('\n');
      buffer = lines.pop() || '';
      
      for (const line of lines) {
        if (line.startsWith('data: ')) {
          const data = line.slice(6);
          if (data === '[DONE]') return;
          
          try {
            const parsed = JSON.parse(data);
            const content = parsed.choices[0]?.delta?.content || '';
            if (content) yield content;
          } catch (e) {
            console.error('Error parsing SSE message:', e);
          }
        }
      }
    }
  }
  
  return streamMessages();
} 