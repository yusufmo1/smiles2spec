const API_URL = '/api'; 
//const API_URL = 'http://192.168.1.211:5050'; 

/**
 * @typedef {Object} SpectrumPredictionResponse
 * @property {number[]} x - The x-axis values (m/z values)
 * @property {number[]} y - The y-axis values (intensity values)
 * @property {string} [name] - The chemical name
 * @property {Array<{mz: number, intensity: number}>} [peaks] - Peak data
 * @property {string} [png] - Structure PNG data
 */

/**
 * Predicts a mass spectrum from a SMILES string
 * @param {string} smiles - The SMILES string to predict from
 * @returns {Promise<SpectrumPredictionResponse>} The predicted spectrum data
 */
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

/**
 * Exports a SMILES structure as MSP format
 * @param {string} smiles - The SMILES string to export
 * @returns {Promise<Blob>} MSP data as text/plain blob
 */
export async function exportMsp(smiles) {
  const res = await fetch(`${API_URL}/api/export_msp`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ smiles })
  });
  if (!res.ok) throw new Error("MSP export failed");
  return res.blob();                       // text/plain blob
}

/**
 * Exports multiple SMILES structures as MSP format
 * @param {string[]} smilesList - List of SMILES strings to export
 * @returns {Promise<Blob>} MSP data as text/plain blob
 */
export async function exportMspBatch(smilesList) {
  const res = await fetch(`${API_URL}/api/export_msp_batch`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ smiles_list: smilesList })
  });
  if (!res.ok) throw new Error("Batch MSP export failed");
  return res.blob();                       // text/plain blob
}

/**
 * @typedef {Object} BulkUploadResponse
 * @property {string[]} smiles - Array of SMILES strings
 */

/**
 * Uploads a file containing SMILES structures
 * @param {File} file - The file to upload
 * @returns {Promise<BulkUploadResponse>} The parsed SMILES strings
 */
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
 * @typedef {Object} Message
 * @property {string} role - The role (user or assistant)
 * @property {string} content - The message content
 */

/**
 * @typedef {Object} ChatResponse
 * @property {string} message - The assistant's response
 */

/**
 * Sends messages to the chat API
 * @param {Message[]} messages - Array of message objects with role and content
 * @param {Object} [options] - Additional options for the chat
 * @param {string} [options.smiles] - SMILES string to include with the chat
 * @param {function} [options.onChunk] - Callback for streaming chunks
 * @returns {Promise<void>} A promise that resolves when the stream is complete
 */
export async function chatWithSpectrum(messages, options = {}) {
  const { smiles, onChunk } = options;
  
  const payload = { 
    messages,
    stream: true
  };
  
  // Include SMILES if provided for spectrum visualization
  if (smiles) {
    payload.smiles = smiles;
  }
  
  try {
    const response = await fetch(`${API_URL}/api/chat`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload)
    });
    
    if (!response.ok) {
      throw new Error(`Chat request failed: ${response.status}`);
    }
    
    // Setup event stream processing
    const reader = response.body.getReader();
    const decoder = new TextDecoder();
    let buffer = '';
    
    // Return a promise that resolves when the stream is complete
    return new Promise((resolve, reject) => {
      async function processStream() {
        try {
          while (true) {
            const { done, value } = await reader.read();
            
            if (done) {
              // Flush any remaining content in buffer
              if (buffer.trim() && typeof onChunk === 'function') {
                try {
                  const parsed = JSON.parse(buffer.trim());
                  if (parsed.chunk) onChunk(parsed.chunk);
                } catch (e) {
                  // Skip invalid JSON
                }
              }
              resolve();
              break;
            }
            
            // Decode chunk and add to buffer
            buffer += decoder.decode(value, { stream: true });
            
            // Process complete lines from buffer
            let lineEnd;
            while ((lineEnd = buffer.indexOf('\n')) !== -1) {
              const line = buffer.slice(0, lineEnd).trim();
              buffer = buffer.slice(lineEnd + 1);
              
              if (line.startsWith('data: ')) {
                const data = line.slice(6);
                
                if (data === '[DONE]') {
                  resolve();
                  break;
                }
                
                try {
                  const parsed = JSON.parse(data);
                  if (parsed.chunk && typeof onChunk === 'function') {
                    onChunk(parsed.chunk);
                  }
                } catch (e) {
                  // Skip invalid JSON
                }
              }
            }
          }
        } catch (error) {
          reject(error);
        }
      }
      
      processStream();
    });
  } catch (error) {
    console.error("Streaming error:", error);
    throw error;
  }
}

/**
 * @typedef {Object} GenerateSmilesOptions
 * @property {number} [count=1] - Number of SMILES to generate
 * @property {string} [description=''] - Description to use for generation
 */

/**
 * @typedef {Object} GenerateSmilesResponse
 * @property {string[]} smiles - Generated SMILES strings
 */

/**
 * Generates SMILES strings using LLM
 * @param {GenerateSmilesOptions} [options={}] - Options for generation
 * @returns {Promise<GenerateSmilesResponse>} Generated SMILES strings
 */
export async function generateSmiles(options = {}) {
  const { count = 1, description = '' } = options;
  
  const response = await fetch(`${API_URL}/api/generate_smiles`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
    },
    body: JSON.stringify({ count, description }),
  });
  
  if (!response.ok) {
    throw new Error('Failed to generate SMILES');
  }
  
  return response.json(); // { smiles: ["...", "..."] }
}