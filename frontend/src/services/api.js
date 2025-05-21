const API_URL = '/api'; 

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
 * @returns {Promise<ChatResponse>} The response from the API
 */
export async function chatWithSpectrum(messages) {
  const res = await fetch(`${API_URL}/api/chat`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ messages })
  });
  if (!res.ok) throw new Error("Chat request failed");
  return res.json();        // { message: "response text" }
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
 * Generates SMILES strings
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