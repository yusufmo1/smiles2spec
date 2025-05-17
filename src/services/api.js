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