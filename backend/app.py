"""
Flask API for mass spectrometry prediction service.
"""

from flask import Flask, request, jsonify, render_template
from flask_cors import CORS
import os
import json
from prediction_service import PredictionService
from config import API_CONFIG
from utils import logger, convert_np_to_list, smiles_to_png_base64

# Initialize Flask app
app = Flask(__name__)
CORS(app)

# Initialize prediction service
prediction_service = PredictionService()

@app.route('/api/predict', methods=['POST'])
def predict():
    """API endpoint for spectrum prediction."""
    try:
        data = request.json
        smiles = data.get('smiles', '')
        
        if not smiles:
            return jsonify({"error": "No SMILES provided"}), 400
        
        result = prediction_service.predict_spectrum_from_smiles(smiles)
        
        if "error" in result:
            return jsonify(result), 400
        
        return jsonify(result)
    except Exception as e:
        logger.error(f"API error: {str(e)}")
        return jsonify({"error": str(e)}), 500

@app.route('/api/export_msp', methods=['POST'])
def export_msp():
    """
    Returns a downloadable .msp for the requested SMILES.
    Body: { "smiles": "CCC..." }
    """
    from flask import make_response
    from utils import peaks_to_msp
    
    data = request.json or {}
    smiles = data.get("smiles", "")
    if not smiles:
        return jsonify({"error": "No SMILES provided"}), 400

    pred = prediction_service.predict_spectrum_from_smiles(smiles)
    if "error" in pred:
        return jsonify(pred), 400

    msp_txt, fname = peaks_to_msp(pred)
    resp = make_response(msp_txt)
    resp.headers["Content-Type"] = "text/plain; charset=utf-8"
    resp.headers["Content-Disposition"] = f'attachment; filename="{fname}.msp"'
    return resp

@app.route('/api/export_msp_batch', methods=['POST'])
def export_msp_batch():
    """
    Returns a downloadable .msp file containing multiple entries.
    Body: { "smiles_list": ["CCC...", "CCCC...", ...] }
    """
    from flask import make_response
    from utils import peaks_to_msp
    import time
    
    data = request.json or {}
    smiles_list = data.get("smiles_list", [])
    if not smiles_list:
        return jsonify({"error": "No SMILES provided"}), 400
    
    # Generate MSP for each SMILES
    msp_entries = []
    timestamp = int(time.time())
    
    for smiles in smiles_list:
        try:
            pred = prediction_service.predict_spectrum_from_smiles(smiles)
            if "error" not in pred:
                msp_txt, _ = peaks_to_msp(pred)
                msp_entries.append(msp_txt)
        except Exception as e:
            logger.error(f"Error processing SMILES {smiles}: {str(e)}")
            # Continue with the next one even if this one fails
    
    # Combine all entries
    full_msp = "\n\n".join(msp_entries)
    
    resp = make_response(full_msp)
    resp.headers["Content-Type"] = "text/plain; charset=utf-8"
    resp.headers["Content-Disposition"] = f'attachment; filename="batch_export_{timestamp}.msp"'
    return resp

@app.route('/api/structure', methods=['GET'])
def get_structure():
    """API endpoint for getting molecular structure as PNG."""
    try:
        smiles = request.args.get('smiles', '')
        
        if not smiles:
            return jsonify({"error": "No SMILES provided"}), 400
        
        png_base64 = smiles_to_png_base64(smiles)
        
        if not png_base64:
            return jsonify({"error": "Invalid SMILES or unable to generate structure"}), 400
            
        return jsonify({"png_base64": png_base64})
    except Exception as e:
        logger.error(f"Structure API error: {str(e)}")
        return jsonify({"error": str(e)}), 500

@app.route('/api/structure/png', methods=['GET'])
def get_structure_png():
    """API endpoint for getting molecular structure as PNG."""
    try:
        smiles = request.args.get('smiles', '')
        
        if not smiles:
            return jsonify({"error": "No SMILES provided"}), 400
        
        png_base64 = smiles_to_png_base64(smiles)
        
        if not png_base64:
            return jsonify({"error": "Invalid SMILES or unable to generate structure"}), 400
            
        return jsonify({"png_base64": png_base64})
    except Exception as e:
        logger.error(f"PNG Structure API error: {str(e)}")
        return jsonify({"error": str(e)}), 500

@app.route('/api/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    status = {
        "status": "healthy",
        "model_loaded": prediction_service.model_loaded
    }
    return jsonify(status)

@app.route('/api/smiles_bulk', methods=['POST'])
def smiles_bulk():
    """
    Accepts a multipart/form-data file (txt or csv, one SMILES per line).
    Returns { smiles: [ "...", ... ] }
    """
    import csv
    from io import StringIO
    
    f = request.files.get("file")
    if f is None:
        return jsonify({"error": "no file"}), 400
    
    content = f.stream.read().decode("utf-8", errors="ignore")
    filename = f.filename.lower()
    
    # Process based on file type
    if filename.endswith('.csv'):
        # Parse CSV content
        smiles_list = []
        reader = csv.reader(StringIO(content))
        for row in reader:
            if row and row[0].strip():  # Use first column for SMILES
                smiles_list.append(row[0].strip())
        
        # Skip header row if it doesn't look like a SMILES string
        if smiles_list and not smiles_list[0].replace('[', '').replace(']', '').replace('(', '').replace(')', '').replace('.', '').strip().isalnum():
            smiles_list = smiles_list[1:]
        
        return jsonify({"smiles": smiles_list})
    else:
        # Default TXT processing - one SMILES per line
        smiles = [ln.split()[0].strip() for ln in content.splitlines() if ln.strip()]
        return jsonify({"smiles": smiles})

@app.route('/', methods=['GET'])
def home():
    """Simple frontend for testing."""
    return render_template('index.html')

# Add a simple HTML template for testing
os.makedirs('templates', exist_ok=True)
with open('templates/index.html', 'w') as f:
    f.write("""
<!DOCTYPE html>
<html>
<head>
    <title>Mass Spectrometry Prediction</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        #container { max-width: 800px; margin: 0 auto; }
        .input-group { margin-bottom: 15px; }
        label { display: block; margin-bottom: 5px; }
        input[type="text"] { width: 80%; padding: 8px; }
        button { padding: 8px 15px; background: #4CAF50; color: white; border: none; cursor: pointer; }
        #plot { width: 100%; height: 400px; }
        #results { margin-top: 20px; }
    </style>
</head>
<body>
    <div id="container">
        <h1>Mass Spectrometry Prediction</h1>
        <div class="input-group">
            <label for="smiles">SMILES String:</label>
            <input type="text" id="smiles" placeholder="Enter SMILES string" value="CC(=O)OC1=CC=CC=C1C(=O)O">
            <button onclick="predict()">Predict</button>
        </div>
        <div id="plot"></div>
        <div id="results"></div>
    </div>
    
    <script>
        function predict() {
            const smiles = document.getElementById('smiles').value;
            if (!smiles) {
                alert('Please enter a SMILES string');
                return;
            }
            
            fetch('/api/predict', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({ smiles: smiles })
            })
            .then(response => response.json())
            .then(data => {
                if (data.error) {
                    alert('Error: ' + data.error);
                    return;
                }
                
                // Create plot
                const layout = {
                    title: 'Predicted Mass Spectrum',
                    xaxis: { title: 'm/z' },
                    yaxis: { title: 'Relative Intensity' }
                };
                
                Plotly.newPlot('plot', [{
                    x: data.spectrum.x,
                    y: data.spectrum.y,
                    type: 'bar'
                }], layout);
                
                // Display results
                const resultsDiv = document.getElementById('results');
                resultsDiv.innerHTML = `
                    <h3>Results</h3>
                    <p><strong>SMILES:</strong> ${data.smiles}</p>
                    <p><strong>Molecular Weight:</strong> ${data.molecular_weight.toFixed(4)}</p>
                    <p><strong>Exact Mass:</strong> ${data.exact_mass.toFixed(4)}</p>
                    <p><strong>Number of Peaks:</strong> ${data.peaks.length}</p>
                `;
                
                if (data.warning) {
                    resultsDiv.innerHTML += `<p><strong>Warning:</strong> ${data.warning}</p>`;
                }
            })
            .catch(error => {
                console.error('Error:', error);
                alert('Error communicating with the server');
            });
        }
    </script>
</body>
</html>
    """)

if __name__ == '__main__':
    app.run(
        host=API_CONFIG.get('host', '0.0.0.0'),
        port=API_CONFIG.get('port', 5050),
        debug=API_CONFIG.get('debug', True)
    )