"""
Flask API for mass spectrometry prediction service.
"""

from flask import Flask, request, jsonify, render_template, Response
from flask_cors import CORS
import os
import json
import requests
from .prediction_service import PredictionService
from .config import API_CONFIG
from .utils import logger, convert_np_to_list, smiles_to_png_base64
from .llm_integration import generate_chat_response, generate_smiles
from .model_downloader import initialize_models

# Try to download models if they don't exist
try:
    initialize_models()
except Exception as e:
    logger.error(f"Failed to download models: {str(e)}")
    logger.warning("Continuing without model download. Some features may not work correctly.")

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
    from .utils import peaks_to_msp
    
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
    from .utils import peaks_to_msp
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

@app.route('/api/chat', methods=['POST'])
def chat():
    """
    API endpoint for chat functionality.
    
    Accepts:
    {
        "messages": [{"role": "user"|"assistant", "content": "message"}],
        "stream": true|false,
        "smiles": "CCC" (optional)
    }
    """
    try:
        data = request.json
        messages = data.get('messages', [])
        stream = data.get('stream', False)
        smiles = data.get('smiles')
        
        logger.debug(f"Chat API request: stream={stream}, smiles={smiles}, messages_count={len(messages)}")
        
        if not messages:
            return jsonify({"error": "No messages provided"}), 400
        
        # Get spectrum data if SMILES provided
        spectrum_data = None
        if smiles:
            try:
                spectrum_data = prediction_service.predict_spectrum_from_smiles(smiles)
                logger.info(f"Generated spectrum data for SMILES: {smiles}")
            except Exception as e:
                logger.error(f"Error generating spectrum data: {str(e)}", exc_info=True)
                # Continue without spectrum data
            
        if stream:
            def generate():
                try:
                    for chunk in generate_chat_response(messages, spectrum_data, stream=True):
                        yield f"data: {json.dumps({'chunk': chunk})}\n\n"
                    yield "data: [DONE]\n\n"
                except Exception as e:
                    logger.error(f"Error in streaming response: {str(e)}", exc_info=True)
                    yield f"data: {json.dumps({'chunk': f'Error: {str(e)}'})}\n\n"
                    yield "data: [DONE]\n\n"
                
            return Response(generate(), mimetype='text/event-stream')
        else:
            # For non-streaming requests, we need to get the full response as a string
            message = generate_chat_response(messages, spectrum_data, stream=False)
            # If it's a generator (which it shouldn't be for non-streaming), get the content
            if hasattr(message, '__iter__') and not isinstance(message, (str, dict, list)):
                message = "".join(message)
            logger.debug(f"Chat API response: {message[:50]}...")
            return jsonify({"message": message})
            
    except Exception as e:
        logger.error(f"Chat API error: {str(e)}", exc_info=True)
        return jsonify({"error": str(e)}), 500

@app.route('/api/generate_smiles', methods=['POST'])
def generate_smiles():
    """API endpoint for generating SMILES strings."""
    try:
        data = request.json or {}
        count = max(1, min(10, data.get('count', 1)))  # Limit between 1-10
        description = data.get('description', '')
        
        logger.info(f"Generate SMILES API request: count={count}, description='{description}'")
        
        # Import here to avoid circular imports
        from .llm_integration import generate_smiles as generate_smiles_func
        
        # Try up to 3 times to get valid results
        max_attempts = 3
        for attempt in range(max_attempts):
            logger.info(f"Attempt {attempt+1}/{max_attempts} for SMILES generation")
            smiles_list = generate_smiles_func(count=count, description=description)
            
            # If we got non-fallback results, return them
            if smiles_list and smiles_list[0] != "C":
                logger.info(f"Generate SMILES API response (attempt {attempt+1}): {smiles_list}")
                return jsonify({"smiles": smiles_list})
            
            logger.warning(f"Got fallback SMILES on attempt {attempt+1}, retrying...")
        
        # If all attempts failed, try a more specific description
        logger.warning("All attempts failed with original description, trying with more specific prompt")
        enhanced_description = f"a {description} with explicit chemical structure"
        smiles_list = generate_smiles_func(count=count, description=enhanced_description)
        
        logger.info(f"Generate SMILES API response (enhanced prompt): {smiles_list}")
        return jsonify({"smiles": smiles_list})
        
    except Exception as e:
        logger.error(f"SMILES generation error: {str(e)}", exc_info=True)
        # Return a simple molecule if generation fails
        return jsonify({"smiles": ["C"] * max(1, min(10, request.json.get('count', 1)))}), 500

@app.route('/', methods=['GET'])
def home():
    """Simple frontend for testing."""
    return render_template('index.html')

if __name__ == '__main__':
    app.run(
        host=API_CONFIG.get('host', '0.0.0.0'),
        port=API_CONFIG.get('port', 5050),
        debug=API_CONFIG.get('debug', True)
    )