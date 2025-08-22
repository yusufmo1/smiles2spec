"""Chat routes."""
import json
from flask import Blueprint, request, jsonify, Response
from ...integrations.llm.services.chat_service import ChatService
from ...integrations.llm.services.smiles_service import SMILESService
from ...services.prediction_service import PredictionService
from ...utils import logger
from ..schemas.chat import ChatRequest, GenerateSMILESRequest
from ..middleware.validation import validate_json

chat_bp = Blueprint('chat', __name__, url_prefix='')
chat_service = ChatService()
smiles_service = SMILESService()
prediction_service = PredictionService()

@chat_bp.route('/chat', methods=['POST'])
@validate_json(ChatRequest)
def chat(validated_data: ChatRequest):
    """
    Chat with AI assistant about mass spectrometry and molecular structures.
    
    This endpoint provides an interactive chat interface powered by GPT-4, with optional
    context awareness of predicted mass spectra when a SMILES string is provided.
    
    Args:
        validated_data: ChatRequest containing:
            - messages: List of chat messages with role and content
            - smiles (optional): SMILES string for spectrum context
            - stream (optional): Whether to stream the response
            
    Returns:
        JSON response with:
            - message: AI assistant's response (non-streaming mode)
        OR
        Server-Sent Events stream (streaming mode)
        
    Raises:
        400: No valid messages provided
        500: Error generating response
    """
    try:
        # Convert Pydantic messages to dict format, filtering empty messages
        messages = []
        for msg in validated_data.messages:
            if msg.content.strip():  # Only include non-empty messages
                messages.append({"role": msg.role, "content": msg.content})
        
        if not messages:
            return jsonify({"error": "No valid messages provided"}), 400
        
        # Get spectrum context if SMILES provided
        spectrum_context = None
        if validated_data.smiles:
            try:
                spectrum_context = prediction_service.predict_spectrum(validated_data.smiles)
                logger.info(f"Generated spectrum context for SMILES: {validated_data.smiles}")
            except Exception as e:
                logger.warning(f"Could not get spectrum context: {str(e)}")
        
        if validated_data.stream:
            def generate():
                try:
                    for chunk in chat_service.generate_response(messages, spectrum_context, stream=True):
                        chunk_json = json.dumps({'chunk': chunk}, ensure_ascii=False)
                        yield f"data: {chunk_json}\n\n"
                    yield "data: [DONE]\n\n"
                except Exception as e:
                    logger.error(f"Error in streaming response: {str(e)}")
                    error_json = json.dumps({'chunk': f'Error: {str(e)}'}, ensure_ascii=False)
                    yield f"data: {error_json}\n\n"
                    yield "data: [DONE]\n\n"
            
            response = Response(generate(), mimetype='text/event-stream; charset=utf-8')
            response.headers['Cache-Control'] = 'no-cache'
            response.headers['Connection'] = 'keep-alive'
            return response
        else:
            message = chat_service.generate_response(messages, spectrum_context, stream=False)
            # Handle generator responses
            if hasattr(message, '__iter__') and not isinstance(message, (str, dict, list)):
                message = "".join(message)
            
            logger.debug(f"Chat response: {str(message)[:100]}...")
            return jsonify({"message": message})
            
    except Exception as e:
        logger.error(f"Chat error: {str(e)}")
        return jsonify({"error": str(e)}), 500

@chat_bp.route('/generate_smiles', methods=['POST'])
@validate_json(GenerateSMILESRequest)
def generate_smiles_route(validated_data: GenerateSMILESRequest):
    """
    Generate SMILES molecular structures from natural language descriptions.
    
    Uses GPT-4 to generate chemically valid SMILES strings based on user descriptions.
    Includes retry logic and fallback mechanisms to ensure valid results.
    
    Args:
        validated_data: GenerateSMILESRequest containing:
            - description: Natural language description of desired molecules
            - count: Number of SMILES strings to generate (default: 1)
            
    Returns:
        JSON response containing:
            - smiles: List of generated SMILES strings
            
    Notes:
        - Attempts generation up to 3 times with the original description
        - Falls back to enhanced prompts if initial attempts fail
        - Returns simple molecules (e.g., 'C') as last resort
        
    Raises:
        500: Error during generation (returns fallback molecules)
    """
    try:
        logger.info(f"Generate SMILES request: count={validated_data.count}, description='{validated_data.description}'")
        
        # Try up to 3 times to get valid results
        max_attempts = 3
        for attempt in range(max_attempts):
            logger.info(f"Attempt {attempt+1}/{max_attempts} for SMILES generation")
            smiles_list = smiles_service.generate_smiles(validated_data.description, validated_data.count)
            
            # If we got non-fallback results, return them
            if smiles_list and smiles_list[0] != "C":
                logger.info(f"SMILES generation successful (attempt {attempt+1}): {smiles_list}")
                return jsonify({"smiles": smiles_list})
            
            logger.warning(f"Got fallback SMILES on attempt {attempt+1}, retrying...")
        
        # If all attempts failed, try with enhanced description
        logger.warning("All attempts failed, trying enhanced prompt")
        enhanced_description = f"a {validated_data.description} with explicit chemical structure"
        smiles_list = smiles_service.generate_smiles(enhanced_description, validated_data.count)
        
        logger.info(f"SMILES generation with enhanced prompt: {smiles_list}")
        return jsonify({"smiles": smiles_list})
        
    except Exception as e:
        logger.error(f"SMILES generation error: {str(e)}")
        # Return fallback molecules if generation fails
        fallback_smiles = ["C"] * validated_data.count
        return jsonify({"smiles": fallback_smiles}), 500 