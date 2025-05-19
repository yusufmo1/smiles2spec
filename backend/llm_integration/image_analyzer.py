"""
Image analysis service using LLM with vision capabilities.
"""

import os
import base64
import requests
from openai import OpenAI
from ..utils import logger

def analyze_molecule_image(image_path):
    """
    Analyze a molecule image and extract information about the structure.
    
    Args:
        image_path: Path to the image file
        
    Returns:
        Dictionary with analysis results (description, possible SMILES, etc.)
    """
    try:
        api_key = os.environ.get('OPENROUTER_API_KEY')
        if not api_key or not os.path.exists(image_path):
            return {"error": "API key missing or image not found"}
            
        # Initialize OpenRouter client
        client = OpenAI(
            base_url="https://openrouter.ai/api/v1",
            api_key=api_key
        )
        
        # Convert image to base64 if needed
        with open(image_path, "rb") as image_file:
            image_data = image_file.read()
            
        # Prepare image for API
        if os.path.splitext(image_path)[1].lower() in ['.jpg', '.jpeg', '.png']:
            image_url = encode_image_to_data_url(image_path)
        else:
            return {"error": "Unsupported image format"}
        
        # Create message with the image
        messages = [
            {
                "role": "user",
                "content": [
                    {
                        "type": "text",
                        "text": "Analyze this molecule image. What is the chemical structure? If possible, provide the SMILES notation and name of the compound."
                    },
                    {
                        "type": "image_url",
                        "image_url": {
                            "url": image_url
                        }
                    }
                ]
            }
        ]
        
        # Call OpenRouter API with image
        completion = client.chat.completions.create(
            extra_headers={
                "HTTP-Referer": os.environ.get("SITE_URL", "https://smiles2spec.app"),
                "X-Title": os.environ.get("SITE_NAME", "SMILES2Spec App"),
            },
            model="google/gemini-2.0-flash-exp:free",  # Vision-capable model
            messages=messages
        )
        
        response_text = completion.choices[0].message.content
        
        # Parse the response to extract information
        result = {
            "description": response_text,
            "possible_smiles": extract_smiles_from_text(response_text),
        }
        
        return result
        
    except Exception as e:
        logger.error(f"Image analysis error: {str(e)}")
        return {"error": str(e)}

def encode_image_to_data_url(image_path):
    """
    Encode an image file to a data URL format that can be used in API requests.
    
    Args:
        image_path: Path to the image file
        
    Returns:
        Data URL string
    """
    with open(image_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    
    # Get file extension
    ext = os.path.splitext(image_path)[1].lower()
    mime_type = 'image/jpeg' if ext in ['.jpg', '.jpeg'] else 'image/png'
    
    return f"data:{mime_type};base64,{encoded_string}"

def extract_smiles_from_text(text):
    """
    Extract SMILES notation from text if present.
    
    Args:
        text: The text to parse
        
    Returns:
        List of possible SMILES strings or empty list if none found
    """
    import re
    
    # Basic regex pattern to find potential SMILES strings
    # This is a simplified pattern and might need refinement
    smiles_pattern = r'(?:SMILES|SMILES notation|notation|string)(?:[:\s]+)([A-Za-z0-9\(\)\[\]\.\=\#\-\+\@\*\/\\]+)'
    
    matches = re.findall(smiles_pattern, text, re.IGNORECASE)
    
    # Clean up matches and remove obvious non-SMILES strings
    cleaned_matches = []
    for match in matches:
        # Remove trailing punctuation and whitespace
        cleaned = match.strip('.,;: \t\n')
        if cleaned and len(cleaned) > 1:  # At least 2 characters
            cleaned_matches.append(cleaned)
    
    return cleaned_matches 