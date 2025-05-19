"""
Chat with Spectrum LLM integration service.
"""

import os
import requests
from openai import OpenAI
from ..utils import logger

def generate_chat_response(messages):
    """
    Generates a response for the chat functionality.
    
    Args:
        messages: List of message objects with role and content
        
    Returns:
        The generated message text
    """
    try:
        # Default system message to provide context about the application
        system_message = {
            "role": "system", 
            "content": """You are Spectra, an AI assistant specialized in mass spectrometry, chemistry, 
            and molecular structures. You can help with interpreting mass spectra, explaining SMILES 
            notation, and providing information about chemical compounds. Be helpful, clear, and 
            provide detailed explanations when discussing chemical compounds, SMILES notation, and 
            mass spectrometry data. When you don't know something, admit it rather than making up
            information."""
        }
        
        # Construct the final message list with the system message first
        final_messages = [system_message] + messages
        
        # Use OpenRouter integration
        client = OpenAI(
            base_url="https://openrouter.ai/api/v1",
            api_key=os.environ.get("OPENROUTER_API_KEY")
        )

        completion = client.chat.completions.create(
            extra_headers={
                "HTTP-Referer": os.environ.get("SITE_URL", "https://smiles2spec.app"),
                "X-Title": os.environ.get("SITE_NAME", "SMILES2Spec App"),
            },
            model="google/gemini-2.0-flash-exp:free",
            messages=final_messages
        )
        
        message = completion.choices[0].message.content
        
        # If API call fails, fall back to mock response
        if not message:
            user_message = next((msg["content"] for msg in messages if msg["role"] == "user"), "")
            message = generate_mock_response(user_message)
        
        return message
    except Exception as e:
        logger.error(f"Chat service error: {str(e)}")
        # Fall back to mock response in case of error
        user_message = next((msg["content"] for msg in messages if msg["role"] == "user"), "")
        return generate_mock_response(user_message)

def generate_mock_response(message):
    """
    Generates a mock response for testing purposes.
    In production, replace this with your actual AI model integration.
    """
    if "SMILES" in message or "smiles" in message:
        return "SMILES (Simplified Molecular Input Line Entry System) is a notation that allows you to represent chemical structures as text strings. It's used for storing and sharing molecular information. For example, 'CC' represents ethane, while 'C1=CC=CC=C1' represents benzene."
    elif "mass spec" in message.lower() or "spectrum" in message.lower():
        return "Mass spectrometry works by ionizing chemical samples and measuring the mass-to-charge ratio of the ions. The data is typically presented as a mass spectrum - a plot of ion signal vs. mass-to-charge ratio. The peaks in the spectrum represent different fragments of the original molecule."
    elif "hello" in message.lower() or "hi" in message.lower():
        return "Hello! I'm Spectra, your spectral analysis assistant. How can I help you today with mass spectrometry, SMILES notation, or chemical structures?"
    else:
        return "I understand you're asking about chemistry or spectral analysis. Could you provide more specific details about what you'd like to know? I can help with interpreting mass spectra, explaining SMILES notation, or discussing chemical structures." 