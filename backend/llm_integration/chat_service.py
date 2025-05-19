"""
Chat with Spectrum LLM integration service.
"""

import os
import requests
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
        
        # If you're using OpenAI or similar API directly, uncomment and use the code below
        # response = requests.post(
        #     "https://api.openai.com/v1/chat/completions",
        #     headers={
        #         "Content-Type": "application/json",
        #         "Authorization": f"Bearer {os.environ.get('OPENAI_API_KEY')}"
        #     },
        #     json={
        #         "model": "gpt-3.5-turbo",
        #         "messages": final_messages,
        #         "temperature": 0.7
        #     }
        # )
        # 
        # if response.status_code != 200:
        #     return "Sorry, I encountered an error while processing your request."
        # 
        # result = response.json()
        # message = result["choices"][0]["message"]["content"]
        
        # For now, use a mock response for testing
        user_message = next((msg["content"] for msg in messages if msg["role"] == "user"), "")
        message = generate_mock_response(user_message)
        
        return message
    except Exception as e:
        logger.error(f"Chat service error: {str(e)}")
        return "Sorry, I encountered an error while processing your request."

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