"""
Chat with Spectrum LLM integration service.
"""

import json
import requests
import re
from dotenv import load_dotenv
from ..utils import logger
from .llm_config import (
    get_model_for_task,
    get_system_prompt,
    get_openrouter_headers,
    OPENROUTER_BASE_URL,
)

load_dotenv()  # pick up OPENROUTER_API_KEY et al.


def _post_chat(messages, model):
    """Low-level helper. No streaming yet."""
    payload = {
        "model": model,
        "messages": messages,
    }
    headers = get_openrouter_headers()
    try:
        r = requests.post(
            f"{OPENROUTER_BASE_URL}/chat/completions",
            headers=headers,
            json=payload,
            timeout=30,
        )
        if r.status_code != 200:
            raise RuntimeError(f"{r.status_code}: {r.text[:200]}")
        data = r.json()
        # defensive checks
        choices = data.get("choices", [])
        if not choices or "message" not in choices[0]:
            raise RuntimeError("Unexpected OpenRouter reply shape")
        return choices[0]["message"]["content"]
    except Exception as exc:
        logger.error(f"LLM call failed: {exc}")
        raise


def generate_chat_response(history):
    """
    Generates a response for the chat functionality.
    
    Args:
        history: List of message objects with role and content
        
    Returns:
        The generated message text
    """
    system_prompt = {
        "role": "system",
        "content": get_system_prompt("general"),
    }
    messages = [system_prompt] + history
    model = get_model_for_task("general")

    try:
        answer = _post_chat(messages, model).strip()
        # Very light post-processing: collapse multiple newlines
        answer = re.sub(r"\n{3,}", "\n\n", answer)
        return answer
    except Exception:
        # Fall back to mock response in case of error
        user_message = next((msg["content"] for msg in history if msg["role"] == "user"), "")
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