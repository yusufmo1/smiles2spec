"""Basic helper for OpenRouter API requests."""
import requests
from typing import Any, Dict
from .config import OPENROUTER_BASE_URL, get_headers


def post_chat(messages, model, **kwargs) -> Dict[str, Any]:
    """Send a chat completion request to OpenRouter."""
    url = f"{OPENROUTER_BASE_URL}/chat/completions"
    headers = get_headers()
    payload = {"model": model, "messages": messages}
    payload.update(kwargs)
    response = requests.post(url, headers=headers, json=payload, timeout=30)
    response.raise_for_status()
    return response.json()
