"""Shared configuration for LLM integration."""
import os
from dotenv import load_dotenv

load_dotenv()

OPENROUTER_API_KEY = os.environ.get("OPENROUTER_API_KEY")
OPENROUTER_BASE_URL = "https://openrouter.ai/api/v1"

SITE_URL = os.environ.get("SITE_URL", "https://smiles2spec.app")
SITE_NAME = os.environ.get("SITE_NAME", "SMILES2Spec App")

MODELS = {
    "chat": "anthropic/claude-3-haiku:free",
    "chemistry": "deepseek/deepseek-chat-v3-0324:free",
}

def get_headers():
    """Get standard headers for API requests."""
    if not OPENROUTER_API_KEY:
        raise ValueError("OPENROUTER_API_KEY is not set")
    return {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {OPENROUTER_API_KEY}",
        "HTTP-Referer": SITE_URL,
        "X-Title": SITE_NAME,
    }
