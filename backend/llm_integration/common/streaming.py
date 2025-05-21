"""Streaming implementation for LLM responses."""
import json
import requests
from typing import Dict, Any, Generator, Optional, Callable
from .config import OPENROUTER_BASE_URL, get_headers


def stream_response(
    messages: list,
    model: str,
    callback: Optional[Callable[[str], None]] = None
) -> Generator[str, None, None]:
    """Stream a response from the LLM."""
    url = f"{OPENROUTER_BASE_URL}/chat/completions"
    headers = get_headers()

    payload = {
        "model": model,
        "messages": messages,
        "stream": True
    }

    buffer = ""

    with requests.post(url, headers=headers, json=payload, stream=True) as r:
        if r.status_code != 200:
            error_msg = f"API request failed: {r.status_code} - {r.text}"
            if callback:
                callback(error_msg)
            yield error_msg
            return

        for chunk in r.iter_content(chunk_size=1024, decode_unicode=True):
            if not chunk:
                continue

            buffer += chunk

            while True:
                line_end = buffer.find('\n')
                if line_end == -1:
                    break

                line = buffer[:line_end].strip()
                buffer = buffer[line_end + 1:]

                if line.startswith('data: '):
                    data = line[6:]
                    if data == '[DONE]':
                        return

                    try:
                        data_obj = json.loads(data)
                        content = data_obj["choices"][0]["delta"].get("content", "")
                        if content:
                            if callback:
                                callback(content)
                            yield content
                    except json.JSONDecodeError:
                        # Skip invalid JSON (like comments)
                        pass


def complete_response(messages: list, model: str) -> str:
    """Get a complete response (non-streaming)."""
    url = f"{OPENROUTER_BASE_URL}/chat/completions"
    headers = get_headers()

    payload = {
        "model": model,
        "messages": messages
    }

    response = requests.post(url, headers=headers, json=payload)

    if response.status_code != 200:
        raise Exception(f"API request failed: {response.status_code} - {response.text}")

    result = response.json()

    return result["choices"][0]["message"]["content"]
