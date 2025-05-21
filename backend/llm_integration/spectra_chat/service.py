"""Chat service for Spectra assistant."""
import logging
from typing import List, Dict, Any, Generator, Optional

from ..common.config import MODELS
from ..common.streaming import stream_response, complete_response
from .prompts import get_system_prompt
from .formatter import create_spectrum_message

logger = logging.getLogger(__name__)


def generate_chat_response(
    messages: List[Dict[str, str]],
    spectrum_data: Optional[Dict[str, Any]] = None,
    stream: bool = False,
    callback: Optional[callable] = None
) -> Any:
    """Generate a chat response from the Spectra assistant."""
    try:
        # Create system message with spectrum context if available
        system_prompt = get_system_prompt(spectrum_data)
        system_message = {"role": "system", "content": system_prompt}

        # Prepare full message list
        full_messages = [system_message]

        # Add spectrum image if available
        if spectrum_data:
            spectrum_message = create_spectrum_message(spectrum_data)
            if spectrum_message:
                full_messages.append(spectrum_message)

        # Add user messages
        full_messages.extend(messages)

        # Use the chat model
        model = MODELS["chat"]

        # Generate response
        if stream:
            return stream_response(full_messages, model, callback)
        else:
            return complete_response(full_messages, model)

    except Exception as e:
        error_msg = f"Chat response generation error: {str(e)}"
        logger.error(error_msg)
        if stream:
            if callback:
                callback(error_msg)
            return (yield error_msg)
        else:
            return error_msg
