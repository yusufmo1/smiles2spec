"""Data utility functions."""
import numpy as np
from typing import Any, Dict, List, Union

def convert_np_to_list(item: Any) -> Any:
    """Recursively convert numpy arrays to lists for JSON serialization."""
    if isinstance(item, np.ndarray):
        return item.tolist()
    elif isinstance(item, dict):
        return {k: convert_np_to_list(v) for k, v in item.items()}
    elif isinstance(item, list):
        return [convert_np_to_list(v) for v in item]
    elif isinstance(item, (np.integer, np.floating)):
        return item.item()
    else:
        return item

def ensure_numpy_array(data: Union[List, np.ndarray]) -> np.ndarray:
    """Ensure data is a numpy array."""
    if isinstance(data, list):
        return np.array(data)
    return data

def safe_float(value: Any, default: float = 0.0) -> float:
    """Safely convert value to float."""
    try:
        return float(value)
    except (ValueError, TypeError):
        return default

def safe_int(value: Any, default: int = 0) -> int:
    """Safely convert value to int."""
    try:
        return int(value)
    except (ValueError, TypeError):
        return default 