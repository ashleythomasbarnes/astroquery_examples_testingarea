import numpy as np
from astropy.table import MaskedColumn

def from_element_to_list(element, element_type=str):
    """
    Converts an input element into a list while ensuring all elements are of the specified type.

    This function ensures that:
    - If `element` is already a list, it is returned as is (after type validation).
    - If `element` is a NumPy array, it is converted to a list.
    - If `element` is an Astropy MaskedColumn, its data is extracted and converted to a list.
    - If `element` is a single value of the expected type, it is wrapped in a list.
    - If `element` is `None`, it returns `None`.

    Args:
        element (any): The input element that will be converted into a list.
        element_type (type, optional): The expected type of elements in the list. Defaults to `str`.

    Returns:
        list or None: A list containing the `element`, or `None` if `element` is `None`.
    """
    if element is None:
        return None

    if isinstance(element, list):
        assert all(isinstance(e, element_type) for e in element), f"All elements must be of type {element_type}"
        return element

    if isinstance(element, np.ndarray):
        element_list = element.tolist()
        assert all(isinstance(e, element_type) for e in element_list), f"All elements must be of type {element_type}"
        return element_list

    if isinstance(element, MaskedColumn):
        element_list = element.data.data.tolist()
        assert all(isinstance(e, element_type) for e in element_list), f"All elements must be of type {element_type}"
        return element_list

    if isinstance(element, element_type):
        return [element]

    raise TypeError(f"Invalid type for: {element} (expected {element_type})")


def from_bytes_to_string(input_in_bytes):
    """
    Converts an input from `bytes` to `str`. 

    This function handles:
    - A single bytes input (decodes to a string).
    - A NumPy array containing bytes (element-wise decoding).
    - A list containing bytes (element-wise decoding).
    - Any other input type remains unchanged.

    Args:
        input_in_bytes (any): The input value, which may be in bytes format.

    Returns:
        any: The input converted to a string (if applicable), otherwise returned unchanged.
    """
    if isinstance(input_in_bytes, bytes):
        return input_in_bytes.decode("utf-8")

    if isinstance(input_in_bytes, np.ndarray) and input_in_bytes.dtype.type is np.bytes_:
        return np.array([x.decode("utf-8") if isinstance(x, bytes) else x for x in input_in_bytes])

    if isinstance(input_in_bytes, list):
        return [x.decode("utf-8") if isinstance(x, bytes) else x for x in input_in_bytes]

    return input_in_bytes  # Return unchanged if not bytes-related