import numpy as np
from astropy.table import MaskedColumn

def from_element_to_list(element, element_type=str):
    r"""Given an element it returns a list containing the element

    It also checks all the elements in the list have the same type defined by `element_type`

    Parameters
    ----------
        element (any): element that will be put in the list
        element_type (any): type of the element that should be contained in the list

    Returns
    -------
        list: list containing `element`

    """
    if element is None:
        return None
    elif isinstance(element, list):
        for element_in_list in element:
            assert isinstance(element_in_list, element_type), r'{} must be a {}'.format(element_in_list, element_type)
        return element
    elif isinstance(element, np.ndarray):
        element_list: list = element.tolist()
        for element_in_list in element_list:
            assert isinstance(element_in_list, element_type), r'{} must be a {}'.format(element_in_list, element_type)
        return element_list
    elif isinstance(element, MaskedColumn):
        element_list = element.data.data.tolist()
        for element_in_list in element_list:
            assert isinstance(element_in_list, element_type), r'{} must be a {}'.format(element_in_list, element_type)
        return element_list
    elif isinstance(element, element_type):
        return [element]
    else:
        print('Not valid type for: {}'.format(element))
    return


def from_bytes_to_string(input_in_bytes):
    r"""Given an input in `bytes` return it the corresponding `str`

    This is mainly to deal with the fact that TAP queries return a list in bytes format that might be annoying. If
    the input is not `bytes` nothing is changed.

    Args:
        input_in_bytes (any): input in bytes

    Returns:
        any: output converted into a string

    """
    # condition for a single entry as bytes
    if isinstance(input_in_bytes, bytes):
        output_as_str = str(input_in_bytes.decode("utf-8"))
    # condition for a np.array containing bytes
    elif isinstance(input_in_bytes, np.ndarray) and len(input_in_bytes.shape) == 1:
        output_as_str = np.empty_like(input_in_bytes)
        for idx in np.arange(0, np.size(output_as_str)):
            if isinstance(input_in_bytes[idx], bytes):
                output_as_str[idx] = str(input_in_bytes[idx].decode("utf-8"))
            else:
                output_as_str[idx] = input_in_bytes[idx]
    # condition for a list containing bytes
    elif isinstance(input_in_bytes, list):
        output_as_str = []
        for element_in_bytes in input_in_bytes:
            if isinstance(element_in_bytes, bytes):
                output_as_str.append(str(element_in_bytes.decode("utf-8")))
            else:
                output_as_str.append(element_in_bytes)
    # return copy of the entry
    else:
        output_as_str = input_in_bytes.copy()
    return output_as_str