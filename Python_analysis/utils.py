"""
These are helper functions that are frequently used
"""
import mat4py
import numpy as np




def get_struct_field(data, field, subfield=None):
    if subfield is None:
        return np.squeeze(data[field][0])
    else:
        return data['data'][field][0, 0][subfield][0, 0][0]

def get_multiple_struct_fields(data, fields):
    """
    :param data: The data structure to read
    :param fields: A list of strings representing fields
    :return: A tuple of extracted fields
    """
    extracted_fields = []
    for field in fields:
        extracted = get_struct_field(data, field)
        extracted_fields.append(extracted)
    return extracted_fields
