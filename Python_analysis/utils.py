"""
These are helper functions that are frequently used
"""
import mat4py
import numpy as np



def get_struct_field(data, field, subfield=None):
    '''
    Extracts a field from the struct returned by scipy.io.loadmat
    :param data: data object returned by scip.io.loadmat
    :param field: field name
    :param subfield: subfield name
    :return: a numpy array corresponding to the extracted field
    '''
    if subfield is None:
        return np.squeeze(data[field][0])
    else:
        return data['data'][field][0, 0][subfield][0, 0][0]

def get_struct_field_mat4py(data, field, iscell=False, trialid=0):
    '''
    Extracts a field from the struct returned by mat4py.loadmat
    :param data: data object returned by mat4py.loadmat
    :param field: field name
    :param iscell: whether the extracted field is a matlab cell
    :param trialid: id of cell to extract (if iscell=True)
    :return: a np array corresponding to the extracted field
    '''
    if not iscell:
        return np.array(data[field])
    else:
        temp = data[field][trialid][0]
        return np.array(temp)


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

def get_multiple_struct_fields_mat4py(data, fields):
    '''
    Get multiple fields from the object returned by mat4py.loadmat
    :param data: The data structure to read
    :param fields: A list of strings representing fields
    :return: A tuple of extracted fields
    '''
    extracted_fields = []
    for field in fields:
        extracted = get_struct_field_mat4py(data, field)
        extracted_fields.append(extracted)
    return extracted_fields
