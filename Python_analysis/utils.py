"""
These are helper functions that are frequently used
"""
import mat4py
import numpy as np
import matplotlib.pyplot as plt


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


# Plot only for individual sessions
def plot_sessions(neuron_group):
    sessions = neuron_group.get_session_list()
    for i in range(np.max(sessions)):
        to_plot_list = np.where(sessions == i)
        mean_act = neuron_group.plot_all_means(sort=True, plotid=to_plot_list, style='heatmap')
        #plt.close('all')


# Compare neurons
def plot_neurons(subgroup1, subgroup2):
    session_id = 0
    for i in range(subgroup1.n_neurons):
        plt.figure()
        plt.subplot('121')
        subgroup1.neurons[i].plot_all_trials()
        plt.subplot('122')
        subgroup2.neurons[i].plot_all_trials()


# Decide if a neuron is c or i-responsive
def plot_example_neurons(neuron_group, id_list):
    """
    Plot the activity of the neurons specified in id_list, from neuron_group
    :param neuron_group: a NeuronGroup object
    :param id_list: a list of neurons to plot
    :return: nothing
    """
    dim = int(np.sqrt(len(id_list)))
    for i in id_list:
        plt.subplot(dim, dim, i + 1)
        neuron = neuron_group.neurons[i]
        plt.errorbar(np.arange(len(neuron.mean_activity)), neuron.mean_activity, neuron.stderr_activity)
        plt.title(str(i))

