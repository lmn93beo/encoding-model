import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import pandas as pd
import utils
import mat4py
import data_classes

# Load data
raw_behavior_summary = mat4py.loadmat('TB41_behavior_summary.mat')
raw_encoding_struct = mat4py.loadmat('TB41_encoding_structs.mat')

# Extract structures from the raw data
one, two, three, four, correct, incorrect, left, right, prev_right, prev_wrong, \
    stim_onset_per_trial = utils.get_multiple_struct_fields_mat4py(raw_behavior_summary,
                        ['one', 'two', 'three', 'four', 'correct', 'incorrect',
                         'left', 'right', 'prev_right', 'prev_wrong', 'stim_onset_per_trial'])

neural_matmat, predall_matmat = \
    utils.get_multiple_struct_fields_mat4py(raw_encoding_struct, ['neural_matmat', 'predall_matmat'])

# Package a single neuron into an object
def make_neuron_obj(rawdata, field, cellid):
    '''
    Create a list of neurons based on raw data
    :param rawdata: raw data returned by mat4py.loadmat
    :param field: field containing neural_act_mat
    :return: a neuron object
    '''
    neuron_activity = []
    ntrials = len(rawdata[field])
    for i in range(ntrials):
        trial_activity = utils.get_struct_field_mat4py(rawdata, 'neural_act_mat', True, i)
        neuron_activity.append(trial_activity[:, cellid])
        neuron = data_classes.Neuron(cellid, neuron_activity)
    return neuron

def make_neuron_list(rawdata, field, cell_lst):
    '''
    Create a list of neurons
    :param rawdata: raw data returned by mat4py.loadmat
    :param field: field containing neural_act_mat
    :param cell_lst: lst of ints, cell ids
    :return: a list of neuron objects
    '''
    neurons = []
    for cellid in cell_lst:
        print('Making neuron # ', cellid, '...')
        neuron = make_neuron_obj(rawdata, field, cellid)
        neurons.append(neuron)
    return neurons


neuron = make_neuron_obj(raw_encoding_struct, 'neural_act_mat', 0)
neuron.plot_all_trials()
#neurons = make_neuron_list(raw_encoding_struct, 'neural_act_mat', np.arange(180))


# For collecting all the neurons
def collect_neuron_aligned(raw_data, field, onset_arr, indices, window, neur_id):
    '''
    :param raw_data: neural activity matrix, ntrials cells, each is T x n_neurons
    :param field: name of field containing neural activity matrix
    :param onset_arr: a list, indices of t = 0 time points
    :param indices: a list indicating trial indices to extract from
    :param window: a tuple (start, end) window to extract
    :param neur_id: neuron to extract from
    :return: array of ntrials x T corresponding to extracted trials of the
    '''
    arr = []
    for i in range(len(indices)):
        trial_ctr = indices[i]
        startT = onset_arr[trial_ctr]
        trial_act = 1
        #arr.append(act_mat[trial_ctr][])


correct = utils.get_struct_field(raw_behavior_summary, 'correct')

