import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import pandas as pd
import utils
import mat4py

# Load data
raw_behavior_summary = mat4py.loadmat('TB41_behavior_summary.mat')
raw_encoding_struct = mat4py.loadmat('TB41_encoding_structs.mat')

# Extract structures from the raw data
one, two, three, four, correct, incorrect, left, right, prev_right, prev_wrong, \
    stim_onset_per_trial = utils.get_multiple_struct_fields(raw_behavior_summary,
                        ['one', 'two', 'three', 'four', 'correct', 'incorrect',
                         'left', 'right', 'prev_right', 'prev_wrong', 'stim_onset_per_trial'])

left_onsetCells, neural_act_mat, pred_allmat, rewardsCell, right_onsetCells = \
    utils.get_multiple_struct_fields(raw_encoding_struct, ['left_onsetCells', 'neural_act_mat', 'pred_allmat',
                                                           'rewardsCell', 'right_onsetCells'])

# For collecting all the neurons
def collect_neuron_aligned(act_mat, onset_arr, indices, window, neur_id):
    '''
    :param act_mat: neural activity matrix, ntrials cells, each is T x n_neurons
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
        #arr.append(act_mat[trial_ctr][])


correct = utils.get_struct_field(raw_behavior_summary, 'correct')

