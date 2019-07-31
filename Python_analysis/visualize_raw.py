import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import pandas as pd
import utils
import mat4py
import neuron_group_operations
import data_classes

# Load data
raw_behavior_summary = mat4py.loadmat('TB41_behavior_summary.mat')
raw_encoding_struct = mat4py.loadmat('TB41_encoding_structs.mat')

# Extract structures from the raw data
one, two, three, four, correct, incorrect, left, right, prev_right, prev_wrong, \
    stim_onset_per_trial = utils.get_multiple_struct_fields_mat4py(raw_behavior_summary,
                        ['one', 'two', 'three', 'four', 'correct', 'incorrect',
                         'left', 'right', 'prev_right', 'prev_wrong', 'stim_onset_per_trial'])

# Being careful about one-indexing
one -= 1
two -= 1
three -= 1
four -= 1
correct -= 1
incorrect -= 1
prev_right -= 1
prev_wrong -= 1
left -= 1
right -= 1

neural_matmat, predall_matmat = \
    utils.get_multiple_struct_fields_mat4py(raw_encoding_struct, ['neural_matmat', 'predall_matmat'])

# Compare reward and non-reward for cells in class 2
to_plot = [ 6,   8,  13,  14,  16,  17,  24,  29,  33,  41,  43,  46,  59,
        60,  61,  63,  65,  69,  71,  76,  82,  83,  85,  87,  91,  96,
        97, 108, 109, 110, 114, 117, 118, 119, 120, 122, 132, 135, 139,
       140, 141, 143, 145, 148, 151, 156, 160, 165, 177]
'''
for i, neuron_id in enumerate(to_plot):
    plt.figure()
    neuron_group_operations.compare_conditions(raw_encoding_struct, 'neural_act_mat', stim_onset_per_trial, [-5, 18],
                                               neuron_id, correct, incorrect)
    plt.title('Neuron %d' % neuron_id)
'''

#neurons = neuron_group_operations.make_neuron_list(raw_encoding_struct, 'neural_act_mat', np.arange(180))
# Create a neuron group, align and classify
print('Aligning and classifying...')
neuron_group = neuron_group_operations.NeuronGroup(neurons)
neuron_group.align_all(stim_onset_per_trial, [-5, 18])
classes = neuron_group.classify_all()
print('%d neurons in class 0, %d neurons in class 1, %d neurons in class 2 ' % (np.sum(classes == 0),
                                                                                np.sum(classes == 1), np.sum(classes == 2)))
