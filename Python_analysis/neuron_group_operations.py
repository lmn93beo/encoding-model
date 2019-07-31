import utils
import data_classes
import numpy as np
import matplotlib.pyplot as plt

class NeuronGroup(object):
    '''
    A class for a group of Neuron objects
    '''
    def __init__(self, neuron_lst):
        self.neurons = neuron_lst
        self.n_neurons = len(neuron_lst)
        self.classes = []
        self.class0cells = []
        self.class1cells = []
        self.class2cells = []

    def align_all(self, tpoints, window):
        '''
        Align all neurons in group
        :param tpoints: time points to align to (one-indexed)
        :param window: a tuple (start, end), range of time points for cropping
        :return: nothing
        '''
        for neuron in self.neurons:
            neuron.align_activity(tpoints, window)

    def classify_all(self):
        '''
        Classify all neurons in group
        Keeping a list of neuron classes
        :return: an np array giving the classes for all neurons
        '''
        classes = []
        for neuron in self.neurons:
            neuron.classify()
            classes.append(neuron.neuron_class)
        self.classes = np.array(classes)
        self.class0cells = np.where(self.classes == 0)[0]
        self.class1cells = np.where(self.classes == 1)[0]
        self.class2cells = np.where(self.classes == 2)[0]
        return np.array(classes)



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

def compare_conditions(rawdata, field, tpoints, window, cellid, condA, condB):
    '''
    Plot graphs to compare between condition A and B
    :param cellid: cell number
    :param condA: np array to indicate trials belonging to condition A
    :param condB: np array to indicate trials belonging to condition B
    :return:
    '''
    neuron = make_neuron_obj(rawdata, field, cellid)
    neuron.align_activity(tpoints, window)
    plt.subplot(121)
    neuron.plot_all_trials(condA)
    plt.ylim((-2, 2))
    plt.subplot(122)
    neuron.plot_all_trials(condB)
    plt.ylim((-2, 2))


