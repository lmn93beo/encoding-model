import numpy as np
import matplotlib.pyplot as plt

class Neuron(object):
    '''
    A class for a neuron in ACC
    '''
    def __init__(self, cellid, activity):
        '''
        :param cellid: neuron number
        :param activity: A list of ntrials np arrays, each array is 1-d activity of that trial
        '''
        self.id = cellid
        self.activity = activity
        self.ntrials = len(activity)
        self.aligned_activity = []

    def plot_all_trials(self):
        '''
        For plotting of all trials
        :return: nothing
        '''
        for i in range(self.ntrials):
            print('Plotting trial', i)
            plt.plot(self.activity[i], 'b', alpha=0.2)
        plt.show()

    def align_activity(self, tpoints, window):
        '''
        Align neural activity to time points specified in tpoints
        :param tpoints: time points to align to (one-indexed)
        :param window: a tuple (start, end), range of time points for cropping
        :return: an np array of size ntrials x T corresponding to aligned activity
        '''
        arr = []
        for i in range(self.ntrials):
            startT = tpoints[i]
            single_trial = self.activity[i][startT + window[0] - 1 : startT + window[1]]
            arr.append(single_trial)
        self.aligned_activity = np.array(arr)
        return np.array(arr)



class Experiment(object):
    '''
    A class of trackball experiment
    '''
    def __init__(self, rawdata):
        self.rawdata = rawdata
