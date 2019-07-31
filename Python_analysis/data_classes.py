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
        self.neuron_class = -1

    def plot_all_trials(self, trials=None):
        '''
        For plotting of all trials
        :param trials: list of trials, if none, will plot all the trials
        :return: nothing
        '''
        assert self.aligned_activity != []

        if trials is None:
            trials = np.arange(self.ntrials)

        extracted_trials = self.aligned_activity[trials, :]
        print(extracted_trials.shape)
        plt.plot(extracted_trials.T, 'b', alpha=0.1)
        plt.plot(np.mean(extracted_trials, axis=0), 'r')
        #for i in trials:
        #    print('Plotting trial', i)
        #    plt.plot(self.aligned_activity[], 'b', alpha=0.2)
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

    def classify(self):
        '''
        Classify the neuron based on the peak of trial-averaged activity
        Peak time: class 0 (0 - 9 frames), class 1 (10-15 frames), class 2 (>= 16 frames)
        :return: class of the neuron
        '''
        # Find peak of mean activity
        assert(self.aligned_activity != [])
        mean_activity = np.mean(self.aligned_activity, axis=0)
        t_max_activity = np.argmax(mean_activity)
        if t_max_activity < 10:
            self.neuron_class = 0
        elif t_max_activity < 16:
            self.neuron_class = 1
        else:
            self.neuron_class = 2
        return self.neuron_class


class Experiment(object):
    '''
    A class of trackball experiment
    '''
    def __init__(self, rawdata):
        self.rawdata = rawdata
