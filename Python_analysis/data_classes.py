import numpy as np
import matplotlib.pyplot as plt

class Neuron(object):
    '''
    A class for a neuron in ACC
    '''
    def __init__(self, cellid, activity):
        '''
        :param cellid: neuron number
        :param activity: A list of np arrays, each array is 1-d activity of that trial
        '''
        self.id = cellid
        self.activity = activity
        self.ntrials = len(activity)

    def plot_all_trials(self):
        '''
        For plotting of all trials
        :return: nothing
        '''
        for i in range(self.ntrials):
            print('Plotting trial', i)
            plt.plot(self.activity[i], 'b', alpha=0.2)
        plt.show()



class Experiment(object):
    '''
    A class of trackball experiment
    '''
    def __init__(self, rawdata):
        self.rawdata = rawdata
