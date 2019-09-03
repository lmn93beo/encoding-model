import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import pandas as pd
import utils
import mat4py
import neuron_group_utils
import neuron_utils


def test_load_data():
    '''
    Test that data is loaded correctly
    '''
    # Load data
    raw_behavior_summary = mat4py.loadmat('TB41_behavior_summary.mat')
    raw_encoding_struct = mat4py.loadmat('TB41_encoding_structs.mat')
    assert len(raw_behavior_summary.keys()) == 11
    assert len(raw_encoding_struct.keys()) == 7
    assert len(raw_encoding_struct['left_onsetCells']) == 173
    assert len(raw_encoding_struct['right_onsetCells']) == 173
    assert len(raw_encoding_struct['rewardsCell']) == 173
    assert len(raw_behavior_summary['correct']) == 153
    assert len(raw_behavior_summary['incorrect']) == 20
    assert len(raw_behavior_summary['left']) == 79



