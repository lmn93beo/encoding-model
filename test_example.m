var1 = {[1 2 3 4 5]', [1 5 19 -2 3]', [5 9 3 2 1]'};
var2 = {[0 0 0 1 0], [0 0 0 1 0], [0 0 0 1 0]};
var3 = [1 -1 0];

[pred_allmat,pred_inds_cell,grouped_var_types] = make_predictor_matrix_generalcase({var1, var2, var3}, {'continuous', 'event', 'whole-trial'});

neuron = {[0 0 0 1 0; 0 0 0 1 0; 0 0 0 1 0], };

