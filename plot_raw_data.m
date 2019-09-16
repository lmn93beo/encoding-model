load TB41_behavior_summary.mat
load TB41_encoding_structs.mat

%% Plot first neuron
neur_ctr = 5;
collector = [];
n_neurons = size(neural_matmat, 2);

prev_right_acts = neural_act_mat(prev_right);
prev_wrong_acts = neural_act_mat(prev_wrong);

prev_right_arr = [];
prev_wrong_arr = [];

for i = 1:numel(prev_right)
    trial_ctr = prev_right(i);
    startT = stim_onset_per_trial(trial_ctr);
    prev_right_arr(i,:) = prev_right_acts{i}(startT - 5 : startT + 10,neur_ctr);
end

for i = 1:numel(prev_wrong)
    trial_ctr = prev_wrong(i);
    startT = stim_onset_per_trial(trial_ctr);
    prev_wrong_arr(i,:) = prev_wrong_acts{i}(startT - 5 : startT + 10,neur_ctr);
end

%plot(prev_right_arr', 'b')

%%
prev_right_arr2 = collect_neuron_aligned(neural_act_mat, stim_onset_per_trial,...
    1:173, [-5 18], 1);

[meansL, stdL] = collect_multiple_neurons(neural_act_mat, stim_onset_per_trial,...
    1:173, [-5, 18], 1:n_neurons);

[meansR, stdR] = collect_multiple_neurons(neural_act_mat, stim_onset_per_trial,...
    right, [-5, 18], 1:n_neurons);

%% Sort by peak activity
pk_locsL = nan(1, n_neurons);
for i = 1:n_neurons
    [~, loc] = max(meansL(i,:));
    pk_locsL(i) = loc;
    
    [~, loc] = max(meansR(i,:));
    pk_locsR(i) = loc;
end

% Sort
[pk_locs_sortedL, ids] = sort(pk_locsL);
mean_sortedL = meansL(ids, :);
mean_sortedR = meansR(ids, :);

color_axis = [-1 1];

figure;
subplot(1,2,1);
imagesc(mean_sortedL);
rdbl = redblue(100);
colormap(rdbl);
caxis(color_axis)
title('Left stim')
colorbar;
ylabel('Neuron');
xlabel('Time');

subplot(1,2,2);
imagesc(mean_sortedR);
colormap(rdbl);
caxis(color_axis)
colorbar
title('Right stim')
xlabel('Time')

%% Effect of trial difficulty
left_one = intersect(left, one);
left_two = intersect(left, two);
left_three = intersect(left, three);
left_four = intersect(left, four);

blues = [198,219,239;
    107,174,214;
    33,113,181;
    8,48,107]/255;

for neur_ctr = 1:n_neurons

    [meansL_one, stdL_one] = collect_multiple_neurons(neural_act_mat, stim_onset_per_trial,...
        left_one, [-5, 17], neur_ctr);
    [meansL_two, stdL_two] = collect_multiple_neurons(neural_act_mat, stim_onset_per_trial,...
        left_two, [-5, 17], neur_ctr);
    [meansL_three, stdL_three] = collect_multiple_neurons(neural_act_mat, stim_onset_per_trial,...
        left_three, [-5, 17], neur_ctr);
    [meansL_four, stdL_four] = collect_multiple_neurons(neural_act_mat, stim_onset_per_trial,...
        left_four, [-5, 17], neur_ctr);
    
    figure;
    hold on;
    errorbar(meansL_one, stdL_one / sqrt(numel(left_one)), 'Color', blues(1,:));
    errorbar(meansL_two, stdL_two / sqrt(numel(left_two)), 'Color', blues(2,:));
    errorbar(meansL_three, stdL_three / sqrt(numel(left_three)), 'Color', blues(3,:));
    errorbar(meansL_four, stdL_four / sqrt(numel(left_four)), 'Color', blues(4,:));
    legend({'one', 'two', 'three', 'four'});
end





%%
function arr = collect_neuron_aligned(act_mat, onset_arr, indices, window, neur_ctr)
% Inputs:
% act_mat: neural activity matrix, ntrials cells, each is T x n_neurons
% onset_arr: indices of t = 0 time points
% indices: trial indices to extract from
% neur_ctr: neuron to extract from
% Returns: array of ntrials x T corresponding to extracted trials of the
% neuron

arr = [];
for i = 1:numel(indices)
    trial_ctr = indices(i);
    startT = onset_arr(trial_ctr);
    arr(i,:) = act_mat{trial_ctr}(startT + window(1) : startT + window(2), neur_ctr);
end
end

function [means, stds] = collect_multiple_neurons(act_mat, onset_arr, indices, window,...
    neuron_lst)
% Inputs: same as collect_neuron_aligned
% Returns: array of n_neurons x T corresponding to trial-averaged activity
% of all the neurons specified

for i = 1:numel(neuron_lst)
    neuron_arr = collect_neuron_aligned(act_mat, onset_arr, indices, ...
        window, neuron_lst(i));
    means(i,:) = mean(neuron_arr, 1);
    stds(i,:) = std(neuron_arr, [], 1);
end
end






