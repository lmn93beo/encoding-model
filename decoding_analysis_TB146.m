% Here, we try to plot signals for consecutive trials (i.e. L-L, L-R, R-L,
% R-R etc), and see if we can 'decode' the next action from what happened
% in the previous trial (aka, if the previous action was L, is the animal
% more likely to choose L?)

options.f_folder_name = 'C:\Users\Sur lab\Dropbox (MIT)\Sur\Physiology_analysis\146\20190807\suite2p\plane0';
options.b_file_name = 'C:\Users\Sur lab\Dropbox (MIT)\Sur\ExternalCode\encoding-model/20190807_trackball_0146.mat';
root_filename = 'TB146';
encoding_struct_fname = [root_filename '_encoding_structs'];
behavior_summary_fname = [root_filename '_behavior_summary'];
savefile = 1;
options.f_file_name = 'F_wNeuropil_partial_tb41_03032017.txt';

options.neuropil = 1;
options.neuropil_subt = 1;
options.dt = [-3 5];
options.suite2p = 1;

[trials_dff, trials_z_dff, dff, z_dff, frametimes, ix, ixCue] = getTrials_tb(options);
ncells = size(z_dff, 1);

load TB146_epochs.mat

class0 = [ 1,  5, 10, 13, 18, 25, 26, 30, 43, 65, 70, 71];
   
class1 = [ 2,  3,  7,  8,  9, 14, 15, 16, 20, 21, 22, 24, 28, 33, 35, 39, 47,...
       49, 59, 64, 69, 72, 73, 74, 75];
   
class2 = [4,  6, 11, 12, 17, 19, 23, 27, 29, 31, 32, 34, 36, 37, 38, 40, 41,...
       42, 44, 45, 46, 48, 50, 51, 52, 53, 54, 55, 56, 57, 58, 60, 61, 62,...
       63, 66, 67, 68, 76];
   
   
   



%% Load and preprocess behavior file
load(options.b_file_name)

%1 = 1, r = 2 stim location
ntrials = sum(ixStart < 9000);
choice = data.response.choice(1:ntrials);

opp_contrasts = data.stimuli.opp_contrast(1:ntrials);
loc = data.stimuli.loc(1:ntrials);
contrast = data.params.contrast;

% Create the stim vectors                                                  
lstim = ones(1, ntrials) * contrast;
rstim = ones(1, ntrials) * contrast;

% TODO: check if 1 means L or R
rstim(loc == 1) = opp_contrasts(loc == 1);
lstim(loc == 2) = opp_contrasts(loc == 2);
difficulty = abs(lstim - rstim);

modified_RT = data.response.reward * 0;
modified_RT(data.response.reward > 0) = data.response.rewardtimes;
dt_rew = modified_RT - data.response.earlyCueTime;
dt_rew(~data.response.reward) = nan;
dt_rew_frames = floor(dt_rew / 0.2);

%% Break into trials
prevchoice = [nan choice(1:end-1)];
ntrials = sum(ixOutcome < 9000);
rewards = data.response.reward(1:ntrials);
ltrials = find(choice == 1);
ltrials = ltrials(1:end-1);
rtrials = find(choice == 2);

n = min(numel(ltrials), numel(rtrials));
fprintf('Number of trials in each condition is %d\n', n);
ltrials = randsample(ltrials, n);
rtrials = randsample(rtrials, n);




% Get the ix of the appropriate trials
Accuracies = [];
trange = -50:100;
for i = 1:100
    disp(i);
    accuracies = [];
    dff0 = dff(class0,:);
    for t = trange
        acc = get_decoding_acc(dff0, ixStart, ltrials, rtrials, t, 0);
        accuracies(t - trange(1) + 1) = acc;
    end
%     plot(trange * 0.2, accuracies);
%     hold on;
%     vline(0, 'k--');
%     hline(0.5, 'k--')
    Accuracies(i,:) = accuracies;
end
%%
l0 = plot(trange, mean(Accuracies, 1));
hold on;
stderr = std(Accuracies, [], 1) / sqrt(size(Accuracies, 1));
plot(trange, mean(Accuracies, 1) + stderr, 'b--');
plot(trange, mean(Accuracies, 1) - stderr, 'b--');


%% Stim decoding
% ltrials = find(choice == 1);
% ltrials = ltrials(2:end-1);
% rtrials = find(choice == 2);

% Get the ix of the appropriate trials
Accuracies2 = [];
for i = 1:100
    disp(i);
    accuracies = [];
    dff1 = dff(class1,:);
    for t = trange
        acc = get_decoding_acc(dff1, ixStart, ltrials, rtrials, t, 0);
        accuracies(t - trange(1) + 1) = acc;
    end
%     plot(trange * 0.2, accuracies);
%     hold on;
%     vline(0, 'k--');
%     hline(0.5, 'k--')
    Accuracies2(i,:) = accuracies;
end
%%
l1 = plot(trange, mean(Accuracies2, 1), 'r');
hold on;
stderr_stim = std(Accuracies2, [], 1) / sqrt(size(Accuracies2, 1));
plot(trange, mean(Accuracies2, 1) + stderr_stim, 'r--');
plot(trange, mean(Accuracies2, 1) - stderr_stim, 'r--');

%% Stim decoding
% ltrials = find(choice == 1);
% ltrials = ltrials(2:end-1);
% rtrials = find(choice == 2);

% Get the ix of the appropriate trials
Accuracies3 = [];
for i = 1:100
    disp(i);
    accuracies = [];
    dff2 = dff(class2,:);
    for t = trange
        acc = get_decoding_acc(dff2, ixStart, ltrials, rtrials, t, 0);
        accuracies(t - trange(1) + 1) = acc;
    end
%     plot(trange * 0.2, accuracies);
%     hold on;
%     vline(0, 'k--');
%     hline(0.5, 'k--')
    Accuracies3(i,:) = accuracies;
end
%%
l2 = plot(trange, mean(Accuracies3, 1), 'k');
hold on;
stderr_stim = std(Accuracies3, [], 1) / sqrt(size(Accuracies3, 1));
plot(trange, mean(Accuracies3, 1) + stderr_stim, 'k--');
plot(trange, mean(Accuracies3, 1) - stderr_stim, 'k--');

%%
legend([l0, l1, l2], {'Group 1', 'Group2', 'Group3'});
title('Choice decoding');



function acc = get_decoding_acc(dff, ix, ltrials, rtrials, t, bootstrap)
% bootstrap: if True, perform sampling with replacement
    dff_left = dff(:, ix(ltrials) + t);
    dff_right = dff(:, ix(rtrials) + t);
    dff_all = [dff_left dff_right];
    
    labels = [ones(1, size(dff_left, 2)) ones(1, size(dff_right, 2)) * 2];
    
    if bootstrap
        order_boot = 1:numel(labels);
        order_boot = datasample(order_boot, numel(labels));
        dff_all = dff_all(:, order_boot);
        labels = labels(order_boot);
    end
    
    ntrials = size(dff_all, 2);
    order = randperm(ntrials);

    labels_p = labels(order);
    dff_p = dff_all(:, order);


    % Split into training and test set: 80% train, 20% test
    ntrain = floor(0.8 * ntrials);
    labels_train = labels_p(1:ntrain);
    labels_test = labels_p(ntrain + 1:end);
    dff_train = dff_p(:, 1:ntrain);
    dff_test = dff_p(:, ntrain + 1:end);


    %% Train
    % Average
    dff_train_L = dff_train(:, labels_train == 1);
    dff_train_R = dff_train(:, labels_train == 2);

    mean_L = mean(dff_train_L, 2);
    mean_R = mean(dff_train_R, 2);

    % Test
    corr_L = corr(dff_test, mean_L);
    corr_R = corr(dff_test, mean_R);
    prediction = (corr_R > corr_L) + 1;
    acc = sum(prediction == labels_test') / numel(labels_test);
    %fprintf('Accuracy: %.2f\n', acc);

end



