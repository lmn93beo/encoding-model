% Here, we try to plot signals for consecutive trials (i.e. L-L, L-R, R-L,
% R-R etc), and see if we can 'decode' the next action from what happened
% in the previous trial (aka, if the previous action was L, is the animal
% more likely to choose L?)

options.f_folder_name = 'C:\Users\Sur lab\Dropbox (MIT)\Sur\ExternalCode\encoding-model';
options.b_file_name = 'C:\Users\Sur lab\Dropbox (MIT)\Sur\ExternalCode\encoding-model/behavior_file_TB41.mat';

options.neuropil = 1;
options.neuropil_subt = 1;
options.dt = [-3 5];

[trials_dff, trials_z_dff, dff, z_dff, ~, ~, ~] = getTrials_tb(options);
ncells = size(z_dff, 1);

load TB41_epochs.mat

class0 = [2,   3,   5,   8,  10,  11,  16,  19,  20,  21,  22,  23,  24,...
        26,  27,  28,  29,  31,  32,  33,  36,  37,  38,  40,  45,  48,...
        50,  51,  53,  54,  55,  57,  58,  63,  71,  75,  76,  78,  79,...
        80,  82,  85,  87,  95,  96,  99, 100, 104, 105, 112, 122, 124,...
       128, 129, 130, 134, 135, 137, 139, 143, 145, 147, 150, 151, 158,...
       162, 163, 164, 165, 167, 168, 169, 170, 171, 172, 174, 175, 176,...
       179];
   
class1 = [ 1,   4,   6,  12,  13,  35,  39,  41,  43,  46,  49,  52,  56,...
        59,  65,  67,  68,  69,  73,  74,  81,  89,  90,  91,  93,  94,...
       101, 102, 103, 106, 107, 108, 113, 114, 116, 117, 125, 126, 127,...
       131, 132, 138, 148, 153, 154, 155, 156, 159, 160, 173, 177, 180];
   
class2 = [ 7,   9,  14,  15,  17,  18,  25,  30,  34,  42,  44,  47,  60,...
        61,  62,  64,  66,  70,  72,  77,  83,  84,  86,  88,  92,  97,...
        98, 109, 110, 111, 115, 118, 119, 120, 121, 123, 133, 136, 140,...
       141, 142, 144, 146, 149, 152, 157, 161, 166, 178];


%% Load and preprocess behavior file
load behavior_file_TB41.mat

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

% Get the ix of the appropriate trials
Accuracies = [];
trange = -30:60;
for i = 1:100
    disp(i);
    accuracies = [];
    dff0 = dff(class0,:);
    for t = trange
        acc = get_decoding_acc(dff0, ixOutcome, ltrials, rtrials, t, 0);
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
        acc = get_decoding_acc(dff1, ixOutcome, ltrials, rtrials, t, 0);
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
        acc = get_decoding_acc(dff, ixOutcome, ltrials, rtrials, t, 0);
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
legend([l0, l1, l2], {'Group 1', 'Group2', 'All'});
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



