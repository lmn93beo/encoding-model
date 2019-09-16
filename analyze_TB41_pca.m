options.f_folder_name = 'C:\Users\Sur lab\Dropbox (MIT)\Sur\ExternalCode\encoding-model';
options.b_file_name = 'C:\Users\Sur lab\Dropbox (MIT)\Sur\ExternalCode\encoding-model\behavior_file_TB41.mat';

options.neuropil = 1;
options.neuropil_subt = 1;
options.dt = [-3 5];

[trials_dff, trials_z_dff, dff, z_dff, ix, ixCue] = getTrials_tb(options);
ncells = size(z_dff, 1);

%% Visualize
imagesc(dff);
%plot(dff(1,:));
hold on;
for i = 1:numel(ix)
   plot([ix(i) ix(i)], [-100 700], 'b');
   plot([ixCue(i) ixCue(i)], [-100 700], 'r');
end

%% Do a PCA on the neurons
% Find correlation matrix
demean = dff - mean(dff, 2); 
corr = demean * demean'; 
% SVD
[U, S, V] = svd(corr);

varexp = cumsum(diag(S));
varexp = varexp ./ sum(diag(S));
%stem(varexp)

% Project
pcs = U(:,1:3);
pcs_proj = demean' * pcs;

y1 = highpass(pcs_proj(:,1), 0.04, 5);
y2 = highpass(pcs_proj(:,2), 0.04, 5);
y3 = highpass(pcs_proj(:,3), 0.04, 5);

% Visualize first 3 pcs
figure;
plot(y1);
hold on;
plot(y2);
plot(y3);

for i = 1:numel(ix)
   plot([ix(i) ix(i)], [-700 700], 'b');
   plot([ixCue(i) ixCue(i)], [-700 700], 'r');
end


%% Load and preprocess behavior file
load behavior_file_TB41.mat
%1 = r, 2 = l stim location

ntrials = sum(ix < 9000);
choice = data.response.choice(1:ntrials);
%ntrials = numel(choice);

opp_contrasts = data.stimuli.opp_contrast(1:ntrials);
loc = data.stimuli.loc(1:ntrials);
contrast = data.params.contrast;

% Create the stim vectors                                                  
lstim = ones(1, ntrials) * contrast;
rstim = ones(1, ntrials) * contrast;

% TODO: check if 1 means L or R
lstim(loc == 1) = opp_contrasts(loc == 1);
rstim(loc == 2) = opp_contrasts(loc == 2);
difficulty = abs(lstim - rstim);

modified_RT = data.response.reward * 0;
modified_RT(data.response.reward > 0) = data.response.rewardtimes;
dt_rew = modified_RT - data.response.earlyCueTime;
dt_rew(~data.response.reward) = 0;
dt_rew_frames = floor(dt_rew / 0.2);

% For ball movement
% samples_start = cell2mat(data.response.samples_start);
%dX = data.response.dx;
ballT = data.response.time;
%ballT = ballT - ballT(1); % So that ballT starts at 0
fs = (length(ballT) - 1) / (ballT(end) - ballT(1));
dx_filt = lowpass(data.response.dx, 0.2);

% Convert cells to mat
precueS = cell2mat(data.response.precue_samples_start);
startS = cell2mat(data.response.samples_start);
rewardS = cell2mat(data.response.samples_reward);

% plot
ball_precueT = data.response.time(precueS);
ball_startT = data.response.time(startS);
ball_rewardT = data.response.time(rewardS);

%% Align on trialStart
pc_trials = {};
figure;
hold on;
nplot = 10;
ixValid = ix(ix < 9000-50);
pcs_filt = [y1 y2 y3];
ctrL = 1;
ctrR = 1;
ctr_rew = 1;
ctr_norew = 1;
ctr1 = 1;
ctr2 = 1;
ctr3 = 1;
ctr4 = 1;
ctr5 = 1;

pc_stim1 = {};
pc_stim2 = {};
pc_stim3 = {};
pc_stim4 = {};
pc_stim5 = {};

for i = 1:numel(ixValid)
    pc_trials{i} = pcs_filt(ix(i) - 10 : ix(i) + 10, :);
    if rstim(i) == 0
        pc_stim1{ctr1} = pcs_filt(ix(i) - 10 : ix(i) + 10, :);
        ctr1 = ctr1 + 1;
        %plot3(y1(ix(i)-20 : ix(i) + 30), y2(ix(i)-20 : ix(i) + 30), y3(ix(i)-20 : ix(i) + 30), 'b');
    elseif rstim(i) == 0.04
        pc_stim2{ctr2} = pcs_filt(ix(i) - 10 : ix(i) + 10, :);
        ctr2 = ctr2 + 1;
        %plot3(y1(ix(i)-20 : ix(i) + 30), y2(ix(i)-20 : ix(i) + 30), y3(ix(i)-20 : ix(i) + 30), 'r');
    elseif rstim(i) == 0.08
        pc_stim3{ctr3} = pcs_filt(ix(i) - 10 : ix(i) + 10, :);
        ctr3 = ctr3 + 1;
        
    elseif rstim(i) == 0.32
        pc_stim4{ctr4} = pcs_filt(ix(i) - 10 : ix(i) + 10, :);
        ctr4 = ctr4 + 1;
        
    elseif rstim(i) == 0.64
        pc_stim5{ctr5} = pcs_filt(ix(i) - 10 : ix(i) + 10, :);
        ctr5 = ctr5 + 1;
        
        
        %plot3(y1(ix(i)-20 : ix(i) + 30), y2(ix(i)-20 : ix(i) + 30), y3(ix(i)-20 : ix(i) + 30), 'k');
    end
    
    if data.response.reward(i) > 0
        pc_reward{ctr_rew} = pcs_filt(ix(i) - 10 : ix(i) + 10, :);
        ctr_rew = ctr_rew + 1;
    else
        pc_noreward{ctr_norew} = pcs_filt(ix(i) - 10 : ix(i) + 10, :);
        ctr_norew = ctr_norew + 1;
    end
end
plot3(y1(ix(1:nplot)), y2(ix(1:nplot)), y3(ix(1:nplot)), 'ro')

%% Plot
figure;
[pc1_arrR, pc2_arrR, pc3_arrR] = make_combined_arr(pc_trialsR);
[pc1_arrL, pc2_arrL, pc3_arrL] = make_combined_arr(pc_trialsL);
[pc1_reward, pc2_reward, pc3_reward] = make_combined_arr(pc_reward);
[pc1_noreward, pc2_noreward, pc3_noreward] = make_combined_arr(pc_noreward);
[pc1_arr1, pc2_arr1, pc3_arr1] = make_combined_arr(pc_stim1);
[pc1_arr2, pc2_arr2, pc3_arr2] = make_combined_arr(pc_stim2);
[pc1_arr3, pc2_arr3, pc3_arr3] = make_combined_arr(pc_stim3);
[pc1_arr4, pc2_arr4, pc3_arr4] = make_combined_arr(pc_stim4);
[pc1_arr5, pc2_arr5, pc3_arr5] = make_combined_arr(pc_stim5);


subplot(3,2,1);
taxis = ((1:21) - 10) / 5;
plot(taxis, pc1_arrR, 'Color', [0 0 1 0.2]);
hold on;
plot(taxis, pc1_arrL, 'Color', [1 0 0 0.2]);
%plot(taxis, mean(pc1_arrR'), 'r');
vline(0, 'k--');
title('PC1, right choices aligned on trial onset');
ylim([-600 600])

subplot(3,2,2)
errorbar(taxis, mean(pc1_arr1'), std(pc1_arr1')/sqrt(size(pc1_arr1, 2)), 'b')
hold on;
errorbar(taxis, mean(pc1_arr2'), std(pc1_arr2')/sqrt(size(pc1_arr2, 2)), 'r')
errorbar(taxis, mean(pc1_arr3'), std(pc1_arr3')/sqrt(size(pc1_arr3, 2)), 'k')
errorbar(taxis, mean(pc1_arr4'), std(pc1_arr4')/sqrt(size(pc1_arr4, 2)), 'g')
errorbar(taxis, mean(pc1_arr5'), std(pc1_arr5')/sqrt(size(pc1_arr5, 2)), 'b')
legend({'Right choice', 'Left choice'});
ylim([-600 600])
vline(0, 'k--');

subplot(3,2,3);
plot(taxis, pc2_arrR, 'Color', [0 0 1 0.2]);
hold on;
plot(taxis, pc2_arrL, 'Color', [1 0 0 0.2]);
%plot(taxis, mean(pc1_arrR'), 'r');
title('PC2');
ylim([-600 600])
vline(0, 'k--');

subplot(3,2,4)
errorbar(taxis, mean(pc2_arrR'), std(pc2_arrR')/sqrt(size(pc2_arrR, 2)), 'b')
hold on;
errorbar(taxis, mean(pc2_arrL'), std(pc2_arrL')/sqrt(size(pc2_arrL, 2)), 'r')
ylim([-600 600])
vline(0, 'k--');

subplot(3,2,5);
plot(taxis, pc3_arrR, 'Color', [0 0 1 0.2]);
hold on;
plot(taxis, pc3_arrL, 'Color', [1 0 0 0.2]);
%plot(taxis, mean(pc1_arrR'), 'r');
title('PC3');
ylim([-600 600])
vline(0, 'k--');

subplot(3,2,6)
errorbar(taxis, mean(pc3_arrR'), std(pc3_arrR')/sqrt(size(pc3_arrR, 2)), 'b')
hold on;
errorbar(taxis, mean(pc3_arrL'), std(pc3_arrL')/sqrt(size(pc3_arrL, 2)), 'r')
ylim([-600 600])
vline(0, 'k--');

%set(h2b.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
%%
% Subsample the trials with opp_contrast == 0.32 and choice ~= 5
goodtrials = (choice(2:end) ~= 5) & (choice(1:end-1) ~= 5);
goodtrials = logical([0 goodtrials]);
goodtrialsID = find(goodtrials);
%goodtrials = (choice ~= 5);
%goodtrials(1) = 0;
%goodtrials(258) = 0; %exclude this trial because neural activity has nan's
%TODO: consider not excluding this

% lstim = lstim(goodtrials);
% rstim = rstim(goodtrials);
% loc = loc(goodtrials);
% choice = choice(goodtrials);
% opp_contrasts = opp_contrasts(goodtrials);

% prevchoice = data.response.choice(goodtrialsID - 1);


% samples_start = samples_start(goodtrials);



%% Make the predictors
disp('Making predictor matrix...');
ntrialsgood = sum(goodtrials);
trial_lengths = diff(ixCue);

% Kernel for cue onset
for i = 1:ntrialsgood
    T = trial_lengths(goodtrialsID(i));
    cue_onset = zeros(T, 1);
    cue_onset(1) = 1;
    cue_onsetCells{i, 1} = cue_onset;
end

% Kernel for L-stim and R-stim onset
for i = 1:ntrialsgood
    T = trial_lengths(goodtrialsID(i));
    left_onset = zeros(T, 1);
    right_onset = zeros(T, 1);
    if loc(goodtrialsID(i)) == 2 % L-stim
        left_onset(ix(goodtrialsID(i)) - ixCue(goodtrialsID(i)) + 1) = 1;
    elseif loc(goodtrialsID(i)) == 1 % R-stim
        right_onset(ix(goodtrialsID(i)) - ixCue(goodtrialsID(i)) + 1) = 1;
    else
        error('Invalid trial. Only l/r choices are acceptable');
    end
    
    left_onsetCells{i, 1} = left_onset;
    right_onsetCells{i, 1} = right_onset;
end

% Reward times
rewards = zeros(ntrialsgood, T);
for i = 1:ntrialsgood
    T = trial_lengths(goodtrialsID(i));
    reward_arr = zeros(T, 1);
    if dt_rew_frames(goodtrialsID(i)) > 0
        reward_arr(dt_rew_frames(goodtrialsID(i)) + 1) = 1;
    end
    rewardsCell{i, 1} = reward_arr; 
end

% %%    
% Build the time information
balldata = cell(ntrialsgood, 1);
speedCell = cell(ntrialsgood, 1);
for i = 1:ntrialsgood
    %figure;
    T = trial_lengths(goodtrialsID(i));
    startT = ball_precueT(goodtrialsID(i));
    endT = ball_precueT(goodtrialsID(i) + 1);
    tpoints = linspace(startT, endT, T);
    
    start_id = find(ballT > startT, 1, 'first');
    end_id = find(ballT > endT, 1, 'first');
    
    %clean up multiple samples
    ballT_trial = ballT(start_id : end_id);
    [uniqueT, ic, ~] = unique(ballT_trial);
    dx_trial = data.response.dx(start_id : end_id);
    dx_unique = dx_trial(ic);
    
    resamp = interp1(uniqueT, dx_unique, tpoints);
    resamp(1) = dx_unique(1);
    
%     plot(ballT_trial, dx_trial);
%     hold on;
%     plot(tpoints, resamp);
   
    
    
    % Find the points of dx_filt that are closest to tpoints
    %[~,ixdx] = min(abs(tpoints - ballT),[],1);
    %balldata{i} = dx_filt(ixdx(1:end-1));
    
    balldata{i} = resamp';
    speedCell{i} = abs(resamp)';
end   

%% Balldata visualization
% start = ixdx(1);
% final = ixdx(end);
% figure;
% plot(start:final, dx_filt(start:final));
% hold on;
% plot(ixdx, dx_filt(ixdx), 'o');

%%
figure;
plot(ballT, dx_filt);
hold on;
for i = 1:ntrials
    start_id = ball_startT(i);
    cue_id = ball_precueT(i);
    plot([cue_id, cue_id], [-20 20], 'k');
    if choice(i) == 1
        plot([start_id, start_id], [-20 20], 'r');
    elseif choice(i) == 2
        plot([start_id, start_id], [-20 20], 'b');
    else
        plot([start_id, start_id], [-20 20], 'g');
    end
    
end
plot([1 ballT(end)], [-15 -15], 'k--');
plot([1 ballT(end)], [15 15], 'k--');


%% Visualize events
figure;
hold on;
for i = 1:ntrialsgood
    T = numel(cue_onsetCells{i});
    plot((1:T) + 50 * i, cue_onsetCells{i} * 20, 'g');
    plot((1:T) + 50 * i, left_onsetCells{i} * 20, 'b');
    plot((1:T) + 50 * i, right_onsetCells{i} * 20, 'r');
    plot((1:T) + 50 * i, rewardsCell{i} * 20, 'k');
    plot((1:T) + 50 * i, speedCell{i}, 'k--');
end

ylim([-30 30])

%% Make whole-trial variables
choiceGood = choice(goodtrialsID);
prevchoice = choice(goodtrialsID - 1);
difficultyGood = difficulty(goodtrialsID);
 

%% Make the predictor matrix
pred_types_cell = {'event', 'event', 'event', 'event', 'whole-trial', 'whole-trial', 'whole-trial', 'continuous'};
groupings = {1, 2, 3, 4, 5, 6, 7, 8};

[pred_allmat,pred_inds_cell,grouped_var_types] = make_predictor_matrix_generalcase({cue_onsetCells, left_onsetCells, right_onsetCells,...
     rewardsCell, choiceGood, prevchoice, difficultyGood, balldata}, ...
    pred_types_cell, groupings);

%% Make the neural matrix (response)
% Make the activity matrix
goodtrialsID = find(goodtrials);
neural_act_mat = cell(ntrialsgood, 1);
for i = 1:ntrialsgood
    T = trial_lengths(goodtrialsID(i));
    neural_act_mat{i} = z_dff(:, ixCue(goodtrialsID(i)) : ixCue(goodtrialsID(i) + 1) - 1)';
end
%%

neural_matmat = cell2mat(neural_act_mat);
predall_matmat = cell2mat(pred_allmat);

figure;
cellID = 1;
plot(neural_matmat(:,cellID));
hold on;
for i = 1:180
    plot([i * 41, i * 41], [-2 10], '--');
end

%% Do encoding model
approach = 'norefit';

%pred_types_cell_group = {'event', 'whole-trial', 'whole-trial', 'whole-trial', 'continuous'};
disp('Fitting encoding model...')
[relative_contrib,~,r2val] = process_encoding_model(pred_allmat, pred_inds_cell, neural_act_mat, pred_types_cell,approach);

%% Visualize
means = nanmean(relative_contrib, 1);
stds = nanstd(relative_contrib, 1) / sqrt(size(relative_contrib, 1));
figure;
bar(means);
hold on;
errorbar(1:size(relative_contrib, 2), means, stds, '.', 'Color', 'b');
%hold on;
xticks(1:size(relative_contrib, 2));

xticks(1:size(relative_contrib, 2));
factors = {'Cue', 'LeftStim', 'RightStim', 'Reward', 'Choice', 'PrevChoice', 'Difficulty', 'Movement'};
ylabel('Relative Contribution');

xticklabels(factors);
set(gca, 'FontSize', 16);


function [pc1arr, pc2arr, pc3arr] = make_combined_arr(raw_arr)

for i = 1:length(raw_arr)
    pc1arr(:,i) = raw_arr{i}(:,1);
    pc2arr(:,i) = raw_arr{i}(:,2);
    pc3arr(:,i) = raw_arr{i}(:,3);
end

end

