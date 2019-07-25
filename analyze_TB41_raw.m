options.f_folder_name = '/Users/tanyayang/Desktop/urop (summer19)/code/encoding-model';
options.b_file_name = '/Users/tanyayang/Desktop/urop (summer19)/code/encoding-model/behavior_file_TB41.mat';

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

% csvwrite('cue_onset.csv', cue_onsetCells)
% csvwrite('left_onset.csv', left_onsetCells)
% csvwrite('right_onset.csv', right_onsetCells)
% csvwrite('rewards.csv', rewardsCell)
% csvwrite('choices.csv', choiceGood)
% csvwrite('previous_choices.csv', prevchoice)
% csvwrite('difficulty.csv', difficultyGood)
% csvwrite('movement.csv', balldata)

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
% csvwrite('neural_activity.csv', neural_act_mat);
%parameter breakdown
%process_encoding_model is our function

% arguments: pred_allmat      - (our x value!) cell array correspoding to a matrix of behavioral predictors, each term is a trial and contains a matrix where rows are timepoints and columns are behavioral predictors.
%            pred_inds_cell   - cell array where each term has a vector of indices of the predictors that belong to a specific behavioral variable
%            neural_act_mat   - (our y value!) cell array correspoding to a matrix of activity traces, each term is a trial and contains a matrix where rows are timepoints and columns are traces correspoding to different
%                               neurons. Timepoints where the neuronal activity is not defined (e.g. a the imaging became unstable) are filled with NaNs (Currently it is assumed that before the first NaN
%                               activity was always defined, and after the first NaN activity is never defined).
%            pred_types_cell  - cell array where each terms indicates the type of behavioral variable ('event', 'whole-trial', or 'continuous').
%            approach         - 'norefit': calculate regression weights with the full model, then zero the weights correspoding to the predictors being dropped. 'refit' : calculate regression weights
%                               without the weights correspoding to the predictors being dropped (partial model).
%
% outputs:   relative_contrib - a matrix where rows are behavioral variables and columns are neurons. each term is the relative contribution of the behavioral variable to the neural activity.
%            Fstat_mat        - a matrix where rows are behavioral variables and columns are neurons. each term is the F-statistic associated with the nested model comparison where the predicted variable is the
%                               activity of the neuron, the full model is the model containing the predictors from al behavioral variables, and the partial model is the model containing the predictors from all
%                               variables except the one being tested. The value of this statistic shoudl be compared to a distirbution of statistic values obtained by erforing the same oprateion on shuffled 
%                               data, directly using the p-value assocaited with the staitstic is not valid given the autocorrelations in the data
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

%csvwrite('relative_contributions_raw.csv',relative_contrib)

