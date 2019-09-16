clear all
close all

options.f_folder_name = 'C:\Users\Sur lab\Dropbox (MIT)\Sur\Physiology_analysis\146\20190802\suite2p\plane0';
options.b_file_name = 'C:\Users\Sur lab\Dropbox (MIT)\Nhat\mouse 0146/20190802_trackball_0146.mat';
root_filename = 'TB146_20190802';
encoding_struct_fname = [root_filename '_encoding_structs'];
behavior_summary_fname = [root_filename '_behavior_summary'];
epoch_filename = [root_filename '_epochs'];
savefile = 1;

options.neuropil = 1;
options.neuropil_subt = 1;
options.dt = [-3 5];
options.suite2p = 1;

[trials_dff, trials_z_dff, dff, z_dff, frametime, ix, ixCue] = getTrials_tb(options);
ncells = size(z_dff, 1);

% High-pass filter
% fprintf('Doing high pass filter...\n');
% for i = 1:ncells
%     original = z_dff(i,:);
%     filtered = highpass(original, 0.1, 5);
%     env = abs(hilbert(filtered));
%     env = lowpass(env, 0.001, 5);
%     z_dff(i,:) = filtered ./ env;
% end

%% Visualize
imagesc(z_dff);
%plot(dff(1,:));
hold on;
for i = 1:numel(ix)
   plot([ix(i) ix(i)], [-100 700], 'b');
   plot([ixCue(i) ixCue(i)], [-100 700], 'r');
end


%% Load and preprocess behavior file
load(options.b_file_name);
%1 = 1, r = 2 stim location

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
ballT = data.response.time;
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
 
    balldata{i} = resamp';
    speedCell{i} = abs(resamp)';
end   

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
pred_types_cell = {'event', 'event', 'event', 'event', 'continuous'};
base_vars = {cue_onsetCells, left_onsetCells, right_onsetCells,...
     rewardsCell, balldata};
spline_filename = 'Splines/spline_basis_9_order5_30pts.mat';

[pred_allmat,pred_inds_cell,grouped_var_types] = make_predictor_matrix_generalcase(base_vars, pred_types_cell, spline_filename);

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

if savefile
    save(encoding_struct_fname, 'neural_matmat', 'predall_matmat', 'neural_act_mat',...
         'pred_allmat', 'left_onsetCells', 'right_onsetCells', 'rewardsCell');
end

%% Do encoding model
approach = 'norefit';

%pred_types_cell_group = {'event', 'whole-trial', 'whole-trial', 'whole-trial', 'continuous'};
disp('Fitting encoding model...')
[relative_contrib,~,r2val] = process_encoding_model(pred_allmat, pred_inds_cell, neural_act_mat, pred_types_cell,approach);
fprintf('Mean R^2 value = %.2f%% \n', mean(r2val) * 100);

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

%% reorganization with cue_cells
tfs = 15; %tfs = timeframes
cue_cells = cell(tfs, 1); %within each cell is an added double of 173 values
left = []; %only giving us the i (index) of trials that have left stimulus
right = []; %same idea with right
     
stim_onset_per_trial = [];

correct = [];
incorrect = [];
prev_right = [];
prev_wrong = [];

one = [];
two = [];
three = [];
four = [];

for i = 1:length(neural_act_mat) %looping through all the trials
    
    for j = 1:tfs
       cue_cells{j} = [cue_cells{j}, neural_act_mat{i}(j)] ; %adding the neural activity that corresponds to each cue onset
    end
    
    %finding out when the onset occurs
    if find(left_onsetCells{i})
        mid = find(left_onsetCells{i});
        left = [left, i];
    elseif find(right_onsetCells{i})
        mid = find(right_onsetCells{i});
        right = [right, i];
    end
    
    %tracking which trials are correct
    if any(rewardsCell{i}) %if a reward is presented
        correct = [correct, i];
        if ~(i + 1 == 174);
            prev_right = [prev_right, i + 1];
        end
    else %when there's no reward
        incorrect = [incorrect, i];
        if i~=173;
            prev_wrong = [prev_wrong, i + 1];
        end
    end
    
     %finding out the level of difficulty
    val = difficultyGood(i);
    %organizing the various levels of difficulty
    if val == 0.3200
        one = [one, i];
    elseif val == 0.5600
        two = [two, i];
    elseif val == 0.6000
        three = [three, i];
    elseif val == 0.6400
        four = [four, i];
    end
end

if savefile
    save(behavior_summary_fname, 'one', 'two', 'three', 'four', 'correct',...
        'incorrect', 'left', 'right', 'prev_right', 'prev_wrong');
end

%% Define important epochs
load(options.b_file_name);
earlyCue = data.response.earlyCueTime';
ts = data.response.trialstart';
rewardtime = nan(size(ts));
rewardtime(data.response.reward > 0) = data.response.rewardtimes';
endtime = data.response.trialtime';
responsetime = cellfun(@(x) x(end), data.response.timePC)';
outcometime = ts + responsetime;

% Convert to frame number
ixCue = find_ix_frames(frametime, earlyCue);
ixStart = find_ix_frames(frametime, ts);
ixReward = find_ix_frames(frametime, rewardtime);
ixReward(isnan(rewardtime)) = nan;
ixOutcome = find_ix_frames(frametime, outcometime);
ixEnd = find_ix_frames(frametime, endtime);

% Find good trials
ntrials = sum(ix < 9000);
choice = data.response.choice(1:ntrials);
prevchoice = [nan choice(1:end-1)];
goodtrials = (choice ~= 5) & (prevchoice ~= 5);
goodtrialsID = find(goodtrials);



%% Save
save(epoch_filename, 'ixCue', 'ixStart', 'ixReward', 'ixOutcome', 'ixEnd', 'goodtrials');



%% Helper functions
function ix = find_ix_frames(frametime, tlist)
% frametime: array of frame times
% tlist: array of time points
% Returns the index of frames specified in tlist
    dum = repmat(frametime,numel(tlist), 1);
    [~, ix] = min(abs(dum - tlist),[], 2);
end

