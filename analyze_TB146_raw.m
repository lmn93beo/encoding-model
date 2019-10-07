clear all
close all

choice = input('Warning! Will overwrite existing encoding structs. Press 1 to continue, 0 to cancel \n >>');
if choice ~= 1
    error('User cancelled');
end


date = '20190807';
options.f_folder_name = ['C:\Users\Sur lab\Dropbox (MIT)\Sur\Physiology_analysis\146\' date '\suite2p\plane0'];
options.b_file_name = ['C:\Users\Sur lab\Dropbox (MIT)\trackball-behavior\Data\146_all/' date '_trackball_0146.mat'];
root_filename = ['Python_analysis\TB146_' date];
encoding_struct_fname = [root_filename '_encoding_structs'];
behavior_summary_fname = [root_filename '_behavior_summary'];
epoch_filename = [root_filename '_epochs'];
savefile = 1;

if strcmp(date, '20190724')
    options.special = 1;
else
    options.special = 0;
end

options.neuropil = 1;
options.neuropil_subt = 1;
options.dt = [-2 3];
options.suite2p = 1;
options.align_by = 'trial_start';

[trials_dff, trials_z_dff, dff, z_dff, frametime, ix, ixCue] = getTrials_tb(options);
ncells = size(z_dff, 1);
save(['C:\Users\Sur lab\Dropbox (MIT)\Sur\evidence_analysis\TB146_' date '_trials_dff.mat'], 'trials_dff', 'trials_z_dff');

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

ntrials = sum(ix < max(ix)) - 1;
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

%% Make the neural matrix (response)
% Make the activity matrix
for i = 1:numel(choice)
    neural_act_mat{i} = dff(:, ixCue(i) : ixCue(i + 1) - 1)';
end

neural_act_mat = neural_act_mat';
neural_matmat = cell2mat(neural_act_mat);

if savefile
    save(encoding_struct_fname, 'neural_matmat', 'neural_act_mat');
end


%% reorganization with cue_cells
tfs = 15; %tfs = timeframes
cue_cells = cell(tfs, 1); %within each cell is an added double of 173 values
left = []; %only giving us the i (index) of trials that have left stimulus
right = []; %same idea with right
     
stim_onset_per_trial = [];

correct = find(choice == loc);
incorrect = find(choice ~= loc);
left_stim = find(loc == 1);
right_stim = find(loc == 2);
left_choice = find(choice == 1);
right_choice = find(choice == 2);


if savefile
    save(behavior_summary_fname, 'correct', 'incorrect', 'left_stim', 'right_stim',...
        'left_choice', 'right_choice', 'opp_contrasts');
end

%% Define important epochs
load(options.b_file_name);
earlyCue = data.response.earlyCueTime(1:ntrials)';
ts = data.response.trialstart(1:ntrials)';
rewardtime = nan(size(ts));
rewardtime(data.response.reward > 0) = data.response.rewardtimes';
rewardtime = rewardtime(1:ntrials);
endtime = data.response.trialtime(1:ntrials)';
responsetime = cellfun(@(x) x(end), data.response.timePC)';
outcometime = ts + responsetime(1:ntrials);

% Convert to frame number
ixCue = find_ix_frames(frametime, earlyCue);
ixStart = find_ix_frames(frametime, ts);
ixReward = find_ix_frames(frametime, rewardtime);
ixReward(isnan(rewardtime)) = nan;
ixOutcome = find_ix_frames(frametime, outcometime);
ixEnd = find_ix_frames(frametime, endtime);


%% Save
save(epoch_filename, 'ixCue', 'ixStart', 'ixReward', 'ixOutcome', 'ixEnd');

fprintf('Finished processing date %s\n', date);

%% Helper functions
function ix = find_ix_frames(frametime, tlist)
% frametime: array of frame times
% tlist: array of time points
% Returns the index of frames specified in tlist
    dum = repmat(frametime,numel(tlist), 1);
    [~, ix] = min(abs(dum - tlist),[], 2);
end

