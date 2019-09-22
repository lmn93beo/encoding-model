% This script is for defining important behavioral epochs for each trial,
% specifically:
% - Early cue time
% - Trial start time
% - Decision time (when threshold is reached)
% - Reward time
% - Trial end time

%% Load data
options.f_folder_name = 'C:\Users\Sur lab\Dropbox (MIT)\Sur\ExternalCode\encoding-model';
options.f_file_name = 'F_wNeuropil_partial_tb41_03032017.txt';
options.b_file_name = 'C:\Users\Sur lab\Dropbox (MIT)\Sur\ExternalCode\encoding-model/behavior_file_TB41.mat';
options.suite2p = 0;
options.neuropil = 1;
options.neuropil_subt = 1;
options.dt = [-3 5];

[trials_dff, trials_z_dff, dff, z_dff, frametime, ix, ixCue] = getTrials_tb(options);
ncells = size(z_dff, 1);

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
%save('TB41_epochs.mat', 'ixCue', 'ixStart', 'ixReward', 'ixOutcome', 'ixEnd', 'goodtrials');



%% Helper functions
function ix = find_ix_frames(frametime, tlist)
% frametime: array of frame times
% tlist: array of time points
% Returns the index of frames specified in tlist
    dum = repmat(frametime,numel(tlist), 1);
    [~, ix] = min(abs(dum - tlist),[], 2);
end




