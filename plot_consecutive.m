% Here, we try to plot signals for consecutive trials (i.e. L-L, L-R, R-L,
% R-R etc), and see if we can 'decode' the next action from what happened
% in the previous trial (aka, if the previous action was L, is the animal
% more likely to choose L?)

options.f_folder_name = 'C:\Users\Sur lab\Dropbox (MIT)\Sur\Physiology_analysis\146\20190807\suite2p\plane0';
options.b_file_name = 'C:\Users\Sur lab\Dropbox (MIT)\Sur\ExternalCode\encoding-model/20190807_trackball_0146.mat';
save_filename = 'TB146_epochs.mat';
options.f_file_name = 'F_wNeuropil_partial_tb41_03032017.txt';

options.neuropil = 1;
options.neuropil_subt = 1;
options.dt = [-3 5];
options.suite2p = 1;

[trials_dff, trials_z_dff, dff, z_dff, ~, ix, ixCue] = getTrials_tb(options);
ncells = size(z_dff, 1);

%% Visualize
visualize = 0;
if visualize
    imagesc(z_dff);
    %plot(dff(1,:));
    hold on;
    for i = 1:numel(ix)
       plot([ix(i) ix(i)], [-100 700], 'b');
       plot([ixCue(i) ixCue(i)], [-100 700], 'r');
    end
end

%% Load and preprocess behavior file
load(options.b_file_name)

%1 = 1, r = 2 stim location
ntrials = sum(ix < 9000);
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
% Find trials with L-L, L-R etc
locnext = [loc(2:end) nan];
lltrials = find(loc == 1 & locnext == 1);
lrtrials = find(loc == 1 & locnext == 2);
rltrials = find(loc == 2 & locnext == 1);
rrtrials = find(loc == 2 & locnext == 2);

%% Plot
% For each neuron, average across these trials...
for neur_id = 1:2
    figure;
    window = [-10, 18];
    %subplot(4, 1, 1)
    h1 = plot_trial_pair(lltrials, window, ix, neur_id, z_dff, 'r');
    hold on;
    title('L-L')

    %subplot(4, 1, 2)
    h2 = plot_trial_pair(lrtrials, window, ix, neur_id, z_dff, 'b');
    title('L-R')

    %subplot(4, 1, 3)
    h3 = plot_trial_pair(rltrials, window, ix, neur_id, z_dff, 'g');
    title('R-L')

    %subplot(4, 1, 4)
    h4 = plot_trial_pair(rrtrials, window, ix, neur_id, z_dff, 'k');
    title('R-R')
    legend([h1, h2, h3, h4], {'LL', 'LR', 'RL', 'RR'})
    
    %savefig(gcf, ['figure' num2str(neur_id)])
end

%% Define subfunctions
function h = plot_trial_pair(trial_lst, window, ix, neur_id, z_dff, color)
    traces1 = [];
    traces2 = [];

    for i = 1:numel(trial_lst)
        framestart1 = ix(trial_lst(i));
        framestart2 = ix(trial_lst(i) + 1);
        trace1 = z_dff(neur_id, framestart1 + window(1) : framestart1 + window(2));
        trace2 = z_dff(neur_id, framestart2 + window(1) : framestart2 + window(2));
        traces1(i, :) = trace1;
        traces2(i, :) = trace2;
    end
        
    h = errorbar(window(1) : window(2), mean(traces1, 1), stderr(traces1, 1), color);
    ylim([-0.5, 0.5])
    hold on;   

    errorbar(window(1) + 40 : window(2) + 40, mean(traces2, 1), stderr(traces2, 1), color);
    ylim([-0.5, 0.5]);
    
    
end

function err = stderr(mat, axis)
    err = std(mat, [], axis) / sqrt(size(mat, axis));
end


















