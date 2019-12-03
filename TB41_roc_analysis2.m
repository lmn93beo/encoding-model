options.f_folder_name = 'C:\Users\Sur lab\Dropbox (MIT)\Sur\ExternalCode\encoding-model';
options.b_file_name = 'C:\Users\Sur lab\Dropbox (MIT)\Sur\ExternalCode\encoding-model/behavior_file_TB41.mat';
options.f_file_name = 'F_wNeuropil_partial_tb41_03032017.txt';

options.neuropil = 1;
options.neuropil_subt = 1;
options.dt = [-3 5];
options.suite2p = 0;
options.special = 0;
options.align_by = 'trial_start';

[trials_dff, trials_z_dff, dff, z_dff, frametimes, ix, ixCue] = getTrials_tb(options);
ncells = size(z_dff, 1);

%ntrials = size(fdata, 3);

load behavior_file_TB41.mat

%% Plot an example cell
cellid = 6;
celldata = trials_z_dff{cellid};
meanresp = mean(celldata, 3);
%plot(meanresp)

%% Behavior
% Correct left & correct right trials
ntrials = numel(data.response.choice);
CL_trials = find(data.response.choice(1:257) == 1 & data.stimuli.loc(1:257) == 1);
CR_trials = find(data.response.choice(1:257) == 2 & data.stimuli.loc(1:257) == 2);

celldata_CL = celldata(CL_trials, :);
celldata_CR = celldata(CR_trials, :);

mean_CL = mean(celldata_CL, 1);
mean_CR = mean(celldata_CR, 1);

figure;
plot(mean_CL')
hold on
plot(mean_CR')

%% Plot all cells
for i = 1:100
    subplot(10, 10, i)
    CL_cell = squeeze(CL_group(i, :, :));
    CR_cell = squeeze(CR_group(i, :, :));
    IL_cell = squeeze(IL_group(i, :, :));
    IR_cell = squeeze(IR_group(i, :, :));
    stdshade(CL_cell', 0.5, 'r--', []);
    hold on
    stdshade(CR_cell', 0.5, 'b',  []);
    stdshade(IL_cell', 0.5, 'b--', []);
    stdshade(IR_cell', 0.5, 'r',  []);
    title(i)
end

