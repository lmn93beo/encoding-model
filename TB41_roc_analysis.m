load tb_preprocessed_TB41.mat
fdata = trial_data{2};
T = size(fdata, 2);
id = trial_id{2}; %vis stim
sortid = trial_sortid{2};
curr_dir = pwd;

ntrials = size(fdata, 3);

load behavior_file_TB41.mat

%% Plot an example cell
cellid = 6;
celldata = fdata(cellid, :, 1:257);
meanresp = mean(celldata, 3);
%plot(meanresp)

%% Behavior
% Correct left & correct right trials
ntrials = numel(data.response.choice);
% CL_trials = find(data.response.choice(1:257) == 1 & data.stimuli.loc(1:257) == 1);
% CR_trials = find(data.response.choice(1:257) == 2 & data.stimuli.loc(1:257) == 2);

correct_left = find(data.response.choice(1:257) == 1 & data.stimuli.loc(1:257) == 1);
correct_right = find(data.response.choice(1:257) == 2 & data.stimuli.loc(1:257) == 2);
incorrect_left = find(data.response.choice(1:257) == 1 & data.stimuli.loc(1:257) == 2);
incorrect_right = find(data.response.choice(1:257) == 2 & data.stimuli.loc(1:257) == 1);

%incorrect_easy = intersect(incorrect, easy);
%incorrect_hard = intersect(incorrect, hard);

% Split
CL_group = squeeze(fdata(:,:,correct_left));
CR_group = squeeze(fdata(:,:,correct_right));
IL_group = squeeze(fdata(:,:,incorrect_left));
IR_group = squeeze(fdata(:,:,incorrect_right));

%IE_group = squeeze(trialsmat_cell(:,:,incorrect_easy));
%IH_group = squeeze(trialsmat_cell(:,:,incorrect_hard));


% celldata_CL = celldata(:, :, CL_trials);
% celldata_CR = celldata(:, :, CR_trials);
%%
% mean_CL = mean(celldata_CL, 3);
% mean_CR = mean(celldata_CR, 3);
% 
% figure;
% plot(mean_CL)
% hold on
% plot(mean_CR)


%% Plot all cells
figure;
for i = 1:100
    subplot(10, 10, i)
    CL_cell = squeeze(CL_group(i, :, :));
    CR_cell = squeeze(CR_group(i, :, :));
    IL_cell = squeeze(IL_group(i, :, :));
    IR_cell = squeeze(IR_group(i, :, :));
    stdshade(CL_cell', 0.5, 'r--', []);
    hold on
    stdshade(CR_cell', 0.5, 'b',  []);
    stdshade(IL_cell', 0.5, 'r', []);
    stdshade(IR_cell', 0.5, 'b--',  []);
    title(i)
    %ylim([-100 300])
end

figure;
for i = 101:180
    subplot(10, 10, i - 100)
    CL_cell = squeeze(CL_group(i, :, :));
    CR_cell = squeeze(CR_group(i, :, :));
    IL_cell = squeeze(IL_group(i, :, :));
    IR_cell = squeeze(IR_group(i, :, :));
    stdshade(CL_cell', 0.5, 'r--', []);
    hold on
    stdshade(CR_cell', 0.5, 'b--',  []);
    stdshade(IL_cell', 0.5, 'r', []);
    stdshade(IR_cell', 0.5, 'b',  []);
    title(i)
    %ylim([-100 300])
end
