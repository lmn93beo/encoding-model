load tb_preprocessed_TB41.mat
fdata = trial_data{2};
T = size(fdata, 2);
id = trial_id{2}; %vis stim
sortid = trial_sortid{2};
curr_dir = pwd;

ntrials = size(fdata, 3);

%%
load behavior_file_TB41.mat
%1 = r, 2 = l stim location
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

modified_RT = data.response.reward * 0;
modified_RT(data.response.reward > 0) = data.response.rewardtimes;
dt_rew = modified_RT - data.response.trialstart;
dt_rew(~data.response.reward) = 0;
dt_rew_frames = floor(dt_rew / 0.2);

% For ball movement
samples_start = cell2mat(data.response.samples_start);
dX = data.response.dx;
timepts = data.response.time;
fs = (length(timepts) - 1) / (timepts(end) - timepts(1));
spacings = floor((1:41) * fs);
% Center at frame 15
spacings = spacings - spacings(15);

%%
% Subsample the trials with opp_contrast == 0.32 and choice ~= 5
goodtrials = (choice(2:end) ~= 5) & (choice(1:end-1) ~= 5);
goodtrials = logical([0 goodtrials]);
%goodtrials = (choice ~= 5);
%goodtrials(1) = 0;
%goodtrials(258) = 0; %exclude this trial because neural activity has nan's
%TODO: consider not excluding this

lstim = lstim(goodtrials);
rstim = rstim(goodtrials);
loc = loc(goodtrials);
choice = choice(goodtrials);
opp_contrasts = opp_contrasts(goodtrials);
goodtrialsID = find(goodtrials);
prevchoice = data.response.choice(goodtrialsID - 1);

difficulty = abs(lstim - rstim);
samples_start = samples_start(goodtrials);



%% Make the predictors
ntrialsgood = sum(goodtrials);

% Kernel for L-stim onset
left_onsetCells = {};
for i = 1:10
    left_onset = zeros(ntrialsgood, T);
    left_onset(loc == 2,15 + i) = 1;
    left_onsetCell = mat2cell(left_onset', 41, ones(1, ntrialsgood))';
    left_onsetCells{i} = left_onsetCell;
end

% Kernel for R-stim onset
right_onsetCells = {};
for i = 1:10
    right_onset = zeros(ntrialsgood, T);
    right_onset(loc == 1,15 + i) = 1;
    right_onsetCell = mat2cell(right_onset', 41, ones(1, ntrialsgood))';
    right_onsetCells{i} = right_onsetCell;
end

% Reward times
rewards = zeros(ntrialsgood, T);
for i = 1:ntrialsgood
    if dt_rew_frames(goodtrialsID(i)) > 0
        rewards(i, 15 + dt_rew_frames(goodtrialsID(i))) = 1;
    end    
end
rewardsCell = mat2cell(rewards', 41, ones(1, ntrialsgood))';
    
% Build the time information
balldata = cell(ntrialsgood, 1);
for i = 1:ntrialsgood
    %id = goodtrialsID(i);
    start = samples_start(i);
    samplesID = start + spacings;
    balldata{i} = dX(samplesID);
end    


pred_types_cell = {'event', 'event', 'event', 'whole-trial', 'whole-trial', 'whole-trial', 'continuous'};
groupings = {1, 2, 3, 4, 5, 6, 7};

[pred_allmat,pred_inds_cell,grouped_var_types] = make_predictor_matrix_generalcase({left_onsetCells{1}, right_onsetCells{1},...
     rewardsCell, choice, prevchoice, difficulty, balldata}, ...
    pred_types_cell, groupings);

%% Make the neural matrix (response)
% Make the activity matrix
goodtrialsID = find(goodtrials);
neural_act_mat = cell(ntrialsgood, 1);
for i = 1:ntrialsgood
    neural_act_mat{i} = fdata(:,:,goodtrialsID(i))';
end

neural_matmat = cell2mat(neural_act_mat);
predall_matmat = cell2mat(pred_allmat);

figure;
cellID = 1;
plot(neural_matmat(:,cellID));
hold on;
for i = 1:180
    plot([i * 41, i * 41], [-200 1000], '--');
end





save('TB41_predallmat.mat', 'neural_matmat', 'predall_matmat', 'neural_act_mat', 'pred_allmat');
%% Do encoding model
approach = 'norefit';
pred_types_cell_group = {'event', 'whole-trial', 'whole-trial', 'whole-trial', 'continuous'};

[relative_contrib,~,r2val] = process_encoding_model(pred_allmat, pred_inds_cell, neural_act_mat, pred_types_cell,approach);


%% Visualize
means = nanmean(relative_contrib, 1);
stds = nanstd(relative_contrib, 1) / sqrt(size(relative_contrib, 1));
bar(means);
hold on;
errorbar(1:size(relative_contrib, 2), means, stds, '.', 'Color', 'b');
%hold on;
xticks(1:size(relative_contrib, 2));

xticks(1:size(relative_contrib, 2));
factors = {'LeftStim', 'RightStim', 'Reward', 'Choice', 'PrevChoice', 'Difficulty', 'Movement'};
ylabel('Relative Contribution');


xticklabels(factors);

%% Scatter plot
figure(2);
hold on;
scatter(relative_contrib(:,1), relative_contrib(:,3), 'filled', 'MarkerEdgeColor', 'w');
xlabel('Stimulus contribution')
ylabel('Choice contribution')


%% Cumulative distribution of r2 values
%plot(sort(r2val), (1:numel(r2val)) / numel(r2val));



for i = 1:numel(rewardsCell)
   if sum(rewardsCell{i}) == 0
       fprintf('Not found\n');
   else
       disp(find(rewardsCell{i}));
   end
end