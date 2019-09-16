load behavior_file_TB41.mat

%%
% Convert cells to mat
precueS = cell2mat(data.response.precue_samples_start);
startS = cell2mat(data.response.samples_start);
rewardS = cell2mat(data.response.samples_reward);

% plot
precueT = data.response.time(precueS);
startT = data.response.time(startS);
rewardT = data.response.time(rewardS);

% Plot
plot(startT - precueT)
hold on
plot(data.response.trialstart - data.response.earlyCueTime);

%%
plot(data.response.trialstart - startT')


%% Make reward times
rewardt = zeros(numel(startT), 1);
rewardt(data.response.reward > 0) = data.response.rewardtimes - data.response.trialstart(data.response.reward > 0);
plot(rewardt)
hold on;
plot(rewardT - startT - 4 + 1)

