
%take average across all the trials
%average of each 

%thinking of 173 (every single neural activity) x 10 (each timeframe)
%(taking the mean/variance of each column)

%how would you do this? looping through - 
%preset data you need to worry about - neural_act_matrix, cropped 1:10

%thinking of rewriting code as looping through 173 trials
%neural_act_mat{i}(j, 1:1)
    %looping witihin is gonna be j
    %think about a cell of the 10 rows, adding more to each spot
  
%% reorganization with cue_cells

cue_cells = cell(10, 1); %within each cell is an added double of 173 values
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

for i = 1:173 %looping through all the trials
    
    for j = 1:10
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
    
    stim_onset_per_trial = [stim_onset_per_trial, mid];
    
    %tracking which trials are correct
    if any(rewardsCell{i}) %if a reward is presented
        correct = [correct, i];
        if ~(i + 1 == 174);
            prev_right = [prev_right, i + 1];
        end
    else %when there's no reward
        incorrect = [incorrect, i];
        if i~173;
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

%% difficulty levels

one_array = loop_across(neural_act_mat, one, stim_onset_per_trial);
two_array = loop_across(neural_act_mat, two, stim_onset_per_trial);
three_array = loop_across(neural_act_mat, three, stim_onset_per_trial);
four_array = loop_across(neural_act_mat, four, stim_onset_per_trial);

%% worrying about previously correct (prev_right + prev_wrong)
%timeframe is going to be the beginning

prev_right_array = loop_across(neural_act_mat, prev_right, zeros(length(prev_right), 1) + 6);
prev_wrong_array = loop_across(neural_act_mat, prev_wrong, zeros(length(prev_wrong), 1) + 6);
%%
reward = correct; %trial # that were correct and thus had a reward

%find reward_onset_per_trial through looping through reward D:
reward_onset_per_trial = [];

for i=1:length(reward)
     time_of_reward = find(rewardsCell{reward(i)});
    reward_onset_per_trial = [reward_onset_per_trial, time_of_reward];
end

reward_array = loop_across(neural_act_mat, reward, reward_onset_per_trial);

%% organizing choice

correct_array = loop_across(neural_act_mat, correct, stim_onset_per_trial);
incorrect_array = loop_across(neural_act_mat, incorrect, stim_onset_per_trial);
cue_array = cell2mat(cue_cells);

% calculating the mean - matter of doing mean(cue_array, 2)
% and variances var(cue_array, 2)
 %% for left and right stimulus
     
left_cells = cell(10, 1);
right_cells = cell(10, 1);

for i = 1:length(left)
    mid = find(left_onsetCells{left(i)});
    start = mid - 6;
    for j = 1:10
        left_cells{j} = [left_cells{j}, neural_act_mat{left(i)}(start + j)];
    end
end

for i = 1:length(right)
    mid = find(right_onsetCells{right(i)});
    start = mid - 6;
    for j = 1:10
        right_cells{j} = [right_cells{j}, neural_act_mat{right(i)}(start + j)];
    end
end

left_array = cell2mat(left_cells);
right_array = cell2mat(right_cells);


%% our important array values are
% left_array  x
% right_array x

% correct_array  x
% incorrect_arrayx
% cue_array      x
% reward_array   x
% prev_right_array
% prev_wrong_array

% one_array   x
% two_array   x
% three_array x
% four_array  x

% for each array, the size is going to be 10 x how many trials the variable
% corresponds to 

%PLOTTING THE STIMULUS ONSET
% 6 is the stimulus onset point
%left_array
figure(1)
plot(left_array)
title('Neural Activity from Left Stimulus')
xlabel('Timeframe')
ylabel('Neural Activity')
hold on
plot([6 6],[-2 8])
hold on
errorbar(1:10, mean(left_array, 2), std(left_array, 0, 2), 'r')

%right_array
figure(2)
plot(right_array)
title('Neural Activity from Right Stimulus')
xlabel('Timeframe')
ylabel('Neural Activity')
hold on
plot([6 6],[-2 8])
hold on
errorbar(1:10, mean(right_array, 2), std(right_array, 0, 2), 'r')

%one_array
figure(3)
plot(one_array)
title('Neural Activity from First Level Difficulty')
xlabel('Timeframe')
ylabel('Neural Activity')
hold on
plot([6 6],[-2 8])
hold on
errorbar(1:10, mean(one_array, 2), std(one_array, 0, 2), 'r')

%two_array
figure(4)
plot(two_array)
title('Neural Activity from Second Level Difficulty')
xlabel('Timeframe')
ylabel('Neural Activity')
hold on
plot([6 6],[-2 8])
hold on
errorbar(1:10, mean(two_array, 2), std(two_array, 0, 2), 'r')

%three_array
figure(5)
plot(three_array)
title('Neural Activity from Third Level Difficulty')
xlabel('Timeframe')
ylabel('Neural Activity')
hold on
plot([6 6],[-2 8])
hold on
errorbar(1:10, mean(three_array, 2), std(three_array, 0, 2), 'r')

%four_array
figure(6)
plot(four_array)
title('Neural Activity from Fourth Level Difficulty')
xlabel('Timeframe')
ylabel('Neural Activity')
hold on
plot([6 6],[-2 8])
hold on
errorbar(1:10, mean(four_array, 2), std(four_array, 0, 2), 'r')

%REWARD ONSET
%reward_array
figure(7)
plot(reward_array)
title('Neural Activity from Reward')
xlabel('Timeframe')
ylabel('Neural Activity')
hold on
plot([6 6],[-2 8])
hold on
errorbar(1:10, mean(reward_array, 2), std(reward_array, 0, 2), 'r')

%CUE ONSET
%cue_array
figure(8)
plot(cue_array)
title('Neural Activity from Cue')
xlabel('Timeframe')
ylabel('Neural Activity')
hold on
errorbar(1:10, mean(cue_array, 2), std(cue_array, 0, 2), 'r')

%prev_right_array
figure(9)
plot(prev_right_array)
title('Neural Activity from Previously Correct Trials')
xlabel('Timeframe')
ylabel('Neural Activity')
hold on 
errorbar(1:10, mean(prev_right_array, 2), std(prev_right_array, 0, 2), 'r')

%prev_wrong_array
figure(10)
plot(prev_wrong_array)
title('Neural Activity from Previously Incorrect Trials')
xlabel('Timeframe')
ylabel('Neural Activity')
hold on 
errorbar(1:10, mean(prev_wrong_array, 2), std(prev_wrong_array, 0, 2), 'r')

