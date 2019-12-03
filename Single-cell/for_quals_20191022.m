%% Load the necessary files
folder = 'C:\Users\Sur lab\Dropbox (MIT)\Sur\Physiology_analysis\ForRam';
epochs_filename = 'TB146_20190807_epochs';
f_filename = 'Fall_20190807';
b_filename = 'TB146_20190807_behavior_summary';
load(fullfile(folder, epochs_filename));
load(fullfile(folder, f_filename));
load(fullfile(folder, b_filename));

% Extract trials
dff = F - 0.7 * Fneu;
dff = dff(iscell(:,1) > 0, :);
trange = [-10, 32];
tvals = trange(1) : trange(2);
framesID = tvals + ixStart;
framesID = framesID';
trialsmat = dff(:, framesID(:));

ncells = sum(iscell(:,1));

% Size ncells x T x ntrials
trialsmat_cell = reshape(trialsmat, ncells, 43, []);
[ncells, T, ntrials] = size(trialsmat_cell);

easy = find(opp_contrasts < 0.2);
hard = find(opp_contrasts > 0.2);


%% Four-way split
correct_left = intersect(correct, left_choice);
correct_right = intersect(correct, right_choice);
incorrect_left = intersect(incorrect, left_choice);
incorrect_right = intersect(incorrect, right_choice);

incorrect_easy = intersect(incorrect, easy);
incorrect_hard = intersect(incorrect, hard);

% Split
CL_group = squeeze(trialsmat_cell(:,:,correct_left));
CR_group = squeeze(trialsmat_cell(:,:,correct_right));
IL_group = squeeze(trialsmat_cell(:,:,incorrect_left));
IR_group = squeeze(trialsmat_cell(:,:,incorrect_right));

IE_group = squeeze(trialsmat_cell(:,:,incorrect_easy));
IH_group = squeeze(trialsmat_cell(:,:,incorrect_hard));


%% Plot all cells
%colorCL = 
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
end

% Average
CL_mean = mean(CL_group, 3);
CR_mean = mean(CR_group, 3);
IL_mean = mean(IL_group, 3);
IR_mean = mean(IR_group, 3);
IE_mean = mean(IE_group, 3);
IH_mean = mean(IH_group, 3);

%% Plot one cell
i = 121;
figure;
CL_cell = squeeze(CL_group(i, :, :));
CR_cell = squeeze(CR_group(i, :, :));
IL_cell = squeeze(IL_group(i, :, :));
IR_cell = squeeze(IR_group(i, :, :));
l1 = stdshade(CL_cell', 0.5, 'r', tvals / 12.2, []);
hold on
l2 = stdshade(CR_cell', 0.5, 'b',  tvals / 12.2,[]);
l3 = stdshade(IL_cell', 0.5, 'r--', tvals / 12.2,[]);
l4 = stdshade(IR_cell', 0.5, 'b--',  tvals / 12.2,[]);
xlabel('Time (s)', 'FontSize', 16);
ylabel('% dF/F', 'FontSize', 16);
set(gca, 'FontSize', 16);
%plot([0, 0], [60 200], 'k--', 'LineWidth', 1.5);
legend([l1, l2, l3, l4], {'Correct L', 'Correct R', 'Error L', 'Error R'})


%% Statistics
meanCL_pre = mean(CL_cell(1:10, :), 1);
meanCL_post = mean(CL_cell(11:30, :), 1);
CL_stim = ones(1, size(CL_cell, 2));
CL_choice = ones(1, size(CL_cell, 2));

meanCR_pre = mean(CR_cell(1:10, :), 1);
meanCR_post = mean(CR_cell(11:30, :), 1);
CR_stim = ones(1, size(CR_cell, 2)) * (-1);
CR_choice = ones(1, size(CR_cell, 2)) * (-1);

meanIL_pre = mean(IL_cell(1:10, :), 1);
meanIL_post = mean(IL_cell(11:30, :), 1);
IL_stim = ones(1, size(IL_cell, 2)) * (-1);
IL_choice = ones(1, size(IL_cell, 2));

meanIR_pre = mean(IR_cell(1:10, :), 1);
meanIR_post = mean(IR_cell(11:30, :), 1);
IR_stim = ones(1, size(IR_cell, 2));
IR_choice = ones(1, size(IR_cell, 2)) * (-1);

posts = [meanCL_post, meanCR_post, meanIL_post, meanIR_post]';
stims = [CL_stim, CR_stim, IL_stim, IR_stim]';
choice = [CL_choice, CR_choice, IL_choice, IR_choice]';
tbl = table(posts, stims, choice);
tbl.stims = categorical(tbl.stims);
tbl.choice = categorical(tbl.choice);

mdl = fitlm(tbl, 'posts ~ stims + choice + stims * choice')

% scatter([meanCR_post meanIR_post], [meanIR_post )
% hold on
% plot(meanCL_post, meanCL_post)
%% ANOVA
p = anovan(tbl.posts, {tbl.stims, tbl.choice},'model','interaction','varnames',{'tbl.stims','tbl.choice'});

%% Make cell array for all cells
for i = 1:ncells
    lines = zeros(T, 2);
    lines(:, 1) = IE_mean(i,:);
    lines(:, 2) = IH_mean(i,:);
%     lines(:, 3) = IL_mean(i,:);
%     lines(:, 4) = IR_mean(i,:);
    cellarr{i} = lines;
end
%%
figure;
cellID = 3;
plot(CL_mean(cellID,:))
hold on;
plot(CR_mean(cellID,:))
plot(IL_mean(cellID,:))
plot(IR_mean(cellID,:))

%% Plot all
plotIndividuals(cellarr);












