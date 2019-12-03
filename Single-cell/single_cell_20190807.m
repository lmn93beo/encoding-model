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

% Size ncells x T x ntrials
trialsmat_cell = reshape(trialsmat, 152, 43, 275);
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

% Average
CL_mean = mean(CL_group, 3);
CR_mean = mean(CR_group, 3);
IL_mean = mean(IL_group, 3);
IR_mean = mean(IR_group, 3);
IE_mean = mean(IE_group, 3);
IH_mean = mean(IH_group, 3);

%% Trying stdshade
CL_cell = squeeze(CL_group(1, :, :));
stdshade(CL_cell');


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












