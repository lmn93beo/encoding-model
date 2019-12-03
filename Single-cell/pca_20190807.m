%% Load the necessary files
folder = 'C:\Users\Sur lab\Dropbox (MIT)\Sur\Physiology_analysis\ForRam';
epochs_filename = 'TB146_20190727_epochs';
f_filename = 'Fall_20190727';
b_filename = 'TB146_20190727_behavior_summary';
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

%% Split into conditions
correct_left_easy = intersect(intersect(correct, left_choice), easy);
correct_left_hard = intersect(intersect(correct, left_choice), hard);
correct_right_easy = intersect(intersect(correct, right_choice), easy);
correct_right_hard = intersect(intersect(correct, right_choice), hard);

correct_left = intersect(correct, left_choice);
correct_right = intersect(correct, right_choice);
incorrect_left = intersect(incorrect, left_choice);
incorrect_right = intersect(incorrect, right_choice);

%% Get appropriate arrays
arr_LE = trialsmat_cell(:, :, correct_left);
arr_LH = trialsmat_cell(:, :, correct_right);
arr_RE = trialsmat_cell(:, :, incorrect_left);
arr_RH = trialsmat_cell(:, :, incorrect_right);

% Average
mean_LE = mean(arr_LE, 3);
mean_LH = mean(arr_LH, 3);
mean_RE = mean(arr_RE, 3);
mean_RH = mean(arr_RH, 3);

% Concatenate
mean_all = [mean_LE mean_LH mean_RE mean_RH];

% z-score
mean_cells = mean(mean_all, 2);
std_cells = std(mean_all, [], 2);
mean_z = (mean_all - mean_cells) ./ std_cells;

%% PCA
[U, S, V] = svd(mean_z');

% Get top K components
K = 5;
z_pca = U(:,1:K) * S(1:K,1:K) * V(:, 1:K)';


col1 = V(:, 1);
col2 = V(:, 2);
col3 = V(:, 3);

proj1 = col1' * mean_z;
proj2 = col2' * mean_z;
proj3 = col3' * mean_z;

l1 = plot3(proj1(1:43), proj2(1:43), proj3(1:43), 'b');
hold on;
l2 = plot3(proj1(44:86), proj2(44:86), proj3(44:86), 'r');
l3 = plot3(proj1(87:129), proj2(87:129), proj3(87:129), 'b--');
l4 = plot3(proj1(130:end), proj2(130:end), proj3(130:end), 'r--');

plot3(proj1(1), proj2(1), proj3(1), 'o');
plot3(proj1(44), proj2(44), proj3(44), 'o');
plot3(proj1(87), proj2(87), proj3(87), 'o');
plot3(proj1(130), proj2(130), proj3(130), 'o');

legend([l1, l2, l3, l4], {'mean_CL', 'mean_CR', 'mean_IL', 'mean_IR'})










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
    stdshade(IL_cell', 0.5, 'b--', []);
    stdshade(IR_cell', 0.5, 'r',  []);
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
i = 61;
figure;
CL_cell = squeeze(CL_group(i, :, :));
CR_cell = squeeze(CR_group(i, :, :));
IL_cell = squeeze(IL_group(i, :, :));
IR_cell = squeeze(IR_group(i, :, :));
l1 = stdshade(CL_cell', 0.5, 'r--', tvals / 12.2, []);
hold on
l2 = stdshade(CR_cell', 0.5, 'b',  tvals / 12.2,[]);
l3 = stdshade(IL_cell', 0.5, 'b--', tvals / 12.2,[]);
l4 = stdshade(IR_cell', 0.5, 'r',  tvals / 12.2,[]);
xlabel('Time (s)', 'FontSize', 16);
ylabel('% dF/F', 'FontSize', 16);
set(gca, 'FontSize', 16);
%plot([0, 0], [60 200], 'k--', 'LineWidth', 1.5);
legend([l1, l2, l3, l4], {'Stim L, Choice L', 'Stim R, Choice R', 'Stim R, Choice L', 'Stim L, Choice R'})



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












