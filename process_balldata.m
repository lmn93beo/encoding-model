load behavior_file_TB41.mat
%1 = r, 2 = l stim location
ntrials = 180;
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

%% Plot raw dX traces
plot(dX)
hold on;
plot([1 length(dX)], [15 15], 'k--');
plot([1 length(dX)], [-15 -15], 'k--');

for i = 1:length(samples_start)
    if data.response.choice(i) == 1
        plot([samples_start(i) samples_start(i)], [-15 15], 'r');
    elseif data.response.choice(i) == 2
        plot([samples_start(i) samples_start(i)], [-15 15], 'g');
    elseif data.response.choice(i) == 5
        plot([samples_start(i) samples_start(i)], [-15 15], 'k');
    end
end

locs = data.stimuli.loc(1:length(samples_start));
lstarts = samples_start(locs == 1);
rstarts = samples_start(locs == 2);
plot(lstarts, ones(1, length(lstarts)) * 20, 'ro');    
plot(rstarts, ones(1, length(rstarts)) * 20, 'go'); 
    
  
    