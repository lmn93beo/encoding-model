function [trials_dff, trials_z_dff, dff, z_dff, ix, ixCue] = getTrials_tb(options) 
%ONLY FOR PROTOCOL FILES GENERATED WITH TRACKBALL_BEHAVIOR SCRIPT

%% Load options
f_folder_name = options.f_folder_name; %folder name for flourescence files
b_file_name = options.b_file_name; %full behavior filename, including folder
neuropil = options.neuropil; %set true if ROIs include neuropil (cells are odd ROIs, corresponding neuropil even ROI)
neuropil_subt = options.neuropil_subt; %set true if want to subtract neuropil from cell flourescence
dt = options.dt; %time window for chunking trials

%% Load Fluoresence data, subtract neuropil if specified, make f into DFF
cd(f_folder_name);
f = importdata('F_wNeuropil_partial_tb41_03032017.txt');
f = f.data(:,2:end)';
if neuropil & neuropil_subt
    neuropil_f = f(2:2:end,:);
    raw_f = f(1:2:end,:) - 0.7*neuropil_subt; %%
elseif neuropil
    raw_f = f(1:2:end,:);
else
    raw_f = f;
end

for c = 1:size(raw_f,1)
        [freq xi] = ksdensity(raw_f(c,:));
        [~,idx] = max(freq);
        baseline = xi(idx);
        dff(c,:) = (raw_f(c,:)-baseline)/baseline*100;
        z_dff(c,:) = zscore(dff(c,:));
end

%% Load Prairie XML file to get exact frame times
cfgfiles = [dir('*.xml')];
if isempty(cfgfiles)
    disp('Could not find config file, using default framerate')
    keyboard
else
    raw_xml = fileread(cfgfiles.name);
end

sep = strfind(raw_xml,'<Frame relativeTime');
frameNum = numel(sep);
for i = 1:frameNum
    curr_idx = strfind(raw_xml(sep(i):sep(i)+100),'"');
    toGet = sep(i)+curr_idx(1);
    if i == 1 %because first entry for "relativetime" is 0
        frametime(i) = str2num(raw_xml(1,toGet));
    else
        frametime(i) = str2num(raw_xml(1,toGet:toGet+8));
    end
end

exactPeriod = (frametime(end) - frametime(1))/frameNum;
avgFR = 1/exactPeriod;

%% Split into trials
load(b_file_name);
ts = data.response.trialstart';
earlyCue = data.response.earlyCueTime';

dum = repmat(frametime,numel(ts),1);
[~,ix] = min(abs(dum - ts),[],2);
[~,ixCue] = min(abs(dum - earlyCue), [], 2);
% ix gives the 

dt_frame = round(dt*avgFR);

ix_range = repmat(dt_frame(1):dt_frame(2),numel(ix),1); 
ix_range = ix_range + ix;
ix_range(ix_range < 1) = 1;

nCells = size(dff,1);
nTrials = size(ix_range,1);

for i = 1:nCells
    for ii = 1:nTrials
        try %accounts for stimulus script running longer than imaging session duration
            trials_dff{i}(ii,:) = dff(i,ix_range(ii,:));
            trials_z_dff{i}(ii,:) = z_dff(i,ix_range(ii,:));
        end
    end
end