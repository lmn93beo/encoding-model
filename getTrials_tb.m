function [trials_dff, trials_z_dff, dff, z_dff, frametime, ix, ixCue] = getTrials_tb(options) 
%ONLY FOR PROTOCOL FILES GENERATED WITH TRACKBALL_BEHAVIOR SCRIPT

%% Load options
f_folder_name = options.f_folder_name; %folder name for flourescence files
b_file_name = options.b_file_name; %full behavior filename, including folder
neuropil = options.neuropil; %set true if ROIs include neuropil (cells are odd ROIs, corresponding neuropil even ROI)
neuropil_subt = options.neuropil_subt; %set true if want to subtract neuropil from cell flourescence
dt = options.dt; %time window for chunking trials

%% Load Fluoresence data, subtract neuropil if specified, make f into DFF
if ~options.suite2p
    f = importdata(options.f_file_name);
    f = f.data(:,2:end)';
    if neuropil & neuropil_subt
        neuropil_f = f(2:2:end,:);
        raw_f = f(1:2:end,:) - 0.7*neuropil_subt; %%
    elseif neuropil
        raw_f = f(1:2:end,:);
    else
        raw_f = f;
    end
else
    load([f_folder_name '\Fall.mat'], 'F', 'Fneu', 'iscell');

    raw_f = F - Fneu * 0.7;
    raw_f = raw_f(logical(iscell(:,1)), :);
end

for c = 1:size(raw_f,1)
        [freq xi] = ksdensity(raw_f(c,:));
        [~,idx] = max(freq);
        baseline = xi(idx);
        dff(c,:) = (raw_f(c,:)-baseline)/baseline*100;
        z_dff(c,:) = zscore(dff(c,:));
end

%% Load Prairie XML file to get exact frame times
cfgfiles = [dir([options.f_folder_name '\*.xml'])];
if isempty(cfgfiles)
    disp('Could not find config file, using default framerate')
    keyboard
else
    assert(numel(cfgfiles) == 1);
    raw_xml = fileread([options.f_folder_name '\' cfgfiles.name]);
end

sep = strfind(raw_xml,'<Frame relativeTime');
if options.special
    frameNum = 15000; %numel(sep); % MANUAL! WILL CHANGE
else
    frameNum = numel(sep);
end

for i = 1:frameNum
    curr_idx = strfind(raw_xml(sep(i):sep(i)+100),'"');
    toGet = sep(i)+curr_idx(1);
    if i == 1 %because first entry for "relativetime" is 0
        frametime(i) = str2num(raw_xml(1,toGet));
    else
        if i == 15000
            disp(i);
        end
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