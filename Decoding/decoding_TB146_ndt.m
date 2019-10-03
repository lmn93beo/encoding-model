

%% Add path to ndt toolbox
clear all
pathname = 'C:\Users\Sur lab\Dropbox (MIT)\Sur\ExternalCode\ndt.1.0.4\ndt.1.0.4\';
encodingpath = 'C:\Users\Sur lab\Dropbox (MIT)\Sur\ExternalCode\encoding-model';
addpath(pathname);
addpath(encodingpath);
add_ndt_paths_and_init_rand_generator;

rand('state', 123);



%% Load the data
dates = {'20190724', '20190725', '20190726', '20190727',...
    '20190730', '20190801', '20190806', '20190807'};
decoding_type = 'reward';
align_by = 'outcome';

for i = 1:numel(dates)
    fprintf('Doing file %d of %d...\n', i, numel(dates));
    date = dates{i};
    do_decoding(date, decoding_type, align_by);
end


function do_decoding(date, decoding_type, align_by)
options.f_folder_name = sprintf('C:\\Users\\Sur lab\\Dropbox (MIT)\\Sur\\Physiology_analysis\\146\\%s\\suite2p\\plane0', date);
options.b_file_name = sprintf('C:\\Users\\Sur lab\\Dropbox (MIT)\\trackball-behavior\\Data\\146_all\\%s_trackball_0146.mat', date);
root_filename = 'TB146';
encoding_struct_fname = [root_filename '_encoding_structs'];
behavior_summary_fname = [root_filename '_behavior_summary'];
savefile = 1;
options.f_file_name = 'F_wNeuropil_partial_tb41_03032017.txt';
options.special = 0;
options.freq = 12.2;
options.align_by = align_by;

if strcmp(date, '20190724')
    options.special = 1;
else
    options.special = 0;
end

options.neuropil = 1;
options.neuropil_subt = 1;
options.dt = [-3 5];
options.suite2p = 1;

[trials_dff, ~, dff, ~, ~, ~, ~] = getTrials_tb(options);

load(sprintf('C:\\Users\\Sur lab\\Dropbox (MIT)\\Sur\\ExternalCode\\encoding-model\\Python_analysis\\TB146_%s_epochs.mat', date));


%% Make the binned labels
load(options.b_file_name)
ncells = numel(trials_dff);

trange = [-50, 100];

%1 = 2, r = 1 stim location
validTrials = find(ixStart < size(dff, 2) - trange(2) & ixStart > -trange(1));

[binned_labels, lrvalid] = make_labels(data, decoding_type, ncells, validTrials);


%% Make the binned data
binned_data = cell(1, ncells);

for i = 1:ncells
    binned_data{i} = trials_dff{i}(lrvalid,:);
end

period = 1000 / 12.2; %ms
binned_site_info.binning_parameters.bin_width = period;
binned_site_info.binning_parameters.sampling_interval = period;
binned_site_info.binning_parameters.start_time = -3000 - period/2;
binned_site_info.binning_parameters.end_time = 5000 + period;


binned_data_file_name = strcat('TB146_decoding\TB146_', decoding_type, '_raw_', date);
save(binned_data_file_name, 'binned_data', 'binned_labels', 'binned_site_info');

%%
for i = 0:100
    inds_of_sites_with_at_least_k_repeats = find_sites_with_k_label_repetitions(binned_labels.side, i);
    num_sites_with_k_repeats(i + 1) = length(inds_of_sites_with_at_least_k_repeats);
end



%%  Begin the decoding analysis  %%


%%  6.  Create a datasource object

% we will use object identity labels to decode which object was shown (disregarding the position of the object)
specific_binned_labels_names = 'side';

% use 20 cross-validation splits (which means that 19 examples of each object are used for training and 1 example of each object is used for testing)
% Confirm that we have enough trials
nsites = num_sites_with_k_repeats(1);
num_cv_splits = 20;
while num_sites_with_k_repeats(num_cv_splits + 1) < nsites
    num_cv_splits = ceil(num_cv_splits / 2);
end
fprintf('Num cv splits = %d\n', num_cv_splits);
% create the basic datasource object
ds = basic_DS(binned_data_file_name, specific_binned_labels_names,  num_cv_splits);


%%   7.  Create a feature preprocessor object

% create a feature preprocess that z-score normalizes each feature
the_feature_preprocessors{1} = zscore_normalize_FP;  



%%  8.  Create a classifier object 

% select a classifier
the_classifier = max_correlation_coefficient_CL;


% other useful options:   

% use a poisson naive bayes classifier (note: the data needs to be loaded as spike counts to use this classifier)
%the_classifier = poisson_naive_bayes_CL;  

% use a support vector machine (see the documentation for all the optional parameters for this classifier)
%the_classifier = libsvm_CL;


%%  9.  create the cross-validator 


the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);  

the_cross_validator.num_resample_runs = 2;  % usually more than 2 resample runs are used to get more accurate results, but to save time we are using a small number here


% other useful options:   

% can greatly speed up the run-time of the analysis by not creating a full TCT matrix (i.e., only trainging and testing the classifier on the same time bin)
the_cross_validator.test_only_at_training_times = 1;  




%%  10.  Run the decoding analysis   

% if calling the code from a script, one can log the code so that one can recreate the results 
%log_code_obj = log_code_object;
%log_code_obj.log_current_file; 


% run the decoding analysis 
DECODING_RESULTS = the_cross_validator.run_cv_decoding; 



%%  11.  Save the results

% save the results
save_file_name = strcat('TB146_decoding\TB146' , decoding_type, '_decoding_results_', date);
save(save_file_name, 'DECODING_RESULTS'); 

% if logged the code that was run using a log_code_object, save the code
%LOGGED_CODE = log_code_obj.return_logged_code_structure;
%save(save_file_name, '-v7.3', 'DECODING_RESULTS', 'LOGGED_CODE'); 



%%  12.  Plot the basic results


% which results should be plotted (only have one result to plot here)
result_names{1} = save_file_name;

% create an object to plot the results
plot_obj = plot_standard_results_object(result_names);

% put a line at the time when the stimulus was shown
plot_obj.significant_event_times = 0;   


% optional argument, can plot different types of results
%plot_obj.result_type_to_plot = 2;  % for example, setting this to 2 plots the normalized rank results


plot_obj.plot_results;   % actually plot the results
title(strcat(decoding_type, ' decoding, date:', date));

saveas(gcf, strcat('TB146_decoding\TB146outcome_', decoding_type, '_', date, '.pdf'));
end

   
   






