%% working with 
% all_left 
% all_right 
% all_one
% all_two 
% all_three 
% all_four 
% all_reward 
% all_cue
% all_pr
% all_pw 

function [peak_storage, x_vals, y_vals] = find_peaks(data)
    peak_storage = cell(180, 1);
    x_vals = [];
    y_vals = [];
    for i = 1:length(data)
        t_start = 5;
        t_end = 15;
        timeframe = t_start:t_end;
        [y, x] = findpeaks(data(timeframe, i), timeframe) ;
        
        %x vals are the timeframe, y vals are the neural act
        peak_storage{i} = {x, y};
        x_vals = [x_vals, x];
        y_vals = [y_vals, y'];
    end
    
        