

function [array] = loop_across(neural_mat, list, timeframe_list, neuron)
%cell_container is the cell(10, 1) container that we'll convert into cell2mat 
    cell_container = cell(15, 1);
    for i = 1:length(list)
        trial_num = list(i);
        start = timeframe_list(trial_num) - 6; %if starting at beginning, start = 0
        for j = 1:length(cell_container)
            cell_container{j} = [cell_container{j}, neural_mat{trial_num}(start + j, neuron)];
        end
    end 
    array = cell2mat(cell_container);
end
