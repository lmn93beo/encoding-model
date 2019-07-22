

function [array] = loop_across(neural_mat, list, timeframe_list)
%cell_container is the cell(10, 1) container that we'll convert into cell2mat 
    cell_container = cell(10, 1);
    for i = 1:length(list)
        start = timeframe_list(i) - 6;
        for j = 1:length(cell_container)
            cell_container{j} = [cell_container{j}, neural_mat{list(i)}(start + j)];
        end
    end
    
    array = cell2mat(cell_container);


end
