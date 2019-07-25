

function [array] = loop_across(neural_mat, list, timeframe_list)
%cell_container is the cell(10, 1) container that we'll convert into cell2mat 
    cell_container = cell(15, 1);
    for i = 1:length(list)
        start = timeframe_list(i) - 6;
        for j = 1:length(cell_container)
            cell_container{j} = [cell_container{j}, neural_mat{list(i)}(start + j)];
        end
    end
    
    array = cell2mat(cell_container);


end

function = plot_each(array, timeframe, vline)
    plot(array)
    t = 'Neural Activity from %s'
    a = array
    title('Neural Activity from %s', 
    xlabel('Timeframe')
    ylabel('Neural Activity')
    hold on
    errorbar(1: timeframe, mean(array, 2), std(array, 0, 2), 'r')
end

formatSpec = 'The array is %dx%d.';
A1 = 2;
A2 = 3;
str = sprintf(formatSpec,A1,A2)