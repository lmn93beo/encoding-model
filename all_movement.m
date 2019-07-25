% container for each neuron's movement + neural activity data
all_neural_data = cell(180, 1);

for n = 1:180
    % crop all the edge zeros out ("cleaning the data")
    
    cropped_ball_data = []; % all the ball data here of the first neuron 
    neural_data = []; %corresponding neural data of the ball data

    for i = 1:length(balldata)
        if find(balldata{i}) ~0; %if there's ball data (doing this because 112 has no data)
            frames = find(balldata{i}); % cropping out the edge zeros
            s = frames(1); 
            e = frames(end);
            for j = s:e
                cropped_ball_data = [cropped_ball_data, balldata{i}(j)];
                neural_data = [neural_data, neural_act_mat{i}(j, n)]; % n represents the neurons
            end
        end
    end

    bins = cell(84, 1);
    %sorting values into its bin
       for i = 1: length(cropped_ball_data)
           val = floor(cropped_ball_data(i)) + 45; % adding 45 because indices begin at 1
           bins{val} = [bins{val}, neural_data(i)]; % val serves as the key
       end
       
    %collapsing each of the bin values into its avg val 
    for i = 1:length(bins)
        if length(bins{i}) > 1; %if val exists
            bins{i} = mean(bins{i}); 
        else % accounting for bins with no values (setting automatically to zero)
            bins{i} = 0;
        end

    end
    bin_array = cell2mat(bins);
    all_bins{n} = bin_array; %add this bin to neuron's data
end

X = -42:41; % 84 range, x-axis showing the balldata speed
all_bin_mat = cell2mat(all_bins);
plot(X, all_bin_mat)
title('Averaged Neural Activity Per Ball Velocity Bin')
xlabel('Velocity')
ylabel('Neural Activity')


    
