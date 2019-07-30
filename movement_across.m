
%% this is only for ONE neuron
% using the matlab histogram function yield 84 bins
% the range of the data is 83.8333

% x axis is gonna be the bins of ballspeeds
% y axis is going to be the neural activity
% can choose to either round down or up

%floor of the minimum is -44 + 45 to have it one-indexed
%floor of the maximum is 39 + 45

each_trial = cell(173, 1);

% crop all the edge zeros out ("cleaning the data")

cropped_ball_data = [];
neural_data = [];
for i = 1:length(balldata)
    if find(balldata{i}) ~= 0; %if there's ball data e
        frames = find(balldata{i});
        s = frames(1);
        e = frames(end);
        for j = s:e
            cropped_ball_data = [cropped_ball_data, balldata{i}(j)];
            neural_data = [neural_data, neural_act_mat{i}(j)];
        end
    end
end

bins = cell(84, 1);
%sorting values into its bin
for i = 1: length(cropped_ball_data)
    val = floor(cropped_ball_data(i)) + 45; %adding 45 because indices begin at 1
    bins{val} = [bins{val}, neural_data(i)];
end

%collapsing each of the bin values into its avg val 
for i = 1:length(bins)
    if length(bins{i}) > 1;
        bins{i} = mean(bins{i});
    else % accounting for bins with no values (setting automatically to zero)
        bins{i} = 0;
    end
    
end
X = [-42:41];
bin_array = cell2mat(bins);
plot(X, bin_array)
title('Averaged Neural Activity Per Ball Velocity Bin')
xlabel('Velocity')
ylabel('Neural Activity')


    
