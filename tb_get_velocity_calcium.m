%% Trackball analysis of movement trajectories -- effect of laser

function [ball_velocity, interp_ballX, rt, ballX_pc, choice] = tb_get_velocity_calcium(data);
%%
factor = data.response.mvmt_degrees_per_pixel; % conversion factor
rt = cellfun(@(x) x(end),data.response.timePC);
choice = data.response.choice;

if choice(end) == 5 & (numel(rt) < numel(choice))
    choice = choice(1:end-1);
    rt = rt(1:numel(choice));
end

%idx = find(rt <= 1 & choice < 5);
%n = find(data.response.reward>0,1,'last');
idx = 1:numel(choice);
time_pc = data.response.timePC(idx)';
ballX_pc = cellfun(@(x) -x,data.response.ballX(idx)','uniformoutput',0);
samps_pc = data.response.samps(idx)';
try
    laser = data.stimuli.laser(idx);
catch
    laser = ones(1,numel(idx));
end
loc = data.stimuli.loc(idx);
choice = data.response.choice(idx);
rt = rt(idx);

time_ard = cell(size(time_pc));
ballX_ard = cell(size(ballX_pc));

%%
mult_idx = [];
for i = 1:numel(idx)
    if choice(i) < 6 & rt(i) <= 2.1
        a = data.response.samples_start{idx(i)};
        b = data.response.samples_stop{idx(i)};
        ard_start = data.response.time(a+1);
        stop_ix = find(data.response.time >= ard_start+2,1,'first');
        
        time_ard{i} = data.response.time((a+1):stop_ix);
        ballX_ard{i} = cumsum(data.response.dx((a+1):stop_ix));
        dx_ard{i} = data.response.dx((a+1):stop_ix);
        try
            time_ard{i} = time_ard{i}-time_ard{i}(b-(a+1)); %zero arduino time
        catch
            time_ard{i} = time_ard{i}-time_ard{i}(end);
        end
        
        t_pc = time_pc{i}(find(samps_pc{i}~=0,1,'last'));
        %t_ard = time_ard{i}(b-a);
        time_ard{i} = time_ard{i}+t_pc;
        
        %clean up ard timing:
        %if first time point after 1s time point is >=10ms away, then x
        %position would be zero there
        index = find(time_ard{i} > 2,1,'first'); %means x did not change in the next 10ms, so dx = 0
        if ~isempty(index)
            if time_ard{i}(index) - 2 > 0.01
                time_ard{i}(index:end) = [];
                time_ard{i}(index) = 2;
                ballX_ard{i}(index:end) = [];
                ballX_ard{i}(index) = ballX_ard{i}(end);
            elseif time_ard{i}(index) - 2 <= 0.01 %means x changed in the next 10ms, pull that movement sample in
                time_ard{i}(index:end) = [];
                time_ard{i}(index) = 2;
                ballX_ard{i}(index+1:end) = [];
            end
        else
            time_ard{i}(end+1) = 2;
            ballX_ard{i}(end+1) = ballX_ard{i}(end);
        end
        
        %clean up multiple samples
        mult_idx{i} = find(diff(time_ard{i})==0);
        if ~isempty(mult_idx{i})
            trouble = [];
            for ii = 1:size(mult_idx{i},1)
                curr = mult_idx{i}(ii);
                curr = curr:curr+1;
                [~,index] = min(ballX_ard{i}(curr));
                trouble(ii) = curr(index);
            end
            ballX_ard{i}(trouble) = [];
            time_ard{i}(trouble) = [];
        end
        
        interp_ballX(i,:) = interp1(time_ard{i},ballX_ard{i},0.01:0.01:2);
        interp_ballX(i,isnan(interp_ballX(i,:))) = 0;
        ball_velocity(i,:) = diff(interp_ballX(i,:)*factor(1))./0.01;
    elseif choice(i) == 5 | rt(i) > 2.1
        interp_ballX(i,1:200) = nan;%interp1(time_ard{i},ballX_ard{i},0.01:0.01:2);
        ball_velocity(i,1:199) = nan;
    end
end

%%

% code trials: 1- correct left, nonlaser; 2- correct left, laser; 3:
% incorrect left, nonlaser; 4-incorrect left, laser; 5-incorrect right,
% nonlaser; 6-right incorrect, laser; 7-correct right; nonlaser;
% 8-correct right, laser
% trial_label = zeros(1,size(ball_velocity,1));
% k = 0;
% labels = [];
% for i = 1:2
%     for ii = 1:2
%         for iii = 1:2 %numel(unique(laser)) %velocity(loc)(choice)(laser)
%             k = k+1;
%             temp_idx = (loc == i & choice == ii & laser == iii);
%             trial_label(temp_idx) = k;
%             label_order = [labels; k i ii iii];
%         end
%     end
% end