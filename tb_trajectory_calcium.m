function [all_metric interp_ballX] = tb_trajectory_calcium(data);

plot_flag = 0;
plot_flag2 = 0;

% for f = 1:size(files,2);
%     disp(files{f});
%     side{f} = input('Left or Right: ');
% end
%
% save('side.mat','side','files');

[ball_velocity, interp_ballX, rt, ballX_pc, choice] = tb_get_velocity_calcium(data);


%%
pos_time = 0.01:0.01:2;
vel_time = pos_time(2:end);
accel_time = pos_time(3:end);
cons_pts = 3;%number of consecutive samples with increasing velocity
factor = data.response.mvmt_degrees_per_pixel;


for ii = 1:numel(rt)
    if choice(ii) == 2;
        scal = -1; %to account for right trials having negative velocities
    else
        scal = 1;
    end
    
    if choice(ii) < 5 & rt(ii)  <= 2.1
        curr_pos = interp_ballX(ii,:);
        curr_vel = ball_velocity(ii,:)*scal;
        curr_accel = diff(curr_vel)./accel_time;
        start = find(abs(curr_accel) > 100*factor(1));
        start2 = find(abs(curr_pos) > 5,1,'first');
        start(start < (start2-1)) = [];
        %         for iii = 1:numel(start)
        %             if iii < numel(start)-1 % if 3 consecutive criteria is not met, take the
        %                 first = find(diff(start(iii:iii+1)) == 1,1,'first');
        %                 if numel(first) == 1
        %                     mvmt_start{i}(ii) = accel_time(start(iii)) - 0.02;
        %                     break
        %                 end
        %             end
        %         end
        %
        %         start_idx = start(iii);
        start_idx = start(1);
        mvmt_start(ii) = accel_time(start(1)) - 0.02; % -0.02 to give mvmt_start in ball_position
        %time from movement start until reaching threshold
        mvmt_dur(ii) = rt(ii) - mvmt_start(ii);
        
        %peak velocity from movement start until reaching threshold
        curr_rt = round(rt(ii),2);%round rt
        rt_idx = fix((curr_rt*100)-1);% -1 becuase later you idx this to velocity
        if rt_idx > 199
            rt_idx = 199;
        end
        peak_vel(ii) = max(curr_vel(start_idx:rt_idx));
        
        
        mean_vel(ii) = mean(curr_vel(start:rt_idx)); %mean velocity in the correct direction
        %mean_vel{i}(ii) = mean(curr_vel(start_idx:rt_idx));
        vel_rt(ii) = curr_rt;
        mvmt_rt(ii) = rt(ii);
         plot_yno = 1;
    elseif choice(ii) == 5 | rt(ii) > 2.1
        mvmt_start(ii) = nan;
        mvmt_dur(ii) = nan;
        peak_vel(ii) = nan;
        mean_vel(ii) = nan;
        mvmt_rt(ii) = rt(ii);
        plot_yno = 0;
        
    end
    
    if plot_flag & plot_yno;
        %             subplot(2,1,1);
        %             plot(vel_time,curr_vel,'-k'); hold on
        %             plot(vel_time(start_idx),curr_vel(start_idx),'or');
        %             subplot(2,1,2);
        %if randi(10) < 6 & trial_label{i}(ii) == 1
        plot(pos_time,interp_ballX(ii,:)*factor(1),'-k'); hold on
        plot(mvmt_start(ii),curr_pos(start_idx)*factor,'or'); shg
        plot([0 0.8], [15 15],'--k','LineWidth',1);
%         ylim([-10 20]); xlim([0 0.8]);
        set(gca,'TickDir','out','TickLength',[0.02 0.05],'XTick',[0:0.1:0.8],'YTick',[-10:10:20]);
        box off
        pause
        clf
        %end
        
    end
    
    
end

%%

%% reorganize mean metric matrices based on contra/ipsi. First 4 columns are contra, last 4 are ipsi


all_metric(:,1) = mvmt_start;
all_metric(:,2) = mvmt_dur;
all_metric(:,3) = mean_vel;
all_metric(:,4) = peak_vel;
all_metric(:,5) = mvmt_rt;