frate = 0.034;

%find mean of exc and inh for 10 Hz Burst
% mean_exc = nanmean(all_data{2,1},1);
% mean_inh = nanmean(all_data{2,2},1);

bins = [0:25:540];
%bins = [0,225,540];
loc_e = NaN(100,21);
loc_i = NaN(100,21);
for curr_bin = 1:length(bins)-1
    clear temp temp_i 
    grouped_resps_e = [];
    grouped_resps_i = [];
    %plot average waveform based on distance
    dist_idx_e = active_neur_dist{1,2} > bins(curr_bin) & active_neur_dist{1,2} <= bins(curr_bin+1);
    %dist_idx_e = active_neur_dist{1,2} > 200;
    e_vals = active_neur_vals{1,2}(dist_idx_e,:);
    %mean_exc = nanmean(e_vals,1);
    
    dist_idx_i = active_neur_dist{2,2} > bins(curr_bin) & active_neur_dist{2,2} <= bins(curr_bin+1);
    %dist_idx_i = active_neur_dist{2,2} > 200;
    i_vals = active_neur_vals{2,2}(dist_idx_i,:);
    %mean_inh = nanmean(i_vals,1);

    %group responses of last 6 bursts into one average
    for curr_neuron = 1:size(e_vals,1)
        for i = 1:6
            temp(i,:) = e_vals(curr_neuron,996+(30*(i-1)):996+(30*i));
        end
        grouped_resps_e(curr_neuron,:) = nanmean(temp,1);
    end
    
    for curr_neuron = 1:size(i_vals,1)
        for i = 1:6
            temp_i(i,:) = i_vals(curr_neuron,996+(30*(i-1)):996+(30*i));
        end
        grouped_resps_i(curr_neuron,:) = nanmean(temp_i,1);
    end


%     for i = 1:6
%         grouped_resps_e(i,:) = mean_exc(994+(31*(i-1)):994+(31*i));
%     end

    %take mean across 6 bursts
    if ~isempty(grouped_resps_e)
        mean_grouped_e = nanmean(grouped_resps_e,1);
        [~,loc_e(1:size(grouped_resps_e,1),curr_bin)] = max(grouped_resps_e(:,1:15)'); %limit it to reasonable range
    else
        mean_grouped_e = NaN(1,31);
    end
    
    if ~isempty(grouped_resps_i)
        mean_grouped_i = nanmean(grouped_resps_i,1);
        [~,loc_i(1:size(grouped_resps_i,1),curr_bin)] = max(grouped_resps_i(:,1:15)'); %limit it to reasonable range
    else
        mean_grouped_i = NaN(1,31);
    end
    %[~,loc_i(curr_bin)] = max(mean_grouped_i(:,1:15));
    
    figure
    plot([1:31]*frate,mean_grouped_e./max(mean_grouped_e), 'LineWidth', 2)
    hold on
    plot([1:31]*frate, mean_grouped_i./max(mean_grouped_i), 'LineWidth', 2)
    ylim([0 1.1])
    ylim([0.6 1.1])
    legend('Exc', 'Inh')
    %title({'10HzB last 6 bursts avg for exc vs. inh at ', num2str(bins(curr_bin+1)), ' um'})
    ylabel('Norm dF/F0')
    xlabel('Time (s)')
end
return

figure 
mean_e = nanmean(loc_e(:,5:11)*frate,1);
plot(mean_e, 'LineWidth',2)
hold on
sem_e = nanstd(loc_e(:,5:11)*frate)./sqrt(sum(~isnan(loc_e(:,5:11))))
errorbar([1:7], mean_e, sem_e,'.')
mean_i = nanmean(loc_i(:,5:11)*frate,1);
plot(mean_i, 'LineWidth',2)
sem_i = nanstd(loc_i(:,5:11)*frate)./sqrt(sum(~isnan(loc_i(:,5:11))))
errorbar([1:7], mean_i, sem_i,'.')
xlim([0,8])
ax = gca;
ax.XTick = [1:7];
ax.XTickLabel = {'100-125','125-150','150-175','175-200','200-225','225-250','250-275'};
xlabel('Dist bins')
ylabel('Time to peak (s)')

%find that no individual bin is statistically significant with repeated
%measures, but can combine <225 and >225
temp = loc_i(:,5:9)*frate;
loc_i_225less = temp(~isnan(temp));
loc_i_225less = [loc_i_225less; NaN(1500-length(loc_i_225less),1)];

temp = loc_e(:,5:9)*frate;
loc_e_225less = temp(~isnan(temp));
loc_e_225less = [loc_e_225less; NaN(1500-length(loc_e_225less),1)];

%bar plot
figure
bar_data = [nanmean(loc_e_225less), nanmean(loc_i_225less)];
bar_error = [nanstd(loc_e_225less)./sqrt(sum(~isnan(loc_e_225less))), nanstd(loc_i_225less)./sqrt(sum(~isnan(loc_i_225less)))];
bar(bar_data)
hold on
errorbar([1,2], bar_data, bar_error, '.');
ax = gca;
ax.XTick = [1:2];
ax.XTickLabels = {'Exc','Inh'};
ylabel('Time to peak (s)')
title('<225 \mum from electrode')

temp = loc_i(:,10:11)*frate;
loc_i_225 = temp(~isnan(temp));
loc_i_225 = [loc_i_225; NaN(1500-length(loc_i_225),1)];

temp = loc_e(:,10:11)*frate;
loc_e_225 = temp(~isnan(temp));
loc_e_225 = [loc_e_225; NaN(1500-length(loc_e_225),1)];

%bar plot
figure
bar_data = [nanmean(loc_e_225), nanmean(loc_i_225)];
bar_error = [nanstd(loc_e_225)./sqrt(sum(~isnan(loc_e_225))), nanstd(loc_i_225)./sqrt(sum(~isnan(loc_i_225)))];
bar(bar_data)
hold on
errorbar([1,2], bar_data, bar_error, '.');
ax = gca;
ax.XTick = [1:2];
ax.XTickLabels = {'Exc','Inh'};
ylabel('Time to peak (s)')
title('>225 \mum from electrode')


