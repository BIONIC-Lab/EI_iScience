%Code for plotting calcium data from excitatory/inhibitory data -VGAT5
%multiday
load('F:\ChronicEI\stimdata_VGAT5_multiday.mat')
x_label = 1:1800;
x_label = x_label*0.0338; %converting to correct time domain
all_days = {stimdata.data.day};

%% Get inhibitory and excitatory counts for each mouse
for curr_day = 1:6
    exc_cnt(curr_day) = sum(stimdata.data(curr_day).emask);
    inh_cnt(curr_day) = sum(stimdata.data(curr_day).imask);
    ei_ratio(curr_day) = exc_cnt(curr_day)/inh_cnt(curr_day);
end

%% First plot will just plot mean activity over time
stim_labels = {'10 Hz Uniform', 'Theta Bursts', '10 Hz Burst', '100 Hz'};
cd('F:\ChronicEI\Figures\Mean time plots VGAT5')
%for 10 Hz uniform
for stim_cnt = 1:4
    figure('PaperUnits', 'inches', 'PaperSize', [3.4, 6.8])
    for curr_day = 1:(length(stimdata.data))
        subplot(2,3,curr_day)
        emask = logical(stimdata.data(curr_day).emask);
        imask = logical(stimdata.data(curr_day).imask);
        plot_data_exc = double(mean(stimdata.data(curr_day).estim.dff(:,emask,stim_cnt),2));
        plot_data_inh = double(mean(stimdata.data(curr_day).estim.dff(:,imask,stim_cnt),2));
        plot(x_label, plot_data_exc', 'LineWidth', 2, 'Color', 'r')
        hold on
        plot(x_label', plot_data_inh', 'LineWidth', 2, 'Color', 'b')
        hold off
        xlabel('Time (s)')
        ylabel('dF/F0')
        title(['Mean activity ', stimdata.data(curr_day).day])
        xline(10, 'r:')
        xline(40, 'b:')
        ylim([0 2])
    end
    sgtitle(stim_labels(stim_cnt))
    print([stim_labels{stim_cnt}, '_mean_VGAT5'],'-dpng')
end

%% Can plot individual neurons to see if there 
stim_labels = {'10 Hz Uniform', 'Theta Bursts', '10 Hz Burst', '100 Hz'};
%for 10 Hz uniform
for stim_cnt = 1:4
    figure('PaperUnits', 'inches', 'PaperSize', [3.4, 6.8])
    for curr_day = 1:6
        subplot(2,3,curr_day)
        emask = logical(stimdata.data(curr_day).emask);
        imask = logical(stimdata.data(curr_day).imask);
        plot_data_exc = double(stimdata.data(curr_day).estim.dff(:,emask,stim_cnt));
        plot_data_inh = double(stimdata.data(curr_day).estim.dff(:,imask,stim_cnt));
        plot(x_label, plot_data_exc', 'LineWidth', 1, 'Color', 'r')
        hold on
        plot(x_label', plot_data_inh', 'LineWidth', 1, 'Color', 'b')
        hold off
        xlabel('Time (s)')
        ylabel('dF/F0')
        title(['All activity for ', all_days{curr_day}])
        xline(10, 'r:')
        xline(40, 'b:')
        %ylim([0 2])
    end
    sgtitle(stim_labels(stim_cnt))
end

%% find activity for individual neurons in the first 5-s of stim, the last 5-s of stim, the whole window
stim_labels = {'10 Hz Uniform', '10 Hz Burst', '100 Hz'};
start_stim = find(x_label>10, 1, 'first');
end_stim = find(x_label>40, 1, 'first');
five_in = find(x_label>15, 1, 'first');
five_out = find(x_label>35, 1, 'first');
cd('F:\ChronicEI\Figures\Full stim period norm')
%for 10 Hz uniform
for mouse = 1:4
    figure('PaperUnits', 'inches', 'PaperSize', [3.4, 6.8])
    for stim_cnt = [1,3:4]
        emask = logical(stimdata(mouse).data.emask);
        imask = logical(stimdata(mouse).data.imask);
        %divide freq_data by mouse and cell type
        freq_data{mouse,1}(:,stim_cnt) = mean(double(stimdata(mouse).data.estim.dff(start_stim:end_stim,emask,stim_cnt)),1);
        freq_data{mouse,1}(:,stim_cnt) = freq_data{mouse,1}(:,stim_cnt)./max(freq_data{mouse,1}(:,stim_cnt));
        freq_data{mouse,2}(:,stim_cnt) = mean(double(stimdata(mouse).data.estim.dff(start_stim:end_stim,imask,stim_cnt)),1);
        freq_data{mouse,2}(:,stim_cnt) = freq_data{mouse,2}(:,stim_cnt)./max(freq_data{mouse,2}(:,stim_cnt));
    end
    plot3(freq_data{mouse,1}(:,1),freq_data{mouse,1}(:,2),freq_data{mouse,1}(:,3), 'r.', 'MarkerSize', 12)
    hold on
    plot3(freq_data{mouse,2}(:,1),freq_data{mouse,2}(:,2),freq_data{mouse,2}(:,3), 'b.', 'MarkerSize', 12)
    xlabel(stim_labels{1})
    ylabel(stim_labels{2})
    zlabel(stim_labels{3})
    title(animalids{mouse})
    print([animalids{mouse}, '_fullstim_norm'],'-dpng')
end

%% find the number of neurons that response most strongly to each stim type for each animal and cell type
stim_labels = {'10 Hz Uniform', 'Theta Bursts', '10 Hz Burst', '100 Hz'};
start_stim = find(x_label>10, 1, 'first');
end_stim = find(x_label>40, 1, 'first');
five_in = find(x_label>15, 1, 'first');
five_out = find(x_label>35, 1, 'first');
cd('F:\ChronicEI\Figures\Pie chart full stim raw VGAT5')
%for 10 Hz uniform
for curr_day = 1:6
    figure('PaperUnits', 'inches', 'PaperSize', [3.4, 6.8])
    for stim_cnt = 1:4
        emask = logical(stimdata.data(curr_day).emask);
        imask = logical(stimdata.data(curr_day).imask);
        %divide freq_data by mouse and cell type
        freq_data{curr_day,1}(:,stim_cnt) = mean(double(stimdata.data(curr_day).estim.dff(start_stim:end_stim,emask,stim_cnt)),1);
        %freq_data{mouse,1}(:,stim_cnt) = freq_data{mouse,1}(:,stim_cnt)./max(freq_data{mouse,1}(:,stim_cnt));
        freq_data{curr_day,2}(:,stim_cnt) = mean(double(stimdata.data(curr_day).estim.dff(start_stim:end_stim,imask,stim_cnt)),1);
        %freq_data{mouse,2}(:,stim_cnt) = freq_data{mouse,2}(:,stim_cnt)./max(freq_data{mouse,2}(:,stim_cnt));
    end
    subplot(1,2,1)
    [~,pref_stim] = max(freq_data{curr_day,1},[],2);
    [gc, grps] = groupcounts(pref_stim);
    all_grps = zeros(1,4);
    all_grps(grps) = gc;
    p = pie(all_grps,stim_labels)
    title('Excitatory')
    subplot(1,2,2)
    [~,pref_stim] = max(freq_data{curr_day,2},[],2);
    [gc, grps] = groupcounts(pref_stim);
    all_grps = zeros(1,4);
    all_grps(grps) = gc;
    pie(all_grps,stim_labels)
    title('Inhibitory')
    sgtitle(stimdata.data(curr_day).day)
    print([stimdata.data(curr_day).day, '_pie'],'-dpng')
end

%% can also plot bar plots for all data
x_labels = {'10 Hz Uniform Exc.', '10 Hz Uniform Inh.', 'Theta Bursts Exc.', 'Theta Bursts Inh.', '10 Hz Bursts Exc.', '10 Hz Bursts Inh.', '100 Hz Exc.', '100 Hz Inh.'};
cd('F:\ChronicEI\Figures\Bar chart full stim raw VGAT5')
for curr_day = 1:6
    figure('PaperUnits', 'inches', 'PaperSize', [3.4, 6.8])
    x = [1,2,4,5,7,8,10,11];
    bar_data{curr_day}(1,:) = mean(freq_data{curr_day,1},1);
    bar_data{curr_day}(2,:) = mean(freq_data{curr_day,2},1);
    bar_data{curr_day} = reshape(bar_data{curr_day},1,8);
    bar_sem{curr_day}(1,:) = std(freq_data{curr_day,1},1)./sqrt(size(freq_data{curr_day,1},1));
    bar_sem{curr_day}(2,:) = std(freq_data{curr_day,2},1)./sqrt(size(freq_data{curr_day,2},1));
    bar_sem{curr_day} = reshape(bar_sem{curr_day},1,8);
    bar(x,bar_data{curr_day})
    hold on
    e = errorbar(x, bar_data{curr_day}, bar_sem{curr_day}, 'k');
    e.LineStyle = 'none';
    ax = gca;
    ax.XTickLabels = x_labels;
    ylabel('Mean df/F0')
    title(stimdata.data(curr_day).day)
    print([stimdata.data(curr_day).day, '_bar_fullstim_raw_VGAT5'],'-dpng')
end