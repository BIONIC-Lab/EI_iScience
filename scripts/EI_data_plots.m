%Code for plotting calcium data from excitatory/inhibitory data
load('C:\Users\chugh\Documents\GitHub\ExcitatoryInhibitory\stimdata.mat')
x_label = 1:1800;
x_label = x_label*0.0338; %converting to correct time domain

%% Get inhibitory and excitatory counts for each mouse
for mouse = 1:4
    exc_cnt(mouse) = sum(stimdata(mouse).data.emask);
    inh_cnt(mouse) = sum(stimdata(mouse).data.imask);
    ei_ratio(mouse) = exc_cnt(mouse)/inh_cnt(mouse);
end

%% First plot will just plot mean activity over time
stim_labels = {'10 Hz Uniform', 'Theta Bursts', '10 Hz Burst', '100 Hz'};
%cd('F:\ChronicEI\Figures\Mean time plots')
%for 10 Hz uniform
for stim_cnt = 1:4
    figure('PaperUnits', 'inches', 'PaperSize', [3.4, 6.8])
    for mouse = 1:4
        subplot(2,2,mouse)
        emask = logical(stimdata(mouse).data.emask);
        imask = logical(stimdata(mouse).data.imask);
        plot_data_exc = double(mean(stimdata(mouse).data.estim.dff(:,emask,stim_cnt),2));
        plot_data_inh = double(mean(stimdata(mouse).data.estim.dff(:,imask,stim_cnt),2));
        plot(x_label, plot_data_exc', 'LineWidth', 2, 'Color', 'r')
        hold on
        plot(x_label', plot_data_inh', 'LineWidth', 2, 'Color', 'b')
        hold off
        xlabel('Time (s)')
        ylabel('dF/F0')
        title(['Mean activity over time for ', animalids{mouse}])
        xline(10, 'r:')
        xline(40, 'b:')
        ylim([0 2])
    end
    sgtitle(stim_labels(stim_cnt))
    %print([stim_labels{stim_cnt}, '_mean'],'-dpng')
end

%% Can plot individual neurons to see if there 
stim_labels = {'10 Hz Uniform', 'Theta Bursts', '10 Hz Burst', '100 Hz'};
%for 10 Hz uniform
for stim_cnt = 1:4
    figure('PaperUnits', 'inches', 'PaperSize', [3.4, 6.8])
    for mouse = 1:4
        subplot(2,2,mouse)
        emask = logical(stimdata(mouse).data.emask);
        imask = logical(stimdata(mouse).data.imask);
        plot_data_exc = double(stimdata(mouse).data.estim.dff(:,emask,stim_cnt));
        plot_data_inh = double(stimdata(mouse).data.estim.dff(:,imask,stim_cnt));
        plot(x_label, plot_data_exc', 'LineWidth', 1, 'Color', 'r')
        hold on
        plot(x_label', plot_data_inh', 'LineWidth', 1, 'Color', 'b')
        hold off
        xlabel('Time (s)')
        ylabel('dF/F0')
        title(['All activity for ', animalids{mouse}])
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
    %print([animalids{mouse}, '_fullstim_norm'],'-dpng')
end

%% find the number of neurons that response most strongly to each stim type for each animal and cell type
stim_labels = {'10 Hz Uniform', 'Theta Bursts', '10 Hz Burst', '100 Hz'};
start_stim = find(x_label>10, 1, 'first');
end_stim = find(x_label>40, 1, 'first');
one_in = find(x_label>11, 1, 'first');
five_in = find(x_label>15, 1, 'first');
five_out = find(x_label>35, 1, 'first');
cd('F:\ChronicEI\Figures\Pie chart last 5s raw')
%for 10 Hz uniform
for mouse = 1:4
    figure('PaperUnits', 'inches', 'PaperSize', [3.4, 6.8])
    for stim_cnt = 1:4
        emask = logical(stimdata(mouse).data.emask);
        imask = logical(stimdata(mouse).data.imask);
        %divide freq_data by mouse and cell type
        freq_data{mouse,1}(:,stim_cnt) = mean(double(stimdata(mouse).data.estim.dff(start_stim:end_stim,emask,stim_cnt)),1);
        stimdata_sep{mouse,stim_cnt,1} = double(stimdata(mouse).data.estim.dff(:,emask,stim_cnt)); %to track activation of exc cells
        %freq_data{mouse,1}(:,stim_cnt) = freq_data{mouse,1}(:,stim_cnt)./max(freq_data{mouse,1}(:,stim_cnt));
        freq_data{mouse,2}(:,stim_cnt) = mean(double(stimdata(mouse).data.estim.dff(start_stim:end_stim,imask,stim_cnt)),1);
        stimdata_sep{mouse,stim_cnt,2} = double(stimdata(mouse).data.estim.dff(:,imask,stim_cnt)); %to track activation of exc cells
        %freq_data{mouse,2}(:,stim_cnt) = freq_data{mouse,2}(:,stim_cnt)./max(freq_data{mouse,2}(:,stim_cnt));
    end
    subplot(1,2,1)
    [~,pref_stim{mouse,1}] = max(freq_data{mouse,1},[],2);
    [gc, grps] = groupcounts(pref_stim{mouse,1});
    all_grps = zeros(1,4);
    all_grps(grps) = gc;
    p = pie(all_grps,stim_labels)
    title('Excitatory')
    subplot(1,2,2)
    [~,pref_stim{mouse,2}] = max(freq_data{mouse,2},[],2);
    [gc, grps] = groupcounts(pref_stim{mouse,2});
    all_grps = zeros(1,4);
    all_grps(grps) = gc;
    pie(all_grps,stim_labels)
    title('Inhibitory')
    sgtitle(animalids{mouse})
    %print([animalids{mouse}, '_pie_last5s'],'-dpng')
end

%% can also plot bar plots for all data
x_labels = {'10 Hz Uniform Exc.', '10 Hz Uniform Inh.', 'Theta Bursts Exc.', 'Theta Bursts Inh.', '10 Hz Bursts Exc.', '10 Hz Bursts Inh.', '100 Hz Exc.', '100 Hz Inh.'};
cd('F:\ChronicEI\Figures\Bar chart full stim raw')
for mouse = 1:4
    figure('PaperUnits', 'inches', 'PaperSize', [3.4, 6.8])
    x = [1,2,4,5,7,8,10,11];
    bar_data{mouse}(1,:) = mean(freq_data{mouse,1},1);
    bar_data{mouse}(2,:) = mean(freq_data{mouse,2},1);
    bar_data{mouse} = reshape(bar_data{mouse},1,8);
    bar_sem{mouse}(1,:) = std(freq_data{mouse,1},1)./sqrt(size(freq_data{mouse,1},1));
    bar_sem{mouse}(2,:) = std(freq_data{mouse,2},1)./sqrt(size(freq_data{mouse,2},1));
    bar_sem{mouse} = reshape(bar_sem{mouse},1,8);
    bar(x,bar_data{mouse})
    hold on
    e = errorbar(x, bar_data{mouse}, bar_sem{mouse}, 'k');
    e.LineStyle = 'none';
    ax = gca;
    ax.XTickLabels = x_labels;
    ylabel('Mean df/F0')
    title(animalids{mouse})
    %print([animalids{mouse}, '_bar_fullstim_raw'],'-dpng')
end

%% can take all neurons across all animals that respond most strongly to a specific stim type and plot the mean for each stim type
clear stim_data_pref
cd('F:\ChronicEI\Figures\Mean plot divided by preference')
for pref = 1:4
    figure('PaperUnits', 'inches', 'PaperSize', [3.4, 6.8])
    for stim_cnt = 1:4
        stim_data_pref{pref,stim_cnt,1} = stimdata_sep{1,1,1}(:,pref_stim{1,1} == pref);
        stim_data_pref{pref,stim_cnt,2} = stimdata_sep{1,1,2}(:,pref_stim{1,2} == pref);
        for mouse = 2:4
            stim_data_pref{pref,stim_cnt,1} = [stim_data_pref{pref,stim_cnt,1}, stimdata_sep{mouse,stim_cnt,1}(:,pref_stim{mouse,1} == pref)];
            stim_data_pref{pref,stim_cnt,2} = [stim_data_pref{pref,stim_cnt,2}, stimdata_sep{mouse,stim_cnt,2}(:,pref_stim{mouse,2} == pref)];
        end
        plot_data_1 = double(mean(stim_data_pref{pref,stim_cnt,1},2));
        plot_data_2 = double(mean(stim_data_pref{pref,stim_cnt,2},2));
        subplot(2,2,stim_cnt)
        plot(x_label, plot_data_1, 'LineWidth', 2, 'Color', 'r')
        hold on 
        plot(x_label, plot_data_2, 'LineWidth', 2, 'Color', 'b')
        xlabel('Time (s)')
        ylabel('dF/F0')
        title(stim_labels(stim_cnt))
        xline(10, 'r:')
        xline(40, 'b:')
        ylim([0 2])
    end
    sgtitle(['Preference for ', stim_labels(pref)]);
    %print([stim_labels{pref}, '_raw_preference'],'-dpng')
end
    
%% plot ei ratio vs preference (group down to two preferences) 0 = 100 Hz or burst, 1 = 10 Hz or Theta

for mouse = 1:4
    curr_stim_pref = [pref_stim{mouse,1}; pref_stim{mouse,2}];
    idx1 = curr_stim_pref == 3 | curr_stim_pref == 4;
    idx2 = curr_stim_pref == 1 | curr_stim_pref == 2;
    if sum(idx1) > sum(idx2)
        top_pref(mouse) = 0;
    else
        top_pref(mouse) = 1;
    end
end

figure('PaperUnits', 'inches', 'PaperSize', [3.4, 6.8])
plot(ei_ratio, top_pref, 'o')
    
