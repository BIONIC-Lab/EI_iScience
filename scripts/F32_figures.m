%Code for generating figures for F32
load('C:\Users\chugh\OneDrive\Documents\GitHub\ExcitatoryInhibitory\stimdata.mat')
x_label = 1:1800;
x_label = x_label*0.0338; %converting to correct time domain

%% First plot will just plot mean activity over time
stim_labels = {'10 Hz',[],[],'100 Hz'};
%cd('F:\ChronicEI\Figures\Mean time plots')
%for 10 Hz uniform
for stim_cnt = [1,4]
    figure('PaperUnits', 'inches', 'PaperSize', [3.4, 6.8])
    for mouse = 2 %only EIF02
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
    %sgtitle(stim_labels(stim_cnt))
    %print([stim_labels{stim_cnt}, '_mean'],'-dpng')
end

%% pca analysis
data_100_e = stimdata(2).data.estim.dff(:,emask,4);
data_100_i = stimdata(2).data.estim.dff(:,imask,4);
data_10_e = stimdata(2).data.estim.dff(:,emask,1);
data_10_i = stimdata(2).data.estim.dff(:,imask,1);

plot_mean_activity(data_100_e, 3)

function plot_mean_activity(data, ks)
    [c,s,l] = pca(data);
    
    %remove outliers 
%     sum_idx = sum(c(:,1:3),2);
%     c = c(sum_idx<0.5,:);
    
    %k-means cluster on pca results
    [idx2, cd, sumd] = kmeans(c(:,1:3),ks);

    %plot results
    figure
    for curr_fac = 1:ks
        curr_idx = idx2 == curr_fac;
        plot3(c(curr_idx,1), c(curr_idx ,2), c(curr_idx ,3), '.')
        hold on
    end

    figure
    for curr_fac = 1:ks
        curr_idx = idx2 == curr_fac;
        plot(mean(data(:,curr_idx),2))
        hold on
    end

end