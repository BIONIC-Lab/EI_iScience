function plot_by_label(all_data,labels,label_types)

    %redefine basic variables
    frate = 0.034; %frame rate
    x_label = [1:1800]*frate; %time axis
    onset_start = 294; %approximate frame that ICMS started
    onset_stop = floor(onset_start+3/frate);
    offset_start = 1177;  %approximate frame ICMS ended
    offset_stop = floor(offset_start+3/frate);
    num_animals = 4; %total number of animals
    animals = {'EIM08', 'EIM09', 'VGAT2','VGAT5'}; %animal IDs
    train_types = {'10 Hz', '10 Hz Burst', '100 Hz', 'TBS'}; %different ICMS trains used 
    cell_types = {'Exc','Inh','Mix'};
    
    for curr_train = 1:4
        for curr_cell = 1:3
            legend_text{curr_cell} = [];
            figure 
            hold on
            num_rois = sum(~isnan(labels{curr_train}(curr_cell,:)));
            for curr_label = 1:5
                if sum(labels{curr_train}(curr_cell,:)==curr_label) >= ceil(0.02*num_rois)
                    plot(x_label, nanmean(all_data{curr_train,curr_cell}(labels{curr_train}(curr_cell,:)==curr_label,:),1), 'LineWidth', 2)
                    legend_text{curr_cell} = [legend_text{curr_cell}, label_types(curr_label)];
                end
            end
        %     plot(x_label, nanmean(all_data{curr_cell}(labels(curr_cell,:)==6,:),1))
        %     title(stim_types{curr_cell})
            title(cell_types{curr_cell})
            ylabel('dF/F0')
            xlabel('Time (s)')
            legend(legend_text{curr_cell})
            ylim([-0.4, 1.8])
        end

        % make pie plots for percentage that each train takes up
        for curr_cell = 1:3
            figure
            num_rois = sum(~isnan(labels{curr_train}(curr_cell,:)));
            for curr_label = 1:5
                curr_labs = sum(labels{curr_train}(curr_cell,:) == curr_label);
                percent_all(curr_label) = (curr_labs/num_rois)*100;
            end
            pie(percent_all)
            legend(legend_text{curr_cell})
        end
    end
end