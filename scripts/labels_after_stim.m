function [labels,offset_mean_all, prestim_mean_all] = labels_after_stim(all_data,error_all)

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
    offset_mean_all = [];
    prestim_mean_all = [];

    for curr_train = 1:4
        labels{curr_train} = NaN(3,1500);
        for curr_cell = 1:3
             for curr_neuron = 1:size(all_data{curr_train,curr_cell},1)
                 if all(~isnan(all_data{curr_train,curr_cell}(curr_neuron,:)))
                    %find max_value in first 3 s following ICMS
                    offset_end = max(all_data{curr_train,curr_cell}(curr_neuron,offset_start:offset_stop));
                    %find max_value after first 3 s following ICMS
                    offset_max = max(all_data{curr_train,curr_cell}(curr_neuron,offset_stop+1:end));
                    %find min value after first 3 s following ICMS
                    offset_min = min(all_data{curr_train,curr_cell}(curr_neuron,offset_stop+1:end));
                    offset_mean = nanmean(all_data{curr_train,curr_cell}(curr_neuron,offset_stop+1:end));
                    %baseline max value
                    %base_peak = max(all_data{curr_train,curr_cell}(curr_neuron,1:onset_start));
                    base_mean = nanmean(all_data{curr_train,curr_cell}(curr_neuron,1:onset_start));
                    %max for onset stim period
                    onset_end = max(all_data{curr_train,curr_cell}(curr_neuron,onset_stop:floor(onset_stop+10)));
                    %the standard deviation for the curr neuron
                    error = error_all{curr_train,curr_cell}(curr_neuron);
                    if offset_max > offset_end+error*2
                        labels{curr_train}(curr_cell,curr_neuron) = 4; %post-ICMS rebound
                    elseif offset_mean < -(abs(base_mean))-error*0.5
                        labels{curr_train}(curr_cell,curr_neuron) = 1; %post_ICMS depression
                    elseif offset_mean > (abs(base_mean))+error*0.5
                        labels{curr_train}(curr_cell, curr_neuron) = 3; %post-ICMS elevation
                    else
                        labels{curr_train}(curr_cell,curr_neuron) = 2; %decreases back to zero
                    end  
                    %save values for histogram plot
                    %divide by error to find general std
                    offset_mean_all(end+1) = offset_mean/error;
                    prestim_mean_all(end+1) = base_mean/error;
                 end
            end
        end
    end

end