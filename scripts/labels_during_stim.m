%function for labels during stim
function [labels, max_peak_all, adapt_all] = labels_during_stim(all_data,error_all)

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
    max_peak_all = [];
    adapt_all = [];

    for curr_train = 1:4
        labels{curr_train} = NaN(3,1500);
        for curr_cell = 1:3
             for curr_neuron = 1:size(all_data{curr_train,curr_cell},1)
                 if all(~isnan(all_data{curr_train,curr_cell}(curr_neuron,:)))
                    first_peak = max(all_data{curr_train,curr_cell}(curr_neuron,onset_start:onset_stop));
                    second_peak = max(all_data{curr_train,curr_cell}(curr_neuron,offset_start:offset_stop));
                    max_peak = max(all_data{curr_train,curr_cell}(curr_neuron,onset_start:offset_start));
                    hold_peak = max(all_data{curr_train,curr_cell}(curr_neuron,onset_stop:offset_start));
                    base_peak = max(all_data{curr_train,curr_cell}(curr_neuron,1:onset_start));
                    base_mean = nanmean(all_data{curr_train,curr_cell}(curr_neuron,1:onset_start));
                    onset_end = max(all_data{curr_train,curr_cell}(curr_neuron,onset_stop:onset_stop+10));
                    error = error_all{curr_train,curr_cell}(curr_neuron);
                    if first_peak > hold_peak+(error/2) 
                        if onset_end < (first_peak/2) %separating slowly and rapid adapting
                            labels{curr_train}(curr_cell,curr_neuron) = 1; %rapid adapt
                        else
                            labels{curr_train}(curr_cell,curr_neuron) = 2; %slowly adapt
                        end
                    elseif hold_peak > first_peak+(error*3) 
                        labels{curr_train}(curr_cell,curr_neuron) = 4; %facilitating
                    else
                        labels{curr_train}(curr_cell,curr_neuron) = 3; %steady state (could be further divided by facilitating)
                    end
                    if nanmean(all_data{curr_train,curr_cell}(curr_neuron,onset_start:offset_start)) < (-(abs(base_mean))-(error*0.5))
                        labels{curr_train}(curr_cell,curr_neuron) = 5; %negative modulation
                    end   
                    %for histogram of active neurons
                    max_peak_all(end+1) = max_peak/error;
                    %histogram for slowly adapting and facilitating
                    adapt_all(end+1) = (first_peak-hold_peak)/error;
                 end
            end
        end        
    end

end