function [labels] = labels_during_stim_split(all_data_split,error_all_split)

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
    
    for curr_cell = 1:3
        for curr_train = 1:4
            for curr_animal = 1:num_animals
                for curr_neuron = 1:size(all_data_split{curr_train,curr_animal,curr_cell},1)
                     if all(~isnan(all_data_split{curr_train,curr_animal,curr_cell}(curr_neuron,:)))
                        first_peak = max(all_data_split{curr_train,curr_animal,curr_cell}(curr_neuron,onset_start:onset_stop));
                        second_peak = max(all_data_split{curr_train,curr_animal,curr_cell}(curr_neuron,offset_start:offset_stop));
                        hold_peak = max(all_data_split{curr_train,curr_animal,curr_cell}(curr_neuron,onset_stop:offset_start));
                        base_peak = max(all_data_split{curr_train,curr_animal,curr_cell}(curr_neuron,1:onset_start));
                        base_mean = nanmean(all_data_split{curr_train,curr_animal,curr_cell}(curr_neuron,1:onset_start));
                        onset_end = max(all_data_split{curr_train,curr_animal,curr_cell}(curr_neuron,onset_stop:onset_stop+10));
                        error = error_all_split{curr_train,curr_animal,curr_cell}(curr_neuron);
                        if first_peak > hold_peak+(error/2) 
                            if onset_end < (first_peak/2) %separating slowly and rapid adapting
                                labels{curr_cell,curr_animal}(curr_train,curr_neuron) = 1; %rapid adapt
                            else
                                labels{curr_cell,curr_animal}(curr_train,curr_neuron) = 2; %slowly adapt
                            end
                        elseif hold_peak > first_peak+(error*3) 
                            labels{curr_cell,curr_animal}(curr_train,curr_neuron) = 4; %facilitating
                        else
                            labels{curr_cell,curr_animal}(curr_train,curr_neuron) = 3; %steady state (could be further divided by facilitating)
                        end  
                        if nanmean(all_data_split{curr_train,curr_animal,curr_cell}(curr_neuron,onset_start:offset_start)) < (-(abs(base_mean))-(error*0.5))
                            labels{curr_cell,curr_animal}(curr_train,curr_neuron) = 5; %negative modulation
                        end
                     else
                          labels{curr_cell,curr_animal}(curr_train,curr_neuron) = 0;
                    end
                end
         %label_names{curr_cell,curr_train}(~label_log) = [];
            end
        end
    end

%might need something to save which labels that are active for each

%above should now be split by train but need to fix other code

end