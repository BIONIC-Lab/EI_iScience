%code to calculate correlation indices between active inhibitory neurons
%and all other neurons

%will look at all trains but for now let's start with 100 Hz

%initialize correlation and lag arrays
%can initialize to 3000x3000 NaN (more than enough neurons in both
%directions) and then ignore NaNs later for analysis
onset_start = 294; %approximate frame that ICMS started
onset_stop = floor(onset_start+3/frate);
offset_start = 1177;  %approximate frame ICMS ended
offset_stop = floor(offset_start+3/frate);

clear curr_resp comp_resp
for curr_animal = 1:4
    for curr_cell = 1:3
        r_max{curr_animal,curr_cell} = NaN(3000,3000);
        r_max_base{curr_animal,curr_cell} = NaN(3000,3000);
        r_max_wn{curr_animal,curr_cell} = NaN(3000,3000);
        r_max_base_wn{curr_animal,curr_cell} = NaN(3000,3000);
        opt_lag{curr_animal,curr_cell} = NaN(3000,3000);
        diff_dist{curr_animal,curr_cell} = NaN(3000,3000);
        comp_label{curr_animal,curr_cell} = NaN(3000,3000);
    end
end

curr_train = 3; %only looking at one train right now
base_cell = 2; %only looking at one cell type to compare to the others
for curr_animal = 1:4
    for curr_neuron = 1:size(all_data_split{curr_train,curr_animal,base_cell},1) %pull data of neuron at current animal and selected train+base cell
        if ~isnan(all_data_split{curr_train,curr_animal,base_cell}(curr_neuron,1)) %check if the current neuron is active
            curr_resp{curr_animal,curr_neuron} = all_data_split{curr_train,curr_animal,base_cell}(curr_neuron,:); %extract resp of individual neuron
            base_label(curr_animal,curr_neuron) = labels_split{base_cell,curr_animal}(curr_train,curr_neuron);
            %curr_resp{curr_animal,curr_neuron} = curr_resp{curr_animal,curr_neuron}./max(curr_resp{curr_animal,curr_neuron}); %test for normalization
%             curr_dist = active_neur_dist_split{base_cell,curr_animal}(curr_train,curr_neuron);
            curr_dist = all_rois{curr_animal}(curr_neuron,:); %we make all arrays same size across cell types so no need to specify that 
            %slightly worried that individual neurons will be too noisy but
            %we will see
            for curr_cell = 1:3 %compare to each cell type
                for curr_comp_neuron = 1:size(all_data_split{curr_train,curr_animal,curr_cell},1)
                    if ~isnan(all_data_split{curr_train,curr_animal,curr_cell}(curr_comp_neuron,1)) %check if there is an active comparison neuron
                        %corr with curr inhibitory neuron and comparison
                        %neuron
                        %lose information here because not specifying
                        %curr_cell
                        comp_resp{curr_animal,curr_neuron,curr_cell}(curr_comp_neuron,:) = all_data_split{curr_train,curr_animal,curr_cell}(curr_comp_neuron,:);
                        comp_label{curr_animal,curr_cell}(curr_neuron,curr_comp_neuron) = labels_split{curr_cell,curr_animal}(curr_train,curr_comp_neuron);
                        %comp_resp{curr_animal,curr_neuron}(curr_comp_neuron,:) = comp_resp{curr_animal,curr_neuron}(curr_comp_neuron,:)./max(comp_resp{curr_animal,curr_neuron}(curr_comp_neuron,:)); %test for normalization
                                               
                        %using cross-correlation
                        %[r,lag] = xcorr(curr_resp{curr_animal,curr_neuron},comp_resp{curr_animal,curr_neuron,curr_cell}(curr_comp_neuron,:), 'normalized');
                        %[r_max{curr_animal,curr_cell}(curr_neuron,curr_comp_neuron)] = r(1800); %no lag %max(r);
                        %opt_lag{curr_animal,curr_cell}(curr_neuron,curr_comp_neuron) = lag(loc);
                        
                        %using Pearson's correlation
                        %r1 = [curr_resp{curr_animal,curr_neuron}(1:onset_start-1),curr_resp{curr_animal,curr_neuron}(offset_start:end)];
                        r1 = curr_resp{curr_animal,curr_neuron};
                        %r2 = [comp_resp{curr_animal,curr_neuron,curr_cell}(curr_comp_neuron,1:onset_start-1),comp_resp{curr_animal,curr_neuron,curr_cell}(curr_comp_neuron,offset_start:end)];
                        r2 = comp_resp{curr_animal,curr_neuron,curr_cell}(curr_comp_neuron,:);
                        [r,~] = corrcoef(r1,r2);
                        r_max{curr_animal,curr_cell}(curr_neuron,curr_comp_neuron) = r(1,2);
                        
                        %can also look at baseline correlation
                        [r_base,~] = corrcoef(curr_resp{curr_animal,curr_neuron}(1:onset_start-1),comp_resp{curr_animal,curr_neuron,curr_cell}(curr_comp_neuron,1:onset_start-1));
                        r_max_base{curr_animal,curr_cell}(curr_neuron,curr_comp_neuron) = r_base(1,2);
                        
                        %also noise correlation - for noise correlation
                        %they typically use repeated trials and subtract
                        %the average to find the noise around the
                        %distribution -we don't have this 
                        %we can potentially calculate this kind of thing by
                        %subtracting the baseline by the baseline of the
                        %mean - can try anyway 
                        
                        %can also correlate to white noise for comparison
                        curr_white_noise = smoothdata(wgn(length(curr_resp{curr_animal,curr_neuron}),1,0),'Gaussian',[7,7]);
                        [r_wn,~] = corrcoef(curr_resp{curr_animal,curr_neuron},curr_white_noise);
                        r_max_wn{curr_animal,curr_cell}(curr_neuron,curr_comp_neuron) = r_wn(1,2);
                        
                        curr_white_noise = smoothdata(wgn(length(curr_resp{curr_animal,curr_neuron}(1:onset_start-1)),1,0),'Gaussian',[7,7]);
                        [r_base_wn,~] = corrcoef(curr_resp{curr_animal,curr_neuron}(1:onset_start-1),curr_white_noise);
                        r_max_base_wn{curr_animal,curr_cell}(curr_neuron,curr_comp_neuron) = r_base_wn(1,2);
                        
                        comp_dist = all_rois{curr_animal}(curr_comp_neuron,:);
                        diff_dist{curr_animal,curr_cell}(curr_neuron,curr_comp_neuron) = sqrt((curr_dist(:,2) - comp_dist(:,2)).^2+(curr_dist(:,3) - comp_dist(:,3)).^2);
                        %can also potentially do some correlation of
                        %excitatory neurons with each other and so on
                        comp_dist = all_rois{curr_animal}(curr_comp_neuron,:);
                        diff_dist{curr_animal,curr_cell}(curr_neuron,curr_comp_neuron) = sqrt((curr_dist(:,2) - comp_dist(:,2)).^2+(curr_dist(:,3) - comp_dist(:,3)).^2);
                        %can also potentially do some correlation of
                        %excitatory neurons with each other and so on
                        
                    end
                end
            end
%             r_all = [r_max_vs; r_max_vb; r_max_sb];
%             figure
%             bar(r_all')
%             title('Max correlation')
%             ax = gca;
%             ax.XTickLabels = stim_types;
%             legend('Male vs. Female', 'Male vs. Baseline', 'Female vs. Baseline')
            end
    end
end

%% have to now figure out how to implement code to look at this as a function
%of distance of neurons from each other, can maybe simultaneously calculate
%distance between neurons when calculating correlation

%distance arrays currently only include information about distance from the
%electrode, but we want to know the distance of the neurons from each other
%so will need to reload distance info in for these calculations of distance
%between neurons

%can make a plot showing correlation indices between excitatory and
%inhibitory neurons
%calculate average correlation indices for inhibitory neurons vs exc, inh,
%and mix

%combine correlations across animals (so only divided by cell type)
for curr_cell = 1:3
    temp = [];
    for curr_animal = 1:4
        temp = [temp; nanmean(r_max{curr_animal,curr_cell},2)]; %taking mean across comparisons
    end
    corr_cell(curr_cell,:) = temp;
end
%this was just for potentially generating a plot that shows the average
%correlation of each cell type to all other cells (so we can see if
%selected cell type is more correlated to inhibitory or excitatory neurons)
figure
violin(corr_cell')

%I think unsurprisingly, we find inhibitory neurons are more correlated to
%other inhibitory neurons

%% correlation as a function of distance
%this seems like it might work, but I am realizing I might want to divide
%neurons close to and far from the electrode
%can use just inhibitory neuron to divide in this way
%may want to save the average waveforms for comparison?
clear curr_comp_resps
bins = [0,25,50,100,200,540];
big_bins = [0,200,540];
for curr_big_bin = 1:2
    for curr_cell = 1:3
        corr_dist{curr_cell,curr_big_bin} = NaN(5,3000);
        for curr_bin = 1:5
            temp = [];
            temp2 = [];
            temp3 = [];
            temp4 = [];
            for curr_animal = 1:4
                %curr_comp_resps{curr_neuron,curr_animal,curr_bin,curr_cell} = NaN(1,1800);
                big_bin_idx = active_neur_dist_split{base_cell,curr_train}(:,curr_animal) > big_bins(curr_big_bin) & active_neur_dist_split{base_cell,curr_train}(:,curr_animal) < big_bins(curr_big_bin+1);
                %this will put zeros where there are NaNs but that should
                %be okay
                temp_diff_array = diff_dist{curr_animal,curr_cell};
                temp_diff_array(~big_bin_idx,:) = NaN; %remove ones not in current big bin (<200 or >200)
                %above works on dimension for base neuron, which should be
                %correct (orienting to base neuron)
                dist_idx = temp_diff_array > bins(curr_bin) & temp_diff_array <= bins(curr_bin+1) & temp_diff_array~=0; %CHECK THIS AGAIN
                temp = [temp; r_max{curr_animal,curr_cell}(dist_idx)]; %not specifying dimension here because dist_idx is 2D and same size 3000x3000
                temp2 = [temp2; r_max_base{curr_animal,curr_cell}(dist_idx)];
                temp3 = [temp3; r_max_wn{curr_animal,curr_cell}(dist_idx)];
                temp4 = [temp4; r_max_base_wn{curr_animal,curr_cell}(dist_idx)];
                
                %pull comp responses?
                for curr_neuron = 1:size(comp_resp,2)
                   resp_idx = diff_dist{curr_animal,curr_cell}(curr_neuron,:) > bins(curr_bin) & diff_dist{curr_animal,curr_cell}(curr_neuron,:) <= bins(curr_bin+1) & diff_dist{curr_animal,curr_cell}(curr_neuron,:)~=0;
                   curr_comp_resps{curr_neuron,curr_animal,curr_bin,curr_cell} = [comp_resp{curr_animal,curr_neuron,curr_cell}(resp_idx,:)];
                end
            end
            %first pull previous info from corr_dist
            corr_dist{curr_cell,curr_big_bin}(curr_bin,:) = [temp; NaN(3000-length(temp),1)]; 
            corr_dist_base{curr_cell,curr_big_bin}(curr_bin,:) = [temp2; NaN(3000-length(temp2),1)]; 
            corr_dist_wn{curr_cell,curr_big_bin}(curr_bin,:) = [temp3; NaN(3000-length(temp3),1)]; 
            corr_dist_base_wn{curr_cell,curr_big_bin}(curr_bin,:) = [temp4; NaN(3000-length(temp4),1)]; 
        end
    end
end

%reason that inhibitory neurons have higher xcorr may be because they are
%more facilitating, and facilitating tends to have higher xcorr values with
%all traces, can try to check on this

%can also look at average traces for each subtype in each bin

%% can plot resps against comp resps in different bins to look at how correlation might work
curr_animal = 3; %desired animal
curr_big_bins_int = active_neur_dist_split{base_cell,curr_train}(:,curr_animal); %can look through this to find cells with desired distance bin
curr_neuron = 55; %choose a desired neuron after reviewing curr_big_bins_int
figure, plot(x_label,curr_resp{curr_animal,curr_neuron},'LineWidth',2) %base cell for comparison
hold on
first_bin = 1; %choose the closest bin you want, may be limited by how many active neurons are close to neuron
second_bin = 5;
curr_cell = 2; %choose cell type you want to compare to
plot(x_label,curr_comp_resps{curr_neuron,curr_animal,first_bin,curr_cell},'LineWidth',2)
plot(x_label,nanmean(curr_comp_resps{curr_neuron,curr_animal,second_bin,curr_cell}),'LineWidth',2)
title('Example exc-inh greater than 200 um') %this will depend on comparison being made
legend({'Base','0-25 um', '200+'}) %this will depend on bins used

%can potentially add code to extract correlation values and add to legend?
curr_bin = 5;
resp_idx = diff_dist{curr_animal,curr_cell}(curr_neuron,:) > bins(curr_bin) & diff_dist{curr_animal,curr_cell}(curr_neuron,:) <= bins(curr_bin+1) & diff_dist{curr_animal,curr_cell}(curr_neuron,:)~=0;
curr_mean = nanmean(r_max{curr_animal,curr_cell}(curr_neuron,resp_idx));
%above seems to work
%% potential plotting
cell_types = {'Exc','Inh','Mix'};
base_cell_name = 'Inh'; %this should be moved earlier in code for user to specify
for curr_cell = 1:3
    figure 
    
    %plot signal near electrode
    err = nanstd(corr_dist{curr_cell,1},[],2)./sqrt(sum(~isnan(corr_dist{curr_cell,1}),2));
    errorbar([1:5],nanmean(corr_dist{curr_cell,1},2),err,'LineWidth',2)
    hold on
    
    %plot signal far from electrode
    err2 = nanstd(corr_dist{curr_cell,2},[],2)./sqrt(sum(~isnan(corr_dist{curr_cell,2}),2));
    errorbar([1:5],nanmean(corr_dist{curr_cell,2},2),err2,'LineWidth',2)
    
    %plot base near electrode
    err3 = nanstd(corr_dist_base{curr_cell,1},[],2)./sqrt(sum(~isnan(corr_dist_base{curr_cell,1}),2));
    errorbar([1:5],nanmean(corr_dist_base{curr_cell,1},2),err3,'LineWidth',2)
    
    %plot base far from electrode
    err4 = nanstd(corr_dist_base{curr_cell,2},[],2)./sqrt(sum(~isnan(corr_dist_base{curr_cell,2}),2));
    errorbar([1:5],nanmean(corr_dist_base{curr_cell,2},2),err4,'LineWidth',2)
    
    %plot white noise near electrode
    err = nanstd(corr_dist_wn{curr_cell,1},[],2)./sqrt(sum(~isnan(corr_dist_base_wn{curr_cell,1}),2));
    errorbar([1:5],nanmean(corr_dist_wn{curr_cell,1},2),err,'LineWidth',2)
    
    %plot white noise far from electrode
    err = nanstd(corr_dist_wn{curr_cell,2},[],2)./sqrt(sum(~isnan(corr_dist_base_wn{curr_cell,2}),2));
    errorbar([1:5],nanmean(corr_dist_wn{curr_cell,2},2),err,'LineWidth',2)        
    
    title([base_cell_name,'-',cell_types{curr_cell}, ' xcorr'])
    ylim([-0.1,1])
    xlim([0,6])
    ax = gca;
    ax.XTick = [1:5];
    ax.XTickLabel = {'0-25','25-50','50-100', '100-200', '200+'};
    xlabel('Distance bin')
    ylabel('Correlation value')
    legend('Stim-evoked <200 um from elec', 'Stim-evoked >200 um from elec', 'Baseline <200 um from elec', 'Baseline >200 um from elec', 'White noise <200 um from elec', 'White noise >200 um from elec')
end

%% plot for base
cell_types = {'Exc','Inh','Mix'};
base_cell_name = 'Exc'; %this should be moved earlier in code for user to specify
for curr_cell = 1:3
    figure 
    plot([1:5], nanmean(corr_dist_base{curr_cell,1},2),'LineWidth',2)
    hold on
    err = nanstd(corr_dist_base{curr_cell,1},[],2)./sqrt(sum(~isnan(corr_dist_base{curr_cell,1}),2));
    errorbar([1:5],nanmean(corr_dist_base{curr_cell,1},2),err)
    
    plot([1:5], nanmean(corr_dist_base{curr_cell,2},2),'LineWidth',2)
    err2 = nanstd(corr_dist_base{curr_cell,2},[],2)./sqrt(sum(~isnan(corr_dist_base{curr_cell,2}),2));
    errorbar([1:5],nanmean(corr_dist_base{curr_cell,2},2),err2)
    
    title([base_cell_name,'-',cell_types{curr_cell}, ' xcorr'])
    %ylim([0.5,1])
    xlim([0,6])
    ax = gca;
    ax.XTick = [1:5];
    ax.XTickLabel = {'0-25','25-50','50-100', '100-200', '200+'};
    xlabel('Distance bin')
    ylabel('Correlation value')
    legend('Within 200 um of elec', 'Outside 200 um of elec')
end

%% can maybe extract labels here, but ultimately might need to write a
%separate function to look at correlations
%NEED TO MAKE THIS SECTION WORK
bins = [0,25,50,100,200,540];
big_bins = [0,200,540];
for curr_big_bin = 1:2
    for curr_cell = 1:3
        %label_perc{curr_cell,curr_big_bin} = NaN(5,5);
        for curr_bin = 1:5
            for curr_label = 1:5
                label_idx = base_label == curr_label;
                temp = [];
                for curr_animal = 1:4
                    big_bin_idx = active_neur_dist_split{base_cell,curr_train}(:,curr_animal) > big_bins(curr_big_bin) & active_neur_dist_split{base_cell,curr_train}(:,curr_animal) < big_bins(curr_big_bin+1);
                    %this will put zeros where there are NaNs but that should
                    %be okay
                    temp_diff_array = diff_dist{curr_animal,curr_cell};
                    temp_diff_array(~big_bin_idx,:) = NaN; %remove ones not in current big bin (<200 or >200)
                    dist_idx = temp_diff_array > bins(curr_bin) & temp_diff_array <= bins(curr_bin+1);
                    dist_idx(~label_idx(curr_animal,:),:) = 0; %remove ones that do not correspond to current label
                    temp = [temp, comp_label{curr_animal,curr_cell}(dist_idx)'];
                end
                corr_label{curr_cell,curr_big_bin,curr_label}(curr_bin,:) = [temp, NaN(3000-length(temp),1)']; 
            end
        end
    end
end