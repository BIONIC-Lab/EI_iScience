%set up load directory
base_dir ='C:\Users\chugh\Documents\GitHub\ExcitatoryInhibitory\Averages';
cd(base_dir)
d = dir;

frate = 0.034; %frame rate
x_label = [1:1800]*frate; %time axis
onset_start = 294; %approximate frame that ICMS started
onset_stop = floor(onset_start+3/frate);
s20 = floor(onset_start+20/frate);
offset_start = 1177;  %approximate frame ICMS ended
offset_stop = floor(offset_start+3/frate);
num_animals = 4; %total number of animals
animals = {'EIM08', 'EIM09', 'VGAT2','VGAT5'}; %animal IDs
train_types = {'10 Hz', '10 Hz Burst', '100 Hz', 'TBS'}; %different ICMS trains used 
cell_types = {'Exc','Inh','Mix'};

%initializing data arrays - MAY WANT TO CONSIDER ONLY USING THE SPLIT BY
%ANIMAL
%4 trains, 3 populations (exc, inh, and mix)
all_data = cell(4,3); %will save all data
all_data_norm = cell(4,3); %will save all data normalized (may not need this one)
all_data_split = cell(4,num_animals,3); %save all data split by animal 
error_all = cell(4,3); %save the standard deviation for all data
error_all_split = cell(4,num_animals,3); %save the standard deviation for all data split by animal

act_hist_data = [];

%go through data and load in to data arrays for processing
cd('C:\Users\chugh\Documents\GitHub\ExcitatoryInhibitory\scripts')
load('act_idx_save.mat')
for curr_train = 1:4
    
    %set directory for green channel data
    curr_dir = [base_dir, '\', d(curr_train+2).name];
    cd([curr_dir, '\Green'])
    d_green = dir;

    %set directory for red channel data
    cd([curr_dir, '\Red'])
    d_red = dir;

    %load in all animal data and correct/sort and combine into single arrays
    hist_plot = [];
    
    for curr_animal = 3:length(d_green)
        base_dir2 = 'C:\Users\chugh\Documents\GitHub\ExcitatoryInhibitory';
        stimpatterns = {'10 Hz', '10 Hz Burst', '100 Hz', 'TBS'};
        
        %call to function that labels inhibitory neurons based on
        %regression correction with baseline in green channel
        [exc_idx{curr_animal},inh_idx{curr_animal},oth_idx{curr_animal},hist_temp,meangreen{curr_animal},meanred{curr_animal},~] = geteimask(base_dir2,animals{curr_animal-2},1,onset_start*frate,frate,0,stimpatterns,1800);
%         [mdl_fit,intercept] = fix_green(base_dir2,animals{curr_animal-2},onset_start*frate,frate,stimpatterns);

        %modulate act_idx_save based on e,i,m idx
        act_idx_save{curr_animal-2}(:,inh_idx{curr_animal}) = NaN;
        act_idx_save{curr_animal-2}(:,oth_idx{curr_animal}) = NaN;

        hist_plot = [hist_plot, hist_temp];
        
        %load in green channel for current animal
        cd([curr_dir, '\Green'])
        m_green = readmatrix(d_green(curr_animal).name); 
        m_green = m_green(:,2:end);
        
%         %fix m_green based on model
%         m_green = (m_green-mdl_fit')+intercept;

        %filter data for noise
        m_green = smoothdata(m_green,'Gaussian',[7,7]);

        %read in red channel
        cd([curr_dir, '\Red'])
        m_red = readmatrix(d_red(curr_animal).name);
        m_red = m_red(:,2:end);

        %calculate baseline for green channel
        baseline_green = nanmean(m_green(1:onset_start,:)); %find mean across all ROIs during baseline period

        %calculate dF/F0 of green channel
        m_green_norm = (m_green-baseline_green)./abs(baseline_green);
        
        %only use neurons modulated i.e. with some significant increase of
        %decrease
        curr_error = std(m_green_norm(1:floor((10/frate)),:),[],1); %std for each ROI in first 10 s 
        baseline_mean{curr_train,curr_animal-2} = mean(m_green_norm(1:floor((10/frate)),:),1); %baseline mean

        %active populations will be greater than std*5
%         act_idx = max(abs(m_green_norm(onset_start:offset_start,:))) > abs(baseline_mean{curr_train,curr_animal-2})+curr_error*5;

  
        %can also limit by act_idx_save
        save_idx = act_idx_save{curr_animal-2}(1,:) == 2 & act_idx_save{curr_animal-2}(2,:) == 2 & act_idx_save{curr_animal-2}(3,:) == 2 & act_idx_save{curr_animal-2}(4,:) == 3;
        
        %implement additional criteria to ensure that activity is above for
        %at least 15 consecutive samples - will also do less than in
        %separate analysis
        clear act_idx
        act_idx_pre = abs(m_green_norm(onset_start:offset_start,:)) > abs(baseline_mean{curr_train,curr_animal-2})+curr_error*3;
        act_idx_all = any(act_idx_pre);
        for curr_neuron = 1:size(act_idx_pre,2)
            act_idx(curr_neuron) = numel(strfind(act_idx_pre(:,curr_neuron)',ones(1,15)));
        end
        act_idx(act_idx>1) = 1;
%         act_idx = save_idx; %just use save_idx to focus on specific groups
      
        %if we want to look at non-consecutive active neurons
%         act_idx = ~act_idx & act_idx_all;         
        %if we want to look at inactive neurons
%         act_idx = ~act_idx & ~act_idx_all;
        


%         act_idx_save{1,curr_animal-2}(curr_train,act_idx) = 3;

        act_idx_raster{curr_train,curr_animal-2,1} = abs(m_green_norm(:,act_idx&exc_idx{curr_animal}')) > abs(baseline_mean{curr_train,curr_animal-2}(:,act_idx&exc_idx{curr_animal}'))+curr_error(act_idx&exc_idx{curr_animal}')*3;
        act_idx_raster{curr_train,curr_animal-2,2} = abs(m_green_norm(:,act_idx&inh_idx{curr_animal}')) > abs(baseline_mean{curr_train,curr_animal-2}(:,act_idx&inh_idx{curr_animal}'))+curr_error(act_idx&inh_idx{curr_animal}')*3;
        act_idx_raster{curr_train,curr_animal-2,3} = abs(m_green_norm(:,act_idx&oth_idx{curr_animal}')) > abs(baseline_mean{curr_train,curr_animal-2}(:,act_idx&oth_idx{curr_animal}'))+curr_error(act_idx&oth_idx{curr_animal}')*3;

        
        act_hist_data = [act_hist_data, (max(abs(m_green_norm(onset_start:offset_start,:)))-abs(baseline_mean{curr_train,curr_animal-2}))./curr_error];
        m_green_norm(:,~act_idx) = NaN; %turns inactive neurons into NaNs
        
        %can exclude data based on act_idx across trains to look at
        %specific kinds of neurons
        
        %active for 10 Hz B but inactive for TBS
%         act_idx2 = act_idx_save{1,curr_animal-2}(2,:) == 2 & act_idx_save{1,curr_animal-2}(4,:) == 1;
%         m_green_norm(:,~act_idx2) = NaN; %turns inactive neurons into NaNs
        
        
        m_green_norm2 = m_green_norm./max(m_green_norm); %normalize to max on each ROI
        
%         figure
%         plot(x_label, nanmean(m_red,2),'r','LineWidth',1)
%         hold on
%         plot(x_label, nanmean(m_green,2),'g','LineWidth',1)
%         hold on
%         xline(10,'--','LineWidth',2)
%         xline(40,'--','LineWidth',2)
%         title('Red Channel Intensity')
%         xlim([0 60])
%         %ylim([-0.2 1.2])
%         xlabel('Time (s)')
%         ylabel('Intensity (F)')
        
        %find the difference between mean intensity for exc vs inh and save
        %for each animal after removing inactive neurons
        m_diff(curr_train,curr_animal-2) = sum(nanmean(m_green_norm(:,inh_idx{curr_animal}),2) - nanmean(m_green_norm(:,exc_idx{curr_animal}),2));
        num_inh(curr_train,curr_animal-2) = sum(inh_idx{curr_animal}&act_idx');
        num_exc(curr_train,curr_animal-2) = sum(exc_idx{curr_animal}&act_idx');
        m_diff_e(curr_train,curr_animal-2) = sum(nanmean(m_green_norm(:,exc_idx{curr_animal}),2));
        m_diff_i(curr_train,curr_animal-2) = sum(nanmean(m_green_norm(:,inh_idx{curr_animal}),2));

        %parse data into bigger arrays for all data
        %place excitatory neurons in first slot
        all_data{curr_train,1} = [all_data{curr_train,1}; m_green_norm(:,exc_idx{curr_animal})'];
        all_data_norm{curr_train,1} = [all_data_norm{curr_train,1}; m_green_norm2(:,exc_idx{curr_animal})'];
        all_data_split{curr_train,curr_animal-2,1} = m_green_norm';
        all_data_split{curr_train,curr_animal-2,1}(~exc_idx{curr_animal},:) = NaN;        
        %error is std of baseline only
        error_all{curr_train,1} = [error_all{curr_train,1}; curr_error(exc_idx{curr_animal})'];
        error_all_split{curr_train,curr_animal-2,1} = curr_error';
        error_all_split{curr_train,curr_animal-2,1}(~exc_idx{curr_animal}) = NaN;

        %place inhibitory neurons in second slot
        all_data{curr_train,2} = [all_data{curr_train,2}; m_green_norm(:,inh_idx{curr_animal})'];
        all_data_norm{curr_train,2} = [all_data_norm{curr_train,2}; m_green_norm2(:,inh_idx{curr_animal})'];
        all_data_split{curr_train,curr_animal-2,2} = m_green_norm';
        all_data_split{curr_train,curr_animal-2,2}(~inh_idx{curr_animal},:) = NaN;
        error_all{curr_train,2} = [error_all{curr_train,2}; curr_error(inh_idx{curr_animal})'];
        error_all_split{curr_train,curr_animal-2,2} = curr_error';
        error_all_split{curr_train,curr_animal-2,2}(~inh_idx{curr_animal}) = NaN;

        %place mixed neurons in third slot
        all_data{curr_train,3} = [all_data{curr_train,3}; m_green_norm(:,oth_idx{curr_animal})'];
        all_data_norm{curr_train,3} = [all_data_norm{curr_train,3}; m_green_norm2(:,oth_idx{curr_animal})'];
        all_data_split{curr_train,curr_animal-2,3} = m_green_norm';
        all_data_split{curr_train,curr_animal-2,3}(~oth_idx{curr_animal},:) = NaN;
        error_all{curr_train,3} = [error_all{curr_train,3}; curr_error(oth_idx{curr_animal})'];
        error_all_split{curr_train,curr_animal-2,3} = curr_error';
        error_all_split{curr_train,curr_animal-2,3}(~oth_idx{curr_animal}) = NaN;

        %plot the mean response of each neuron group for each animal on a
        %separate plot
%         figure
%         plot(x_label, nanmean(m_green_norm(:,exc_idx),2),'b','LineWidth',2)
%         hold on
%         plot(x_label, nanmean(m_green_norm(:,inh_idx),2),'r','LineWidth',2)
%         plot(x_label, nanmean(m_green_norm(:,oth_idx),2),'m','LineWidth',2)
%         yline(0, '--', 'LineWidth', 2)
%         xline(10,'--','LineWidth',2)
%         xline(40,'--','LineWidth',2)
%         legend('Exc','Inh')
%         title(animals{curr_animal-2})
%         xlim([0 60])
    %     ylim([-0.2 1.2])


        %look at average baseline to see if there are different
        %should probably make this across all animals
%         figure
%         base_exc = [baseline_green(exc_idx)'; NaN(3000-sum(exc_idx),1)];
%         base_inh = [baseline_green(inh_idx)'; NaN(3000-sum(inh_idx),1)];
%         base_oth = [baseline_green(oth_idx)'; NaN(3000-sum(oth_idx),1)];
%         violin([base_exc,base_inh,base_oth])
%         title('Baseline')
%         ax = gca;
%         ax.XTick = [1,2];
%         ax.XTickLabel = {'Excitatory', 'Inhibitory'};

        %save the total percent of each population that is active in each
        %animal
        perc_act{1}(curr_train,curr_animal-2) = sum(~isnan(m_green_norm(1,exc_idx{curr_animal})))/(sum(exc_idx{curr_animal}));
        perc_act{2}(curr_train,curr_animal-2) = sum(~isnan(m_green_norm(1,inh_idx{curr_animal})))/(sum(inh_idx{curr_animal}));
        perc_act{3}(curr_train,curr_animal-2) = sum(~isnan(m_green_norm(1,oth_idx{curr_animal})))/(sum(oth_idx{curr_animal}));
        
        %save logical indices to find what percent of neurons are not
        %active for any train
        perc_log{1,curr_animal-2}(curr_train,:) = ~isnan(m_green_norm(1,exc_idx{curr_animal}));
        perc_log{2,curr_animal-2}(curr_train,:) = ~isnan(m_green_norm(1,inh_idx{curr_animal}));
        perc_log{3,curr_animal-2}(curr_train,:) = ~isnan(m_green_norm(1,oth_idx{curr_animal}));
    end
    %plot histogram of red/green across animals
    figure
    histogram(hist_plot,50) %looks like division could be around 5.3?

    %plot average profiles for each subpopulation averaged across all
    %animals
    figure
    mean_exc = nanmean(all_data{curr_train,1}(:,onset_start:onset_start+29),1);
    plot(x_label(1:30), mean_exc,'b','LineWidth',2)
    hold on
    mean_inh = nanmean(all_data{curr_train,2}(:,onset_start:onset_start+29),1);
    plot(x_label(1:30), mean_inh,'r','LineWidth',2)
    mean_mix = nanmean(all_data{curr_train,3}(:,onset_start:onset_start+29),1);
    plot(x_label(1:30), mean_mix,'m','LineWidth',2)
    yline(0, '--', 'LineWidth', 2)
    %mark std 
    yline(nanmean(error_all{curr_train,1}), 'b')
    yline(nanmean(error_all{curr_train,2}), 'r')
    yline(nanmean(error_all{curr_train,3}), 'm')
    
    %mark location of peaks
    [~,loc] = max(mean_exc);
    xline(loc*frate,'b')
    [~,loc] = max(mean_inh);
    xline(loc*frate,'r')
    [~,loc] = max(mean_mix);
    xline(loc*frate,'m')

    xline(10,'--','LineWidth',2)
    xline(40,'--','LineWidth',2)
    legend('Exc','Inh', 'Mixed')
    title(train_types(curr_train))
    xlim([0 60])
    %ylim([-0.2 1.2])
    xlabel('Time (s)')
    ylabel('dF/F0')
end

%calculate percent division into each category
total_cnt = zeros(1,3);
for curr_cell = 1:3
    for curr_train = 1:4
        total_cnt(curr_cell) = total_cnt(curr_cell)+sum(~isnan(all_data{curr_train,curr_cell}(:,1)));
    end
end
figure
pie(total_cnt)
legend('Exc','Inh','Mix')

%do stats between cell types for percent active
for curr_train = 1:4
    perc_cell{curr_train}(:,1) = [perc_act{1}(curr_train,:), NaN(1,3000-length(perc_act{1}(curr_train,:)))];
    perc_cell{curr_train}(:,2) = [perc_act{2}(curr_train,:), NaN(1,3000-length(perc_act{2}(curr_train,:)))];
    perc_cell{curr_train}(:,3) = [perc_act{3}(curr_train,:), NaN(1,3000-length(perc_act{3}(curr_train,:)))];
    [p_perc(curr_train),tbl_perc{curr_train},stats_perc{curr_train}] = anova1(perc_cell{curr_train});
    mp_perc{curr_train} = multcompare(stats_perc{curr_train});
end

%make bar plot for percent active for each cell type
for curr_cell = 1:3
    bar_data(:,curr_cell) = nanmean(perc_act{curr_cell},2);
    bar_error(:,curr_cell) = nanstd(perc_act{curr_cell},[],2)/sqrt(4);
end
figure
b = bar(bar_data);
hold on
ngroups = size(bar_data,1);
nbars = size(bar_data,2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars+1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, bar_data(:,i), bar_error(:,i), '.');
end
hold off
ax = gca;
ax.XTick = [1:4];
ax.XTickLabels = train_types;
ylim([0 1])
ax.YTickLabels = {'0','10','20','30','40','50','60','70','80','90','100'};
ylabel('Percent of ROIs active')
legend('Exc','Inh','Mix')

%can also make plot and stats for combining across trains
comb = ([perc_cell{1,1}(:,1)',perc_cell{1,2}(:,1)',perc_cell{1,3}(:,1)',perc_cell{1,4}(:,1)'])';
comb(:,2) = ([perc_cell{1,1}(:,2)',perc_cell{1,2}(:,2)',perc_cell{1,3}(:,2)',perc_cell{1,4}(:,2)'])';
comb(:,3) = ([perc_cell{1,1}(:,3)',perc_cell{1,2}(:,3)',perc_cell{1,3}(:,3)',perc_cell{1,4}(:,3)'])';

bar_data2(1) = nanmean([perc_cell{1,1}(:,1)',perc_cell{1,2}(:,1)',perc_cell{1,3}(:,1)',perc_cell{1,4}(:,1)']);
bar_error2(1) = nanstd([perc_cell{1,1}(:,1)',perc_cell{1,2}(:,1)',perc_cell{1,3}(:,1)',perc_cell{1,4}(:,1)'])/sqrt(12);

bar_data2(2) = nanmean([perc_cell{1,1}(:,2)',perc_cell{1,2}(:,2)',perc_cell{1,3}(:,2)',perc_cell{1,4}(:,2)']);
bar_error2(2) = nanstd([perc_cell{1,1}(:,2)',perc_cell{1,2}(:,2)',perc_cell{1,3}(:,2)',perc_cell{1,4}(:,2)'])/sqrt(12);

bar_data2(3) = nanmean([perc_cell{1,1}(:,3)',perc_cell{1,2}(:,3)',perc_cell{1,3}(:,3)',perc_cell{1,4}(:,3)']);
bar_error2(3) = nanstd([perc_cell{1,1}(:,3)',perc_cell{1,2}(:,3)',perc_cell{1,3}(:,3)',perc_cell{1,4}(:,3)'])/sqrt(12);

[p_comb,tbl_comb,stats_comb] = anova1(comb);
mp_comb = multcompare(stats_comb);


figure
bar(bar_data2)
hold on
e = errorbar([1:3], bar_data2, bar_error2);
e.LineStyle = 'none';
ylim([0 1])
ax = gca;
ax.YTickLabel = {0:10:100};
ax.XTickLabel = cell_types;

%significant between trains?
anova1([perc_cell{1,1}(:,1),perc_cell{1,2}(:,1),perc_cell{1,3}(:,1),perc_cell{1,4}(:,1)]) %excitatory neurons
anova1([perc_cell{1,1}(:,2),perc_cell{1,2}(:,2),perc_cell{1,3}(:,2),perc_cell{1,4}(:,2)]) %inhibitory neurons
anova1([perc_cell{1,1}(:,3),perc_cell{1,2}(:,3),perc_cell{1,3}(:,3),perc_cell{1,4}(:,3)]) %mixed neurons

%determine what percent of neurons are not active for any trains 

%excitatory
comb_array{1} = [perc_log{1,1},perc_log{1,2},perc_log{1,3},perc_log{1,4}];
perc_not_act(1) = sum(sum(comb_array{1})==0)/length(comb_array{1});

%inhibitory
comb_array{2} = [perc_log{2,1},perc_log{2,2},perc_log{2,3},perc_log{2,4}];
perc_not_act(2) = sum(sum(comb_array{2})==0)/length(comb_array{2});

%calculate mean intensity difference in first 10 s and last 20 s across
%subtypes for each train
for curr_train = 1:4
    %mean in first 10 s
    temp = nanmean(all_data{curr_train,1}(:,onset_start:onset_start+round(15/frate)),2);
    mean_10s{curr_train}(:,1) = [temp; NaN(3000-length(temp),1)];
    temp = nanmean(all_data{curr_train,2}(:,onset_start:onset_start+round(15/frate)),2);
    mean_10s{curr_train}(:,2) = [temp; NaN(3000-length(temp),1)];
    temp = nanmean(all_data{curr_train,3}(:,onset_start:onset_start+round(15/frate)),2);
    mean_10s{curr_train}(:,3) = [temp; NaN(3000-length(temp),1)];
    
    %mean in last 20 s
    temp = nanmean(all_data{curr_train,1}(:,onset_start+round(15/frate):offset_start),2);
    mean_30s{curr_train}(:,1) = [temp; NaN(3000-length(temp),1)];
    temp = nanmean(all_data{curr_train,2}(:,onset_start+round(15/frate):offset_start),2);
    mean_30s{curr_train}(:,2) = [temp; NaN(3000-length(temp),1)];
    temp = nanmean(all_data{curr_train,3}(:,onset_start+round(15/frate):offset_start),2);
    mean_30s{curr_train}(:,3) = [temp; NaN(3000-length(temp),1)];
end
%% raster plot 
%can potentially color by neuron type and then separate animals with horz
%line
colors_raster = [0 0 1; 1 0 0; 1 0 1]; %label by animal
%all_raster = [];
for curr_train = 1:4
    figure
    hold on
    all_raster{curr_train} = [];
    y = 0;
    for curr_animal = 1:4
        for curr_cell = 1:3
            [x_curr,y_curr] = find(act_idx_raster{curr_train,curr_animal,curr_cell}==1);
            all_raster{curr_train} = [all_raster{curr_train} , act_idx_raster{curr_train,curr_animal,curr_cell}];
            x = x_curr*frate;
            if ~isempty(y_curr)
                y = y_curr+max(y); %add to yaxis
                scatter(x,y,'.', 'MarkerEdgeColor', colors_raster(curr_cell,:))
            end
        end
        y = y + 1;
        yline(max(y), '--', 'LineWidth', 2)
    end
    xline(10, '--', 'LineWidth', 2)
    xline(40, '--', 'LineWidth', 2)
    title(train_types{curr_train})
    legend(animals)
end

figure
hold on
colors_hist = [255,0,255; 250,200,152; 255,250,160; 119,221,119]./255; 
rast_save = [];
for curr_train = 1:4
    sum_all_rast = sum(all_raster{curr_train},2);
    bin_size = 1;
    binned_all_rast = squeeze(sum(reshape(sum_all_rast,size(sum_all_rast,2),bin_size,[]),2));
    rast_save = [rast_save, binned_all_rast];
    %bar(x_label,binned_all_rast, 'FaceColor', colors_hist(curr_train,:))
    smooth_sum_rast = smoothdata(binned_all_rast,'Gaussian',[30,30]);
    hold on
    plot(x_label,smooth_sum_rast,'LineWidth',2, 'Color', colors_hist(curr_train,:))
end
yyaxis right
bar(x_label,sum(rast_save,2), 'FaceColor', 'k')
xline(10,'--','LineWidth',2)
xline(40,'--','LineWidth',2)
ylabel('Threshold crossing count')
xlabel('Time (s)')
legend(train_types)

k = mean(binned_all_rast);
v = var(binned_all_rast);
opt_val = (2*k*v)/bin_size;

%% can calculate mean intensity in first 3 s and determine preferential response and then mean intensity across whole train
%can switch between using "no pref" by commenting in last section of loop
%this is not currently included in paper draft but I think it probably
%should be, maybe before we break up neurons into individual responses
%maybe also want to look at how individual neuron class relates to
%preferential responses
%also may want to try it only with trains that had activity for all trains
clear mean_3s max_vals std_3s m3s m3s_split
for curr_cell = 1:3
    for curr_train = 1:4
        temp = [];
        temp2 = [];
        for curr_animal = 1:4
            %need to 
            temp = [temp, nanmean(all_data_split{curr_train,curr_animal,curr_cell}(:,onset_start:onset_stop),2)'];
            temp2 = [temp2, nanstd(all_data_split{curr_train,curr_animal,curr_cell}(:,1:onset_start-1),[],2)'];
            mean_3s_split{curr_cell,curr_animal}(curr_train,:) = nanmean(all_data_split{curr_train,curr_animal,curr_cell}(:,onset_start:onset_stop),2);
            %std_split{curr_cell,curr_animal}(curr_train,:) = nanstd(all_data_split{curr_train,curr_animal,curr_cell}(:,1:onset_start-1),[],2); %calculate baseline std to determine preferences
            
            %if a neuron is active is determined by whole ICMS interval so
            %this should be the same for both short and long (unless we
            %want to split it with more code)
            %count_3s{curr_cell}(curr_animal,curr_train) = sum(~isnan((all_data_split{curr_train,curr_animal,curr_cell}(:,1))));
            %realizing this is not really any different than percent active
        end
        mean_3s{curr_cell}(curr_train,:) = temp;
        std_3s{curr_cell}(curr_train,:) = temp2;
        %count_3s(curr_cell,curr_train) = sum(~isnan(temp));
    end
    %get rid of nans
    nan_idx = all(isnan(mean_3s{curr_cell}));  %can eliminate neurons that were not active for all trains by changing "all" to "any"
    mean_3s{curr_cell} = mean_3s{curr_cell}(:,~nan_idx);
    std_3s{curr_cell} = std_3s{curr_cell}(:,~nan_idx);
    
    %find which train each ROI had the max value for
    %currently this includes ROIs that were inactive for one or more
    %trains, which could maybe bias things?
    
    [real_max,temp2] = max(abs(mean_3s{curr_cell}));
    max_vals{curr_cell} = temp2;
    
    %can we make a logical array for max_vals?
    max_idx = abs(mean_3s{curr_cell}) == real_max;
    
    %can instead use std to determine if they are significantly different
%     max_comp = mean_3s{curr_cell};
%     max_comp(max_idx) = NaN; %remove the maxes from each column
%     max_comp = max_comp+std_3s{curr_cell};
%     max_comp(isnan(max_comp)) = 0; %need to change to zeros because it should be greater than inactive neurons
%     max_vals_log{curr_cell} = max_comp<=real_max; %this will compare the max_value to each other response plus its standard deviation
%     max_all_idx = (all(max_vals_log{curr_cell}==1)); %checks if the max train is greater than all of the other trains by at least 1 std
%     max_vals{curr_cell}(max_all_idx) = 5; %change max to 5 if there is no real preference, which we can then label later 
end

%plot violins of average and pies for max 
for curr_cell = 1:3
    figure
    violin(mean_3s{curr_cell}')
    [p_3s(curr_cell),~,stats] = anova1(mean_3s{curr_cell}');
    m3s{curr_cell} = multcompare(stats);
    ax = gca;
    ax.XTick = [1:4];
    ax.XTickLabel = train_types;
    ylabel('Mean intensities first 3 s')
    title(cell_types(curr_cell))
   
    %can make pie charts for max vals
    figure
    for curr_train = 1:4
        max_cnt(curr_train) = sum(max_vals{curr_cell} == curr_train);
        for curr_animal = 1:4 %find percent from split data set
            nan_idx = all(isnan(mean_3s_split{curr_cell,curr_animal}));
            [~,max_vals_split{curr_cell,curr_animal}] = max(mean_3s_split{curr_cell,curr_animal}(:,~nan_idx));
            max_cnt_split{curr_cell}(curr_train,curr_animal) = nansum(max_vals_split{curr_cell,curr_animal}==curr_train);
            %we can also normalize these values based on the number of
            %exc,inh,mix
            max_cnt_split_norm{curr_cell}(curr_train,curr_animal) = nansum(max_vals_split{curr_cell,curr_animal}==curr_train)/length(max_vals_split{curr_cell,curr_animal});
        end
    end
    %stats
    [p_3s_split(curr_cell),~,stats] = anova1(max_cnt_split{curr_cell}');
    m3s_split{curr_cell} = multcompare(stats);
    
    %pie plot
    pie(max_cnt)
    legend([train_types,'No Pref'])
    title(cell_types(curr_cell))
    
    %bar plot
    figure
    bar(nanmean(max_cnt_split{curr_cell},2))
    hold on
    e = errorbar([1:4], nanmean(max_cnt_split{curr_cell},2), nanstd(max_cnt_split{curr_cell},[],2)./sqrt(size(max_cnt_split{curr_cell},2)));
    e.LineStyle = 'none';
    title(cell_types(curr_cell))
    ax = gca;
    ax.XTickLabels = train_types;
    ylabel('Neuron count')
end


%% do same thing for last 10 s
clear mean_30s mean_30s_split max_vals2
for curr_cell = 1:3
    for curr_train = 1:4
        temp = [];
        for curr_animal = 1:4
            %can change to not combine across animals to do stats
            temp = [temp, nanmean(all_data_split{curr_train,curr_animal,curr_cell}(:,s20:offset_start),2)'];
            %can insert labels_split here if we want to look only at
            %RA,SA,etc.t labels_split{curr_cell,curr_animal}(3,:)==1
            %or insert active neur dist
            temp3 = nanmean(all_data_split{curr_train,curr_animal,curr_cell}(:,s20:offset_start),2);
            mean_30s_split{curr_cell,curr_animal}(curr_train,:) = [temp3; NaN(3000-length(temp3),1)];
        end
        temp = temp(active_neur_dist{curr_cell,curr_train}<125);
        mean_30s{curr_cell}(curr_train,:) = [temp, NaN(1,3000-length(temp))];
    end
    nan_idx = all(isnan(mean_30s{curr_cell})); %find nans
    mean_30s{curr_cell} = mean_30s{curr_cell}(:,~nan_idx); %get rid of neurons that did not respond to any train
    [real_max,temp2] = max(mean_30s{curr_cell});
    max_vals2{curr_cell} = temp2;
    
    
    %can instead use std to determine if they are significantly different
%     max_idx = abs(mean_30s{curr_cell}) == real_max;
%     max_comp = mean_30s{curr_cell};
%     max_comp(max_idx) = NaN; %remove the maxes from each column
%     max_comp = max_comp+std_3s{curr_cell}; %std is the same for both versions
%     max_comp(isnan(max_comp)) = 0; %need to change to zeros because it should be greater than inactive neurons
%     max_vals_log{curr_cell} = max_comp<=real_max; %this will compare the max_value to each other response plus its standard deviation
%     max_all_idx = (all(max_vals_log{curr_cell}==1)); %checks if the max train is greater than all of the other trains by at least 1 std
%     max_vals2{curr_cell}(max_all_idx) = 5; %change max to 5 if there is no real preference, which we can then label later 
end

% plot violins of average for each train
%plot violins of average and pies for max 
for curr_cell = 1:3
%     figure
%     violin(mean_30s{curr_cell}')
%     [p_30s(curr_cell),~,stats] = anova1(mean_30s{curr_cell}');
%     m30s{curr_cell} = multcompare(stats);
%     ax = gca;
%     ax.XTick = [1:4];
%     ax.XTickLabel = train_types;
%     ylabel('Mean intensities first 3 s')
%     title(cell_types(curr_cell))
    
    %can make pie charts for max vals
    figure
    for curr_train = 1:4 %this was 5 before?
        max_cnt2(curr_train) = sum(max_vals2{curr_cell} == curr_train);
        for curr_animal = 1:4 %find percent from split data set
            nan_idx = all(isnan(mean_30s_split{curr_cell,curr_animal}));
            
            %can add dist_idx here to separate by dist bins
            % have to get distance from lower code
%             for curr_train_2 = 1:4
%                 dist_idx(curr_train_2,:) = (active_neur_dist_split{curr_cell,curr_train_2}(:,curr_animal) > 125 & active_neur_dist_split{curr_cell,curr_train_2}(:,curr_animal) < 200)';
%             end
            
            [~,max_vals2_split{curr_cell,curr_animal}] = max(mean_30s_split{curr_cell,curr_animal}(:,~nan_idx));
            max_cnt2_split{curr_cell}(curr_train,curr_animal) = nansum(max_vals2_split{curr_cell,curr_animal}==curr_train);
            max_cnt2_split_norm{curr_cell}(curr_train,curr_animal) = nansum(max_vals2_split{curr_cell,curr_animal}==curr_train)/length(max_vals2_split{curr_cell,curr_animal});
        end
    end
    %stats
    [p_30s_split(curr_cell),~,stats] = anova1(max_cnt2_split{curr_cell}');
    m30s_split{curr_cell} = multcompare(stats);
    
    %pie plot
    pie(max_cnt2)
    legend(train_types)
    title(cell_types(curr_cell))
    
    %bar plot
    figure
    bar(nanmean(max_cnt2_split{curr_cell},2))
    hold on
    e = errorbar([1:4], nanmean(max_cnt2_split{curr_cell},2), nanstd(max_cnt2_split{curr_cell},[],2)./sqrt(size(max_cnt2_split{curr_cell},2)));
    e.LineStyle = 'none';
    title(cell_types(curr_cell))
    ax = gca;
    ax.XTickLabels = train_types;
    ylabel('Neuron count')
    
end

%bar plot comparing percent for TBS  exc vs. inh
figure
exc_mean = nanmean(max_cnt2_split{1}(4,:)./sum(max_cnt2_split{1},1));
inh_mean = nanmean(max_cnt2_split{2}(4,:)./sum(max_cnt2_split{2},1));
bar([exc_mean,inh_mean])
hold on
exc_err = nanstd(max_cnt2_split{1}(4,:)./sum(max_cnt2_split{1},1))./sqrt(4);
inh_err = nanstd(max_cnt2_split{2}(4,:)./sum(max_cnt2_split{2},1))./sqrt(4);
e = errorbar([1:2], [exc_mean,inh_mean], [exc_err,inh_err]);
e.LineStyle = 'none';
title('TBS: Exc vs. Inh')
ax = gca;
ax.XTickLabels = {'Exc','Inh'};
ylabel('Percent preference')
anova1([(max_cnt2_split{1}(4,:)./sum(max_cnt2_split{1},1))', (max_cnt2_split{2}(4,:)./sum(max_cnt2_split{2},1))'])

%bar plot comparing percent for 10 Hz Burst exc vs. inh
figure
exc_mean = nanmean(max_cnt2_split{1}(2,:)./sum(max_cnt2_split{1},1));
inh_mean = nanmean(max_cnt2_split{2}(2,:)./sum(max_cnt2_split{2},1));
bar([exc_mean,inh_mean])
hold on
exc_err = nanstd(max_cnt2_split{1}(2,:)./sum(max_cnt2_split{1},1))./sqrt(4);
inh_err = nanstd(max_cnt2_split{2}(2,:)./sum(max_cnt2_split{2},1))./sqrt(4);
e = errorbar([1:2], [exc_mean,inh_mean], [exc_err,inh_err]);
e.LineStyle = 'none';
title('10 Hz Burst: Exc vs. Inh')
ax = gca;
ax.XTickLabels = {'Exc','Inh'};
ylabel('Percent preference')
anova1([(max_cnt2_split{1}(2,:)./sum(max_cnt2_split{1},1))', (max_cnt2_split{2}(2,:)./sum(max_cnt2_split{2},1))'])

%ranova stuff
% Meas = dataset([1 2 3 4]','VarNames',{'Measurements'});
% t = table([1:255]',mean_test(:,1),mean_test(:,2),mean_test(:,3),mean_test(:,4),...
% 'VariableNames',{'neurons','train1','train2','train3','train4'});
% rm = fitrm(t, 'train1-train4~neurons','WithinDesign',Meas);
% tbl = multcompare(rm,'Measurements');

%% calculate time to peak and time to baseline
%should double check to be sure but things look about right
%also should add stats

%first time to peak
for curr_cell = 1:3
    for curr_train = 1:4
        [curr_peak, loc] = max(all_data{curr_train,curr_cell}(:,onset_start:offset_start),[],2); 
        loc(isnan(curr_peak)) = NaN; %otherwise defaults to 1
        time_to_peak{curr_cell}(:,curr_train) = loc; 
    end
    pk_avg(:,curr_cell) = nanmean(time_to_peak{curr_cell},1)*frate;
    pk_err(:,curr_cell) = nanstd(time_to_peak{curr_cell})./(sqrt(sum(~isnan(time_to_peak{curr_cell}))))*frate;
end

%do stats between cell types
for curr_train = 1:4
    pk_cell{curr_train}(:,1) = [time_to_peak{1}(:,curr_train); NaN(3000-length(time_to_peak{1}(:,curr_train)),1)];
    pk_cell{curr_train}(:,2) = [time_to_peak{2}(:,curr_train); NaN(3000-length(time_to_peak{2}(:,curr_train)),1)];
    pk_cell{curr_train}(:,3) = [time_to_peak{3}(:,curr_train); NaN(3000-length(time_to_peak{3}(:,curr_train)),1)];
    [pp(curr_train),tbl_p{curr_train},stats_p{curr_train}] = anova1(pk_cell{curr_train});
    mpp{curr_train} = multcompare(stats_p{curr_train});
end

figure
bar(pk_avg)
hold on
ngroups = size(pk_avg,1);
nbars = size(pk_avg,2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars+1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, pk_avg(:,i), pk_err(:,i), '.');
end
ax = gca;
ax.XTick = [1:4];
ax.XTickLabels = train_types;
ylabel('Time to peak (s)')

%% time to baseline
for curr_cell = 1:3
    clear loc
    for curr_train = 1:4
        baseline_std = nanstd(all_data{curr_train,curr_cell}(:,1:onset_start),[],2);
        baseline_mean = nanmean(all_data{curr_train,curr_cell}(:,1:onset_start),2);
        for curr_idx = 1:size(all_data{curr_train,curr_cell},1)
            if ~isnan(all_data{curr_train,curr_cell}(curr_idx,:))
                curr_loc = find(all_data{curr_train,curr_cell}(curr_idx,offset_start:end)<(0.5*baseline_std(curr_idx)),1,'first');
                if ~isempty(curr_loc)
                    loc(curr_idx) =  curr_loc;
                else
                    loc(curr_idx) = 30/frate; %if it doesn't reach baseline by end of 20 s, set to max
                end
            else
                loc(curr_idx) =  NaN; %if cell was inactive
            end
        end        
        time_to_base{curr_cell}(:,curr_train) = loc;
    end
    base_avg(:,curr_cell) = nanmean(time_to_base{curr_cell},1)*frate;
    base_err(:,curr_cell) = nanstd(time_to_base{curr_cell})./(sqrt(sum(~isnan(time_to_base{curr_cell}))))*frate;
end

%do stats between cell types
for curr_train = 1:4
    base_cell{curr_train}(:,1) = [time_to_base{1}(:,curr_train); NaN(3000-length(time_to_base{1}(:,curr_train)),1)];
    base_cell{curr_train}(:,2) = [time_to_base{2}(:,curr_train); NaN(3000-length(time_to_base{2}(:,curr_train)),1)];
    base_cell{curr_train}(:,3) = [time_to_base{3}(:,curr_train); NaN(3000-length(time_to_base{3}(:,curr_train)),1)];
    [pb(curr_train),tbl_b{curr_train},stats_b{curr_train}] = anova1(base_cell{curr_train});
    mpb{curr_train} = multcompare(stats_b{curr_train});
end


figure
bar(base_avg)
hold on
ngroups = size(base_avg,1);
nbars = size(base_avg,2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars+1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, base_avg(:,i), base_err(:,i), '.');
end
ax = gca;
ax.XTick = [1:4];
ax.XTickLabels = train_types;
ylabel('Time to baseline (s)')

%% pca (not included in paper)
% data = [all_data{3,1}]';
% %labels = [repmat(1,1,size(all_data{3,1},1)), repmat(2,1,size(all_data{3,2},1)), repmat(3,1,size(all_data{3,3},1))];
% ks = 5;
% 
% %have to remove nans
% idx = ~isnan(data);
% data_new = [];
% %labels_new = [];
% for i = 1:size(idx,2)
%     data_new = [data_new, data(idx(:,i),i)];
%     %labels_new = [labels_new, labels(idx(1,i),i)];
% end
% 
% [c,s,l,~,explained] = pca(data_new);
% 
% %remove outliers 
% %     sum_idx = sum(c(:,1:3),2);
% %     c = c(sum_idx<0.5,:);
% 
% %k-means cluster on pca results
% [idx2, cd, sumd] = kmeans(c(:,1:2),ks);
% 
% %plot results
% figure
% for curr_fac = 1:ks
%     curr_idx = idx2 == curr_fac;
%     %plot3(c(curr_idx,1), c(curr_idx ,2), c(curr_idx ,3), '.')
%     scatter(c(curr_idx,1),c(curr_idx,2),'.')
%     hold on
% end
% 
% figure
% for curr_fac = 1:ks
%     curr_idx = idx2 == curr_fac;
%     plot(mean(data_new(:,curr_idx),2))
%     hold on
% end

%% try to classify neurons
%1 - based on stable, depressing, facilitating or decreased activity during
%stim
%2 - based on activity after stim - return to baseline, go below baseline,
%offset activity i.e. increased activity

%% during stim
[labels, max_peak_all, adapt_all] = labels_during_stim(all_data,error_all);
label_types = {'RA','SA','NA','Fac','Decrease'};
plot_by_label(all_data,labels,label_types)

%% should maybe actually try splitting labels by animal
%splitting labels in this way can allow us possibly to do statistics
%(do inhibitory neurons have significantly more Fac, etc.)
%somehow this seems to be resulting in differences from other labels
labels_split = labels_during_stim_split(all_data_split,error_all_split);

%plot average responses - not sure how different this is from plot_by_label
%function
cell_types = {'Excitatory','Inhibitory','Mixed'};
label_types = {'RA','SA','NA','Fac','Decrease'};
for curr_cell = 1:3
    legend_text{curr_cell} = [];
    figure 
    hold on
    for curr_label = 1:5
        temp = [];
        for curr_animal = 1:num_animals
           temp = [temp; all_data_split{curr_animal,curr_cell}(labels_split{curr_cell}(curr_animal,:)==curr_label,:)];
        end
        if ~isempty(temp)
            plot(x_label, nanmean(temp,1), 'LineWidth', 2)
            legend_text{curr_cell} = [legend_text{curr_cell}, label_types(curr_label)];
        end
    end
%     plot(x_label, nanmean(all_data{curr_cell}(labels(curr_cell,:)==6,:),1))
%     title(stim_types{curr_cell})
    title(cell_types{curr_cell})
    ylabel('dF/F0')
    xlabel('Time (s)')
    legend(legend_text{curr_cell})
end

%try to do stats with labels split
for curr_cell = 1:3
    for curr_animal = 1:4
        num_rois = sum((labels_split{curr_cell,curr_animal})~=0,2);
        for curr_label = 1:5
            label_cnt{curr_cell,curr_animal}(:,curr_label) = sum(labels_split{curr_cell,curr_animal}==curr_label,2);
            %turn it into percent
            label_cnt{curr_cell,curr_animal}(:,curr_label) = (label_cnt{curr_cell,curr_animal}(:,curr_label)./num_rois)*100;
        end
    end
end

%do stats for individual comparisons

%this is for percent facilitating for low frequency trains inh vs exc
array_label_exc = [];
array_label_inh = [];
for curr_animal = 1:4
    array_label_exc = [array_label_exc, label_cnt{1,curr_animal}(1,4), label_cnt{1,curr_animal}(2,4), label_cnt{1,curr_animal}(4,4)];
    array_label_inh = [array_label_inh, label_cnt{2,curr_animal}(1,4), label_cnt{2,curr_animal}(2,4), label_cnt{2,curr_animal}(4,4)];
end

%to check if TBS has significantly higher percent facilitating than others
for curr_animal = 1:4
    array_label_inh2(1,curr_animal) = label_cnt{1,curr_animal}(1,4);
    array_label_inh2(2,curr_animal) = label_cnt{1,curr_animal}(2,4);
    array_label_inh2(3,curr_animal) = label_cnt{1,curr_animal}(4,4);
end

%compare facilitating inhibitory neurons at 100 Hz to other populations
for curr_animal = 1:4
    array_label_comp(1,curr_animal) = [label_cnt{1,curr_animal}(3,4)];
    array_label_comp(2,curr_animal) = [label_cnt{2,curr_animal}(3,4)];
end

%compare amount of depressing neurons between excitatory and inhibitory
for curr_animal = 1:4
    array_label_comp2(1,curr_animal) = [label_cnt{1,curr_animal}(3,3)+label_cnt{1,curr_animal}(3,2)];
    array_label_comp2(2,curr_animal) = [label_cnt{2,curr_animal}(3,1)+label_cnt{2,curr_animal}(3,2)];
end

%% comparison of labels across animals for significant differences

%TBS facilitation inh vs. exc
for curr_animal = 1:4
    TBS_exc_cnt(curr_animal) = sum(labels_split{1,curr_animal}(4,:) == 4)./sum(labels_split{1,curr_animal}(4,:)~=0);
    TBS_inh_cnt(curr_animal) = sum(labels_split{2,curr_animal}(4,:) == 4)./sum(labels_split{2,curr_animal}(4,:)~=0);
end

%TBS facilitation vs other trains for inh
for curr_animal = 1:4
    TBS_inh_cnt(curr_animal) = sum(labels_split{2,curr_animal}(4,:) == 4)./sum(labels_split{2,curr_animal}(4,:)~=0);
    Hz10_inh_cnt(curr_animal) = sum(labels_split{2,curr_animal}(1,:) == 4)./sum(labels_split{2,curr_animal}(1,:)~=0);
    Hz10B_inh_cnt(curr_animal) = sum(labels_split{2,curr_animal}(2,:) == 4)./sum(labels_split{2,curr_animal}(2,:)~=0);
end
        

%% after stim
[labels2, offset_mean_all, prestim_mean_all] = labels_after_stim(all_data,error_all);
label_types2 = {'post-ICMS depression', 'post-ICMS baseline', 'post-ICMS excitation', 'post-ICMS rebound'};
plot_by_label(all_data,labels2,label_types2)

%% after stim split
[labels2_split] = labels_after_stim_split(all_data_split,error_all_split);

%% can see if labels are correlated - going to need to do something similar to what we did before
%labels are currently ordered the same way so it should be easy to look at
%correlations or percentages
%if we want to do stats will probably need to break down by animal

%can also limit to specific subpopulation and compare
%should probably only use labels that were plotted i.e. have more than 2%
%population
%this combines across cell types
clear label_cnts label_percents
for curr_train = 1:4
    for curr_label = 1:5
        idx = labels{curr_train}==curr_label; %find neurons that are identified with current during stim label
        labels_curr = labels2{curr_train}(idx); %find the post-stim labelss for the identified neurons
        for curr_label2 = 1:4
            label_cnts{curr_train}(curr_label,curr_label2) = sum(labels_curr == curr_label2, 'omitnan');
        end
        label_percents{curr_train}(curr_label,:) = (label_cnts{curr_train}(curr_label,:)/sum(label_cnts{curr_train}(curr_label,:)))*100;
    end
    %labels_percents{curr_train}= (label_cnts{curr_train}./sum(label_cnts{curr_train}))*100
    %bar(labels_percents{curr_train}')
end

for curr_train = 1:4
    %for curr_plot = 1:5
        figure
        bar(label_percents{curr_train})
        hold on
        title([train_types{curr_train}])
    %end
end

% labels_dur = [labels(1,:),labels(2,:)];
% labels_after = [labels2(1,:),labels2(2,:)];
% 
% for curr_label = 1:5
%     for curr_label2 = 1:4
%         idx = labels_dur == curr_label;
%         cnt_labels(curr_label,curr_label2) = sum(labels_after(idx)==curr_label2);
%     end
% end

%can classify neurons based on before and after and then plot percents of
%each subgroup and average response

%% code for looking at location
%this will evaluate but not sure it is correct, the amount of labels in
%each row did not look like it corresponded to the amount of rois

%went through and double checked and at least the portion looking at
%distance of individual ROIs by label should be correct

cd('C:\Users\chugh\Documents\GitHub\ExcitatoryInhibitory\Locations')
%cd('C:\Users\chugh\Desktop\Locations backup')
D = dir();
active_neur_dist = cell(3,4);
active_neur_dist_split = cell(3,4);
active_neur_vals = cell(3,4);
active_neur_vals_split = cell(3,4);
label_neur_dist = cell(3,4,5);
for curr_animal = 3:2:length(D)
    
    %EIM08 is 1x instead of 1.5 x
    if curr_animal == 3
        pixel2um = 1.61;
    else
        pixel2um = 1.075; %for 1.5x
    end
    
    %pull in electrode location
    all_elecs{ceil((curr_animal-2)/2)} = readmatrix(D(curr_animal).name).*pixel2um; 
    
    %pull in roi locations
    all_rois{ceil((curr_animal-2)/2)} = readmatrix(D(curr_animal+1).name).*pixel2um;

    %remove first column because it is just a count
    curr_roi = all_rois{ceil((curr_animal-2)/2)}(:,2:end);
    
    %calculate the distance of each roi from the electrode
    dist_from_elec = sqrt((curr_roi(:,1) - all_elecs{ceil((curr_animal-2)/2)}(2)).^2+(curr_roi(:,2) - all_elecs{ceil((curr_animal-2)/2)}(3)).^2);

    for curr_cell = 1:3  
        for curr_train = 1:4
        
            %the current index will be where there are NaN values for the
            %current animal and cell division i.e. inactive neurons
            %will also be NaN for neurons that are of different cell
            %division to maintain consistent size of arrays
            curr_idx = isnan(all_data_split{curr_train,ceil((curr_animal-2)/2),curr_cell}(:,1));
            %pull the distance values based on the index
            curr_dist = dist_from_elec;
            curr_dist(curr_idx) = NaN; %change all inactive ones or other cell types into NaN to maintain same size
            %save the distance for the given cell type across animals (will
            %stack so will lose separation between animals here)
            
            %remove data that doesn't have the switch
%             switch_idx = labels_split{curr_cell,ceil((curr_animal-2)/2)}(4,:) == 3 & labels_split{curr_cell,ceil((curr_animal-2)/2)}(2,:) == 4;
%             curr_dist(~switch_idx) = NaN;

%             offset_idx = labels2_split{curr_cell,ceil((curr_animal-2)/2)}(curr_train,:) == 4;
%             curr_dist(~offset_idx) = NaN;
            
            active_neur_dist{curr_cell,curr_train} = [active_neur_dist{curr_cell, curr_train}; curr_dist]; 
            
            active_neur_dist_split{curr_cell,curr_train}(:,ceil((curr_animal-2)/2)) = [curr_dist; NaN(3000-length(curr_dist),1)];
            
            %separate distance based on labels, need split labels for this
            curr_labels = labels_split{curr_cell,ceil((curr_animal-2)/2)}(curr_train,:);
            
            %will likely want to do this for both kinds of labels, but for
            %now only doing during ICMS
            
            %now pull the data for all active neurons
            %for this actually the all_data_split already has data set to
            %NaN if it is inactive or not the right cell type so we can
            %probably just save it directly
            curr_vals = all_data_split{curr_train,ceil((curr_animal-2)/2),curr_cell}; 
                        
            %save the values that correspond to the given distances
            active_neur_vals{curr_cell,curr_train} = [active_neur_vals{curr_cell,curr_train}; curr_vals];
            
            %split to do stats for ei ratio, must take mean across train
            active_neur_vals_split{curr_cell,curr_train}(:,ceil((curr_animal-2)/2)) = [nanmean(curr_vals,2); NaN(3000-size(curr_vals,1),1)];
            
            %go through each label and find the distance based on the
            %labels
            for curr_label = 1:5
                curr_idx = curr_labels == curr_label;
                curr_dist = dist_from_elec(curr_idx);
                %combining across animals here
                label_neur_dist{curr_cell,curr_train,curr_label} = [label_neur_dist{curr_cell,curr_train,curr_label}; curr_dist]; 
                %split by animal
                label_neur_dist_split{curr_cell,curr_train,curr_label}(ceil((curr_animal-2)/2),:) = [curr_dist; NaN(3000-length(curr_dist),1)]; 
            end
        end
    end
end

%% look at shape/size of different neuron types (no significant differences)
for curr_cell = 1:3
    area_rois{curr_cell} = [];
    for curr_animal = 1:4
        if curr_cell == 1
            curr_idx = exc_idx{curr_animal+2};
        elseif curr_cell == 2
            curr_idx = inh_idx{curr_animal+2};
        else
            curr_idx = oth_idx{curr_animal+2};
        end
        curr_area = pi*0.5*(all_rois{curr_animal}(curr_idx,4).*all_rois{curr_animal}(curr_idx,5));
        area_rois{curr_cell}  = [area_rois{curr_cell}; curr_area];
    end
    area_rois{curr_cell} = [area_rois{curr_cell}; NaN(1500-length(area_rois{curr_cell}),1)];
end

%%  violin plot for avg distance
%this will plot the average distance for each cell type
for curr_train = 1:4
    for curr_cell = 1:3
        bar_all{curr_train}(:,curr_cell) = active_neur_dist{curr_cell,curr_train};
    end
    [p_dist(curr_train),tbl_dist,stats_dist] = anova1(bar_all{curr_train});
    m_dist{curr_train} = multcompare(stats_dist);
    figure
    violin(bar_all{curr_train})
    ax = gca;
    ylim([0,max(max(bar_all{curr_train}))]);
    xlim([0,4])
    ax.XTick = [1:3];
    ax.XTickLabels = {'Exc','Inh','Mix'};
    ylabel('Mean distance')
    title(train_types{curr_train})
end

%% violin plot of just exc neurons across trains
for curr_train = 1:4
    bar_all_2(:,curr_train) = active_neur_dist{1,curr_train};
end
% [p_dist(curr_train),tbl_dist,stats_dist] = anova1(bar_all{curr_train});
% m_dist{curr_train} = multcompare(stats_dist);
figure
violin(bar_all_2)
ax = gca;
ylim([0,max(max(bar_all_2))]);
xlim([0,5])
ax.XTick = [1:4];
ax.XTickLabels = {'10 Hz','10 Hz B','100 Hz', 'TBS'};
ylabel('Mean distance')
title('Exc Distance')

%% plot for active neurons of each subtype across distance 25 um distance bins
bins = [0:25:540]; 
for curr_train = 1:4
    figure
    hold on
    for curr_cell = 1:3
        for curr_bin = 1:length(bins)-1
            %plot average waveform based on distance
            dist_idx = active_neur_dist{curr_cell,curr_train} > bins(curr_bin) & active_neur_dist{curr_cell,curr_train} <= bins(curr_bin+1);
            bin_cnt_small{curr_train}(curr_cell,curr_bin) = sum(dist_idx);
            
            %can potentially use the split labels if we want to do any
            %stats
            dist_idx_split = active_neur_dist_split{curr_cell,curr_train} > bins(curr_bin) & active_neur_dist_split{curr_cell,curr_train} <= bins(curr_bin+1);
            bin_cnt_split_small{curr_train,curr_bin}(curr_cell,:) = sum(dist_idx_split,1);
        end            
        %plot sum and divisions
        yyaxis left
        plot(bin_cnt_small{curr_train}(curr_cell,:), 'LineWidth', 2)
        ylabel('Active neurons count')
%         %ylim([0,100])
%         yyaxis right
%         cum_sum = cumsum(bin_cnt_small{curr_train}(curr_cell,:));
%         plot(cum_sum/max(cum_sum), 'LineWidth', 2)
%         ylabel('Cumulative sum')

    end
    xlim([0,22])
    legend(cell_types)
    ax = gca;
    ax.XTick = [1:length(bins)-1];
    ax.XTickLabel = {'0-25','25-50','50-75','75-100','100-125','125-150','150-175','175-200','200-225','225-250','250-275',...
        '275-300','300-325','325-350','350-375','375-400','400-425','425-450','450-475','475-500','500-525'};
    ylabel('Neuron count')
    title(train_types{curr_train})
    xline(5)
    xline(8)
    
    %ei_ratio{curr_train} = all_vals_split{curr_train,1}./all_vals_split{curr_train,2};
    %ei_ratio{curr_train} = bin_cnt{curr_train}(1,:)./bin_cnt{curr_train}(2,:);
end

%can create a bar plot for each big bin
bin_cnt_big(:,3) = sum(bin_cnt_small{1,3}(:,9:end),2);
bin_cnt_big(:,2) = sum(bin_cnt_small{1,3}(:,6:8),2);
bin_cnt_big(:,1) = sum(bin_cnt_small{1,3}(:,1:5),2);

%if we want stats - need to instead break into the bigger bins and keep
%animals separate
for curr_cnt = 1:21
    big_bin_means(:,curr_cnt) = nanmean(bin_cnt_split_small{3,curr_cnt},2);
end
%% can also look at distance bins and how exc-inh ratio varies across distance
%can look at both number active and intensity potentially
%if we want to make bar plots, will need to use active_neur_dist_split
clear all_vals all_vals_split ei_ratio bin_cnt bin_cnt_split mdl
bins = [0,125,200,540]; 
for curr_train = 1:4
    for curr_cell = 1:3
        for curr_bin = 1:3
            %plot average waveform based on distance
            dist_idx = active_neur_dist{curr_cell,curr_train} > bins(curr_bin) & active_neur_dist{curr_cell,curr_train} <= bins(curr_bin+1);
            bin_cnt{curr_train}(curr_cell,curr_bin) = sum(dist_idx);
            
            %split for bar plots
            dist_idx_split = active_neur_dist_split{curr_cell,curr_train} > bins(curr_bin) & active_neur_dist_split{curr_cell,curr_train} <= bins(curr_bin+1);
            bin_cnt_split{curr_train,curr_bin}(curr_cell,:) = sum(dist_idx_split,1)+1;
            %make minimum count 1 for ratio calculation
            
            %calculate mean intensity for all neurons of given type
            all_vals{curr_train,curr_cell}(:,curr_bin) = [nanmean(active_neur_vals{curr_cell,curr_train}(dist_idx,:),2); NaN(3000-sum(dist_idx),1)]; 
            
            %calculte mean intensity split by animal
            temp = active_neur_vals_split{curr_cell,curr_train};
            temp(~dist_idx_split) = NaN;
            all_vals_split{curr_train,curr_bin}(curr_cell,:) = nanmean(temp,1);
        end
        if ~any(all(isnan(all_vals{curr_train,curr_cell}(:,1:3))))
            [p_int(curr_train,curr_cell),tbl_dist,stats_int] = anova1(all_vals{curr_train,curr_cell});
            m_int{curr_train,curr_cell} = multcompare(stats_int);
        
            figure
            %scatter([1:3],nanmean(all_vals{curr_train,curr_cell}(:,1:3),1))
            violin(all_vals{curr_train,curr_cell}(:,1:3),1)
            ylim([-0.4,2])
            hold on

            mdl{curr_train,curr_cell} = fitlm([1:3],nanmean(all_vals{curr_train,curr_cell}(:,1:3),1));
            plot([1:3], mdl{curr_train,curr_cell}.Fitted, 'k', 'LineWidth', 2) 

            %plot cumulative sum and divisions - can do it with error bars if
            %we use split
            yyaxis right
            plot(bin_cnt{curr_train}(curr_cell,[1:3]), 'LineWidth', 2)
            ylabel('Active neurons count')
            ylim([0,100])

            xlim([0,4])
            ax = gca;
            ax.XTick = [1:3];
            ax.XTickLabel = {'0-125','125-200','200+'};
            ylabel('dF/F0')
            title([train_types{curr_train}, ' ', cell_types{curr_cell}])
        end
    end
    %ei_ratio{curr_train} = all_vals_split{curr_train,1}./all_vals_split{curr_train,2};
    %ei_ratio{curr_train} = bin_cnt{curr_train}(1,:)./bin_cnt{curr_train}(2,:);
end
figure
bar(bin_cnt{3}(1,:))

%instead make a bar plot for percents based on reviewer comments
%hard code percents because otherwise we have to change code to have both
%active and inactive at the same time
perc_exc = [75, 84, 42];
figure, bar(perc_exc)

perc_inh = [100, 92, 67];
figure, bar(perc_inh)


figure
bar_mean = (nanmean(bin_cnt_split{3,1},2));
bar(bar_mean)
hold on
bar_error = nanstd(bin_cnt_split{3,1},[],2)./sqrt(4);
e = errorbar([1:3], bar_mean, bar_error);
e.LineStyle = 'none';

figure
bar_mean = bin_cnt{3}(1,:);
bar(bar_mean)
hold on
bar_error = nanstd(bin_cnt_split{3,1},[],2)./sqrt(4);
e = errorbar([1:3], bar_mean, bar_error);
e.LineStyle = 'none';

%% plot ei ratio based on count and split by animal
figure
for curr_bin = 1:3
    bar_ei_ratio(curr_bin,:) = (bin_cnt_split{3,curr_bin}(1,:))./(bin_cnt_split{3,curr_bin}(2,:));
    bar_ei_mean(curr_bin) = nanmean(bin_cnt_split{3,curr_bin}(1,:))./nanmean(bin_cnt_split{3,curr_bin}(2,:));
    curr_error = bin_cnt_split{3,curr_bin}(1,:)./(bin_cnt_split{3,curr_bin}(2,:));
    bar_ei_error(curr_bin) = nanstd(curr_error(isfinite(curr_error)))/sqrt(4);
end
bar(bar_ei_mean)
hold on
e = errorbar([1:3], bar_ei_mean, bar_ei_error);
e.LineStyle = 'none';
ax = gca;
ax.XTickLabel = {'0-125','125-200','200+'};
ylabel('EI Ratio')
xlabel('Distance from electrode (\mum)')

%ei ratio based on intensity
figure
for curr_bin = 1:3
    bar_ei_ratio(curr_bin,:) = (all_vals_split{3,curr_bin}(1,:))./(all_vals_split{3,curr_bin}(2,:));
    bar_ei_mean(curr_bin) = nanmean(all_vals_split{3,curr_bin}(1,:))./nanmean(all_vals_split{3,curr_bin}(2,:));
    bar_ei_error(curr_bin) = nanstd(bar_ei_ratio(curr_bin,(isfinite(bar_ei_ratio(curr_bin,:)))))/sqrt(4);
end
bar(bar_ei_mean)
hold on
e = errorbar([1:3], bar_ei_mean, bar_ei_error);
e.LineStyle = 'none';
ax = gca;
ax.XTickLabel = {'0-125','125-200','200+'};
ylabel('EI Ratio')
xlabel('Distance from electrode (\mum)')

%ei ratio based on count and intensity
figure
for curr_bin = 1:3
    bar_ei_ratio(curr_bin,:) = (bin_cnt_split{3,curr_bin}(1,:).*all_vals_split{3,curr_bin}(1,:))./(bin_cnt_split{3,curr_bin}(2,:).*all_vals_split{3,curr_bin}(2,:));
    bar_ei_mean(curr_bin) = (nanmean(bin_cnt_split{3,curr_bin}(1,:))*nanmean(all_vals_split{3,curr_bin}(1,:)))./(nanmean(bin_cnt_split{3,curr_bin}(2,:))*nanmean(all_vals_split{3,curr_bin}(2,:)));
    bar_ei_error(curr_bin) = nanstd(bar_ei_ratio(curr_bin,(isfinite(bar_ei_ratio(curr_bin,:)))))/sqrt(4);
end
bar(bar_ei_mean)
hold on
e = errorbar([1:3], bar_ei_mean, bar_ei_error);
e.LineStyle = 'none';
ax = gca;
ax.XTickLabel = {'0-125','125-200','200+'};
ylabel('EI Ratio')
xlabel('Distance from electrode (\mum)')


%plot ei ratio based on count and intensity
figure
bar_ei_ratio = (bin_cnt{1,3}(1,:).*nanmean(all_vals{3,1},1))./(bin_cnt{1,3}(2,:).*nanmean(all_vals{3,2},1));
bar(bar_ei_ratio)
ax = gca;
ax.XTickLabel = {'0-125','125-200','200+'};
ylabel('EI Ratio')
xlabel('Distance from electrode (\mum)')

%% plot average trace within each bin
bins = [0,125,200,540]; 
for curr_train = 1:4
    for curr_cell = 1:3
        figure
        for curr_bin = 1:length(bins)-1
            %figure
            curr_resp = NaN(1,1800);
            for curr_animal = 1:4
                dist_idx_split = active_neur_dist_split{curr_cell,curr_train}(:,curr_animal) > bins(curr_bin) & active_neur_dist_split{curr_cell,curr_train}(:,curr_animal) <= bins(curr_bin+1);
                curr_resp = [curr_resp;all_data_split{curr_train,curr_animal,curr_cell}(dist_idx_split,:)];
            end
%             for i = 1:6
%                 mean_resp = nanmean(curr_resp,1);
%                 grouped_resps(i,:) = mean_resp(994+(31*(i-1)):994+(31*i));
%             end
            plot(x_label,nanmean(curr_resp,1),'LineWidth',2)
%             plot([1:32],nanmean(grouped_resps,1),'LineWidth',2)
            hold on
        end
        title([cell_types{curr_cell}, ' ', train_types{curr_train}])
        legend({'Within 50 um', '50-200 um', '200+ um'})
    end
end



%% distance by labels
cell_types = {'Excitatory','Inhibitory','Mixed'};
train_types = {'10 Hz', '10 Hz Burst', '100 Hz', 'TBS'};
for curr_cell = 1:3
    for curr_train = 1:4
        label_names{curr_cell,curr_train} = {'RA','SA','NA','Fac','Dec'};
        num_rois = 0;
        for curr_animal = 1:4
            num_rois = num_rois+sum(sum(labels_split{curr_cell,curr_animal}(curr_train,:)~=0));
        end
        log_idx = zeros(1,5);
        for curr_label = 1:5
            label_cnt = 0;
            for curr_animal = 1:4
                label_cnt = label_cnt + sum(labels_split{curr_cell,curr_animal}(curr_train,:)==curr_label);
            end
            if label_cnt >= ceil(0.02*num_rois)
                bar_all_2{curr_cell,curr_train}(:,curr_label) = [label_neur_dist{curr_cell,curr_train,curr_label}; NaN(3000-length(label_neur_dist{curr_cell,curr_train,curr_label}),1)];
                log_idx(curr_label) = 1;
            else
                bar_all_2{curr_cell,curr_train}(:,curr_label) = NaN(3000,1);
            end
        end
        [p4(curr_cell,curr_train),tbl,stats] = anova1(bar_all_2{curr_cell,curr_train});
        m4{curr_cell,curr_train} = multcompare(stats);
        figure 
        violin(bar_all_2{curr_cell,curr_train}(:,any(~isnan(bar_all_2{curr_cell,curr_train}))))
        y_lim = ylim;
        y_lim(1) = 0;
        ylim(y_lim)
        ax = gca;
        label_names{curr_cell,curr_train}(~log_idx) = [];
        ax.XTick = [1:length(label_names{curr_cell,curr_train})];
        ax.XTickLabel = label_names{curr_cell,curr_train};
        ylabel('Distance (\mum)')
        title([cell_types{curr_cell}, ' ', train_types{curr_train}])
    end
end

%% distance by bins - output not looking right for this should double check
clear labels_percents labels_cnts
bins = [0,125,200,540];
for curr_train = 1:4
    for curr_cell = 1:3
        for curr_bin = 1:3
    %         figure
            for curr_label = 1:5
                %plot bar plot for split of classes
                curr_num(curr_label) = sum(label_neur_dist{curr_cell,curr_train,curr_label} > bins(curr_bin) & label_neur_dist{curr_cell,curr_train,curr_label} <= bins(curr_bin+1));
                curr_num_split =  sum(label_neur_dist_split{curr_cell,curr_train,curr_label} > bins(curr_bin) & label_neur_dist_split{curr_cell,curr_train,curr_label} <= bins(curr_bin+1),2);
                labels_cnts_split{curr_train,curr_bin,curr_label}(curr_cell,:) = curr_num_split; %will be split by animal and label
        %         dist_int(curr_train,curr_bin) = nanmean(active_neur_mean{curr_train}(dist_idx));
        %         dist_int_num(curr_train,curr_bin) = sum(dist_idx);
            end
            labels_cnts{curr_train,curr_cell}(curr_bin,:) = curr_num;
            labels_percents{curr_train,curr_cell}(curr_bin,:) = (curr_num/sum(curr_num))*100;
    %         bar(curr_perc)
    %         title([stim_types{curr_train}, ' ', num2str(bins(curr_bin+1)), ' um'])
        end
        %legend('0-50','50-100','100-200','200-300','300-400','400+')
    end
end

%can calculate number of neurons from each class that make up facilitating
%neurons in each bin

%split version - not statistical significance
% for curr_cell = 1:3
%     for curr_bin = 1:3
%         %for plot
%         bin_cnts_fac{curr_cell}(curr_bin,:) = labels_cnts_split{3,curr_bin,4}(curr_cell,:);
%         %for stats
%         bin_cnts_fac2{curr_bin}(curr_cell,:) = labels_cnts_split{3,curr_bin,4}(curr_cell,:);
%     end
% end
% figure
% plot(nanmean(bin_cnts_fac{1},2))
% hold on
% errorbar([1:3], nanmean(bin_cnts_fac{1},2), nanstd(bin_cnts_fac{1},[],2)./sqrt(4))
% plot(nanmean(bin_cnts_fac{2},2))
% errorbar([1:3], nanmean(bin_cnts_fac{2},2), nanstd(bin_cnts_fac{2},[],2)./sqrt(4))
% % plot(nanmean(bin_cnts_fac{3},2))
% % errorbar([1:3], nanmean(bin_cnts_fac{3},2), nanstd(bin_cnts_fac{3},[],2)./sqrt(4))
% xlim([0,4])
% ylim([0,15])
% ylabel('Neuron count')
% ax = gca
% ax.XTick = [1:3];
% ax.XTickLabel = {'0-100','100-200','200+'};
% legend('Exc','Inh','Mix')
% title('Number of Fac neurons')

%can calculate number of neurons from each class that make up facilitating
%neurons in each bin
for curr_cell = 1:3
    bin_cnts_fac(:,curr_cell) = labels_cnts{3,curr_cell}(:,4);
end
figure
plot(bin_cnts_fac(:,1))
hold on
plot(bin_cnts_fac(:,2))
plot(bin_cnts_fac(:,3))
xlim([0,4])
ylim([0,15])
ylabel('Neuron count')
ax = gca
ax.XTick = [1:3];
ax.XTickLabel = {'0-100','100-200','200+'};
legend('Exc','Inh','Mix')
title('Number of Fac neurons')

for curr_cell = 1:3
    bin_cnts_ra(:,curr_cell) = labels_cnts{3,curr_cell}(:,1);
end
figure
plot(bin_cnts_ra(:,1))
hold on
plot(bin_cnts_ra(:,2))
plot(bin_cnts_ra(:,3))
xlim([0,4])
ylim([0,15])
title('Number of RA neurons')
ylabel('Neuron count')
ax = gca
ax.XTick = [1:3];
ax.XTickLabel = {'0-100','100-200','200+'};
legend('Exc','Inh','Mix')

for curr_cell = 1:3
    bin_cnts_sa(:,curr_cell) = labels_cnts{3,curr_cell}(:,2);
end
figure
plot(bin_cnts_sa(:,1))
hold on
plot(bin_cnts_sa(:,2))
plot(bin_cnts_sa(:,3))
xlim([0,4])
%ylim([0,15])
title('Number of SA neurons')
ylabel('Neuron count')
ax = gca
ax.XTick = [1:3];
ax.XTickLabel = {'0-100','100-200','200+'};
legend('Exc','Inh','Mix')

for curr_cell = 1:3
    bin_cnts_na(:,curr_cell) = labels_cnts{3,curr_cell}(:,3);
end
figure
plot(bin_cnts_na(:,1))
hold on
plot(bin_cnts_na(:,2))
plot(bin_cnts_na(:,3))
xlim([0,4])
%ylim([0,15])
title('Number of NA neurons')
ylabel('Neuron count')
ax = gca
ax.XTick = [1:3];
ax.XTickLabel = {'0-125','125-200','200+'};
legend('Exc','Inh','Mix')
        

% cd('C:\Users\chugh\Documents\GitHub\BiomimeticMice\Figures\Rev Fig 2\Long_for_short')
% for curr_train = 1:4
%     for curr_plot = 1:5
%         figure
%         bar(labels_percents{curr_train}(:,curr_plot))
%         hold on
%         percent_labels = (sum(labels(curr_train,:)==curr_plot))/(sum(~isnan(labels(curr_train,:))))*100; 
%         yline(percent_labels)
%         title([stim_types{curr_train},' Label ', num2str(curr_plot)])
%         ax = gca;
%         ax.XTickLabels = {'0-50','50-100','100-200','200-300','300-400','400+'};
% %         filename = [stim_types{curr_train}, '_label_', num2str(curr_plot)];
% %         print(filename,'-dsvg');
%     end
% end


%% also can divide into bins and look at the percent in each bin
%can only really do stats here if we divide by animal, may want to divide
%by animal earlier instead of using all_data and this will also make it
%easier to do labels by distance probably
bins = [0,50,100,200,300,400,540];
for curr_cell = 1:2
    for curr_bin = 1:6
        %plot average waveform based on distance
        dist_idx = active_neur_dist{curr_cell} > bins(curr_bin) & active_neur_dist{curr_cell} <= bins(curr_bin+1);
        cell_cnt(curr_cell,curr_bin) = sum(dist_idx);
    end
end
perc_inh = cell_cnt(2,:)./(cell_cnt(1,:)+cell_cnt(2,:));

figure
bar(perc_inh)

%% calculate adaptation indices of all neurons
%this is messed up, putting a bunch of zeros
%maybe we should be weighting by the overall intensity
%look into this more
start_range = ceil(onset_start:onset_start+(offset_stop-onset_start)/2);
stop_range = ceil(onset_start+(offset_stop-onset_start)/2:offset_stop);
for curr_train = 1:4
    adapt_idx{curr_train} = NaN(3,3000);
    for curr_cell = 1:3
        for curr_neuron = 1:size(all_data{curr_train,curr_cell},1)
            curr_resp = all_data{curr_train,curr_cell}(curr_neuron,:);
            curr_int = nanmean(curr_resp);
            adapt_idx{curr_train}(curr_cell,curr_neuron) = abs(mean(curr_resp(stop_range))) - abs(mean(curr_resp(start_range)));
            %adapt_idx{curr_train}(curr_cell,curr_neuron) = (abs(mean(curr_resp(start_range))) - abs(mean(curr_resp(stop_range))))/(abs(mean(curr_resp(start_range))) + abs(mean(curr_resp(stop_range))));
            %adapt_idx{curr_train}(curr_cell,curr_neuron) = adapt_idx{curr_train}(curr_cell,curr_neuron)*curr_int;
            %if we want to calculate "percent"
            %adapt_idx{curr_train}(curr_cell,curr_neuron) = (abs(max(curr_resp(start_range))) - abs(max(curr_resp(stop_range))))/(abs(max(curr_resp(start_range))));
        end
        sem_adapt{curr_train}(curr_cell) = nanstd(adapt_idx{curr_train}(curr_cell,:))/sqrt(sum(~isnan(adapt_idx{curr_train}(curr_cell,:))));
    end
    mean_idx = nanmean(adapt_idx{curr_train},2);
    sem_idx = nanstd(adapt_idx{curr_train},[],2)./(sqrt(sum(~isnan(adapt_idx{curr_train}),2)));

    figure
    violin(adapt_idx{curr_train}')
    %swarmchart([1:3], adapt_idx{curr_train}','XJitter','rand')
    ax = gca;
    ax.XTick = [1:3];
    ax.XTickLabels = {'Exc','Inh','Mix'};
    title(train_types(curr_train))
    ylabel('dF/F0')
end

% for curr_train = 1:4
%     [pa(curr_train),tbl_a{curr_train},stats_a{curr_train}] = anova1(adapt_idx{curr_train});
%     mpa{curr_train} = multcompare(stats_a{curr_train});
% end

%can also calculate stats across trains within cell types
for curr_cell = 1:3
    for curr_train = 1:4
        adapt_idx2{curr_cell}(:,curr_train) = adapt_idx{curr_train}(curr_cell,:)';
    end
    [p_trains(curr_cell),~,stats_trains{curr_cell}] = anova1(adapt_idx2{curr_cell});
    m_trains{curr_cell} = multcompare(stats_trains{curr_cell});
end

%% calculate post-stim adaptation
for curr_train = 1:4
    adapt_idx_post{curr_train} = NaN(3,3000);
    for curr_cell = 1:3
        for curr_neuron = 1:size(all_data{curr_train,curr_cell},1)
            curr_resp = all_data{curr_train,curr_cell}(curr_neuron,:);
            %curr_int = nanmean(curr_resp);
            adapt_idx_post{curr_train}(curr_cell,curr_neuron) = min(curr_resp(stop_range:end)) - mean(curr_resp(1:start_range));
            %adapt_idx{curr_train}(curr_cell,curr_neuron) = (abs(mean(curr_resp(start_range))) - abs(mean(curr_resp(stop_range))))/(abs(mean(curr_resp(start_range))) + abs(mean(curr_resp(stop_range))));
            %adapt_idx{curr_train}(curr_cell,curr_neuron) = adapt_idx{curr_train}(curr_cell,curr_neuron)*curr_int;
            %if we want to calculate "percent"
            %adapt_idx{curr_train}(curr_cell,curr_neuron) = (abs(max(curr_resp(start_range))) - abs(max(curr_resp(stop_range))))/(abs(max(curr_resp(start_range))));
        end
    end

    figure
    violin(adapt_idx_post{curr_train}')
    %swarmchart([1:3], adapt_idx{curr_train}','XJitter','rand')
    ax = gca;
    ax.XTick = [1:3];
    ax.XTickLabels = {'Exc','Inh','Mix'};
    title(train_types(curr_train))
    ylabel('dF/F0')
end

%plot them against each other
all_adapt = [];
all_adapt_post = [];
for curr_train = 1:4
    all_adapt = [all_adapt; reshape(adapt_idx{curr_train},9000,1)];
    all_adapt_post = [all_adapt_post; reshape(adapt_idx_post{curr_train},9000,1)];
end

figure 
hold on
mdl = fitlm(all_adapt,all_adapt_post);
scatter(all_adapt,all_adapt_post, 'b')
plot(all_adapt, mdl.Fitted, 'k')
xlim([-0.6,1])
ylim([-0.6,0.4])
xlabel('Depression index pre-ICMS')
ylabel('Depression index post-ICMS')

%% can make code for labeling neurons 

%first load in image
img = imread('C:\Users\chugh\OneDrive\Documents\GitHub\ExcitatoryInhibitory\EIM09_100Hz_Average.png');
figure
imshow(img)
hold on

%then pull in label and distance information
pixel2um = 1.062;
colors = [168,203,254;90,156,254;13,110,253;9,75,172;5,40,91];
colors = colors/255;
locs = all_rois{2};
locs(:,2:3) = locs(:,2:3)./1.07;
locs(:,4:5) = locs(:,4:5)./1.5; %have to make a slight adjustment for conversion

%base neuron is neuron 48 from animal EIM09
curr_roi = 48;
curr_color = 1;
ellipse(locs(curr_roi,4),locs(curr_roi,5),locs(curr_roi,6),locs(curr_roi,2),locs(curr_roi,3),colors(curr_color,:))
text(locs(curr_roi,2)-4,locs(curr_roi,3)-1,'B', 'Color', colors(curr_color,:))

curr_roi = 93;
curr_color = 2;
ellipse(locs(curr_roi,4),locs(curr_roi,5),locs(curr_roi,6),locs(curr_roi,2),locs(curr_roi,3),colors(curr_color,:))

curr_roi = 91;
curr_color = 3;
ellipse(locs(curr_roi,4),locs(curr_roi,5),locs(curr_roi,6),locs(curr_roi,2),locs(curr_roi,3),colors(curr_color,:))

curr_roi = 117;
curr_color = 4;
ellipse(locs(curr_roi,4),locs(curr_roi,5),locs(curr_roi,6),locs(curr_roi,2),locs(curr_roi,3),colors(curr_color,:))

curr_roi = 127;
curr_color = 5;
ellipse(locs(curr_roi,4),locs(curr_roi,5),locs(curr_roi,6),locs(curr_roi,2),locs(curr_roi,3),colors(curr_color,:))
