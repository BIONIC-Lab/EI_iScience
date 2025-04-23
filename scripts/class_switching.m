des_neur = 1; %1 for exc, 2 for inh
order = [3,2,1,4]; %order of trains for plot

idx = ~isnan(labels{1}(des_neur,:)) & ~isnan(labels{2}(des_neur,:)) & ~isnan(labels{3}(des_neur,:)) & ~isnan(labels{4}(des_neur,:));

figure
hold on
colors = hsv(5);
for curr_neuron = 1:length(idx)
    %if idx(curr_neuron)==1
    %iterates through each neuron of a specific type
    %checks if the curr neuron is a nan i.e. is active or not
    curr_idx = [~isnan(labels{order(1)}(des_neur, curr_neuron)), ~isnan(labels{order(2)}(des_neur, curr_neuron)), ~isnan(labels{order(3)}(des_neur, curr_neuron)), ~isnan(labels{order(4)}(des_neur, curr_neuron))];
    %checks if the curr neuron has a label of 4 for any trains
    %curr_idx_2 = [(labels{order(1)}(des_neur, curr_neuron))==4, (labels{order(2)}(des_neur, curr_neuron))==4, (labels{order(3)}(des_neur, curr_neuron))==4, (labels{order(4)}(des_neur, curr_neuron))==4];
    if any(curr_idx)
        first_idx = find(curr_idx==1, 1, 'first');
        curr_color = colors(labels{order(first_idx)}(des_neur,curr_neuron),:); %color by first label for 10 Hz train
        xax = [1,2,3,4];
        xax(~curr_idx) = NaN;
        curr_rand = (0.25-rand(1)*0.5);
        plot(xax,[labels{order(1)}(des_neur,curr_neuron), labels{order(2)}(des_neur,curr_neuron), labels{order(3)}(des_neur,curr_neuron), labels{order(4)}(des_neur,curr_neuron)]+curr_rand, 'Color', curr_color)
        scatter(xax,[labels{order(1)}(des_neur,curr_neuron), labels{order(2)}(des_neur,curr_neuron), labels{order(3)}(des_neur,curr_neuron), labels{order(4)}(des_neur,curr_neuron)]+curr_rand, 'MarkerFaceColor', curr_color, 'MarkerEdgeColor', 'k')
    end
    %end
end

ax = gca;
ax.YTick = [1:5];
ax.YTickLabel = {'RA', 'SA', 'NA' 'Fac', 'Dec'};
% ax.YTickLabel = {'post-ICMS depression', 'post-ICMS baseline', 'post-ICMS excitation', 'post-ICMS rebound'};
xlim([0,5])
ax.XTick = [1:4];
ax.XTickLabel = {'100 Hz', '10 Hz', '10 Hz Burst' 'TBS'};
%ax.XTickLabel = {'TBS', '10 Hz Burst', '10 Hz' '100 Hz'};

%% same logic as above but with thick lines that represent total neurons
des_neur = 2; %1 for exc, 2 for inh
order = [3,1,2,4]; %order of trains for plot

idx = ~isnan(labels{1}(des_neur,:)) & ~isnan(labels{2}(des_neur,:)) & ~isnan(labels{3}(des_neur,:)) & ~isnan(labels{4}(des_neur,:));

figure
hold on
colors = hsv(5);
for curr_train = 1:3
    for curr_label_1 = 1:5
        curr_idx_1 = labels{order(curr_train)}(des_neur,:) == curr_label_1;
        for curr_label_2 = 1:5
            curr_idx_2{curr_train}(curr_label_1,curr_label_2) = sum(labels{order(curr_train+1)}(des_neur,:) == curr_label_2 & curr_idx_1);
            curr_color = colors(curr_label_1,:);
            if curr_idx_2{curr_train}(curr_label_1,curr_label_2) > 0
                plot([curr_train,curr_train+1],[curr_label_1,curr_label_2],'LineWidth',curr_idx_2{curr_train}(curr_label_1,curr_label_2)/20,'Color',curr_color)
                num_save{curr_train}(curr_label_1,curr_label_2) = curr_idx_2{curr_train}(curr_label_1,curr_label_2);
            end
        end
    end
end

ax = gca;
ax.YTick = [1:5];
ax.YTickLabel = {'RA', 'SA', 'NA' 'Fac', 'Dec'};
ax.XTick = [1:4];
ax.XTickLabel = {'100 Hz', '10 Hz', '10 Hz B', 'TBS'}; %{'10 Hz', '10 Hz B', '100 Hz', 'TBS'};

%% can try to make pie plots that show similar thing
des_neur = 1; %1 for exc, 2 for inh

idx = ~isnan(labels{1}(des_neur,:)) & ~isnan(labels{2}(des_neur,:)) & ~isnan(labels{3}(des_neur,:)) & ~isnan(labels{4}(des_neur,:));


for curr_train = 1:4
    label_dist{curr_train} = zeros(5,5);
    for curr_label_1 = 1:5
        curr_idx_1 = labels{curr_train}(des_neur,:) == curr_label_1;
        for curr_train_2 = 1:4
            if curr_train_2 ~= curr_train
                for curr_label_2 = 1:5
                    label_dist{curr_train}(curr_label_2,curr_label_1) = label_dist{curr_train}(curr_label_2,curr_label_1) + sum(labels{curr_train_2}(des_neur,curr_idx_1) == curr_label_2);
                end
            end
        end
    end
end



%% %load in labels to limit to ones that switch
% cd('C:\Users\chugh\OneDrive\Documents\GitHub\ExcitatoryInhibitory')
% load('labels.mat')
% %ones that are fac for TBS but NA for 10 HzB
% switch_idx = labels{2} == 4 & labels{4} == 3;

%% add up percent for each combination
des_neur = 1;
count_idx = zeros(5,5,5,5);
total_neurons = 0;

poss_titles{1} = 'RA';
poss_titles{2} = 'SA';
poss_titles{3} = 'NA';
poss_titles{4} = 'Fac';
poss_titles{5} = 'Dec';

labels_use = [labels{1}(des_neur,:); labels{2}(des_neur,:); labels{3}(des_neur,:); labels{4}(des_neur,:)];


for curr_neuron = 1:size(labels_use ,2)
    curr_idx = labels_use(:,curr_neuron);
    if all(~isnan(curr_idx))
        count_idx(curr_idx(1),curr_idx(2),curr_idx(3),curr_idx(4)) = count_idx(curr_idx(1),curr_idx(2),curr_idx(3),curr_idx(4))+1;
        total_neurons = total_neurons+1;
    end
end


count_idx_2 = reshape(count_idx,[1,625])./total_neurons;

%make ylabels to label groups
for idx1 = 1:5
    for idx2 = 1:5
        for idx3 = 1:5
            for idx4 = 1:5
                titles{idx1,idx2,idx3,idx4} = ['10 Hz ', poss_titles{idx1}, ' 10 Hz B ', poss_titles{idx2}, ' 100 Hz ', poss_titles{idx3}, ' TBS ', poss_titles{idx4}];
            end
        end
    end
end
all_titles = reshape(titles,[1,625]);

%reorder by percent accounted for and reorder titles in same way
[count_idx_3,sort_idx] = sort(count_idx_2);
all_titles_2 = all_titles(sort_idx);

figure
barh(count_idx_3)
ax = gca;
ylim([599,626])
ax.YTick = [600:625];
ax.YTickLabels = all_titles_2(600:625);