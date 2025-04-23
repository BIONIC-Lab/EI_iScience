%first combine acrross animals
comb_act_idx(1,:) = [act_idx_save{1}(3,:),act_idx_save{2}(3,:),act_idx_save{3}(3,:),act_idx_save{4}(3,:)];
comb_act_idx(2,:) = [act_idx_save{1}(1,:),act_idx_save{2}(1,:),act_idx_save{3}(1,:),act_idx_save{4}(1,:)];
comb_act_idx(3,:) = [act_idx_save{1}(2,:),act_idx_save{2}(2,:),act_idx_save{3}(2,:),act_idx_save{4}(2,:)];
comb_act_idx(4,:) = [act_idx_save{1}(4,:),act_idx_save{2}(4,:),act_idx_save{3}(4,:),act_idx_save{4}(4,:)];


figure
hold on
colors = hsv(3);
for curr_train = 1:3
    for curr_label_1 = 1:3
        curr_idx_1 = comb_act_idx(curr_train,:) == curr_label_1;
        for curr_label_2 = 1:3
            curr_idx_2{curr_train}(curr_label_1,curr_label_2) = sum(comb_act_idx(curr_train+1,:) == curr_label_2 & curr_idx_1);
            curr_color = colors(curr_label_1,:);
            if curr_idx_2{curr_train}(curr_label_1,curr_label_2) > 0
                plot([curr_train,curr_train+1],[curr_label_1,curr_label_2],'LineWidth',curr_idx_2{curr_train}(curr_label_1,curr_label_2)/5,'Color',curr_color)
                num_save{curr_train}(curr_label_1,curr_label_2) = curr_idx_2{curr_train}(curr_label_1,curr_label_2);
            end
        end
    end
end

ax = gca;
ax.YTick = [1:3];
ax.YTickLabel = ['Active','Non-cons','Inactive'];
ax.YTickLabel = {'Active','Non-cons','Inactive'};
ax.XTick = [1:4];
ax.XTickLabel = {'100 Hz', '10 Hz', '10 Hz B', 'TBS'}; %{'10 Hz', '10 Hz B', '100 Hz', 'TBS'};

%% similar to above plot but plot each train to other trains separately
%first combine acrross animals
colors = hsv(3);
train_types =  {'100 Hz', '10 Hz', '10 Hz B', 'TBS'};
for curr_train_1 = 1:4
    for curr_train_2 = 1:3
        figure
        hold on
        for curr_label_1 = 1:3
            curr_idx_1 = comb_act_idx(curr_train_1,:) == curr_label_1;
            comb_act_idx_2 = comb_act_idx;
            comb_act_idx_2(curr_train_1,:) = [];
            types_idx = ones(1,4);
            types_idx(curr_train_1) = 0;
            train_types_2 = train_types(logical(types_idx));
            for curr_label_2 = 1:3
                curr_idx_2 = comb_act_idx_2(curr_train_2,:) == curr_label_2 & curr_idx_1;
                curr_color = colors(curr_label_1,:);
                if sum(curr_idx_2) > 0
                    plot([1,2],[curr_label_1,curr_label_2],'LineWidth',sum(curr_idx_2)/3,'Color',curr_color)
                end
            end
        end
        ax = gca;
        ax.YTick = [1:3];
        ax.YTickLabel = {'Active','Non-cons','Inactive'};
        ax.XTick = [1:2];
        ax.XTickLabel = [train_types(curr_train_1), train_types_2(curr_train_2)];
    end
end


%% add up percent for each combination
count_idx = zeros(3,3,3,3);
total_neurons = 0;

poss_titles{1} = 'A';
poss_titles{2} = 'N';
poss_titles{3} = 'I';

for curr_animal = 1:4
    for curr_neuron = 1:size(act_idx_save{curr_animal},2)
        curr_idx = act_idx_save{curr_animal}(:,curr_neuron);
        if all(~isnan(curr_idx))
            count_idx(curr_idx(1),curr_idx(2),curr_idx(3),curr_idx(4)) = count_idx(curr_idx(1),curr_idx(2),curr_idx(3),curr_idx(4))+1;
            total_neurons = total_neurons+1;
        end
    end
end

count_idx_2 = reshape(count_idx,[1,81])./total_neurons;

%make ylabels to label groups
for idx1 = 1:3
    for idx2 = 1:3
        for idx3 = 1:3
            for idx4 = 1:3
                titles{idx1,idx2,idx3,idx4} = ['10Hz ', poss_titles{idx1}, ' 10HzB ', poss_titles{idx2}, ' 100Hz ', poss_titles{idx3}, ' TBS ', poss_titles{idx4}];
            end
        end
    end
end
all_titles = reshape(titles,[1,81]);

%reorder by percent accounted for and reorder titles in same way
[count_idx_3,sort_idx] = sort(count_idx_2);
all_titles_2 = all_titles(sort_idx);

figure
barh(count_idx_3)
ax = gca;
ylim([36,82])
ax.YTick = [37:81];
ax.YTickLabels = all_titles_2(37:81);
xlim([0,0.8])