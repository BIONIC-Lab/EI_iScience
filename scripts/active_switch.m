figure
hold on
colors = hsv(3);
for curr_animal = 1:4
    for curr_neuron = 1:size(act_idx_save{1,curr_animal},2)
        curr_rand = (0.25-rand(1)*0.5);
        curr_color = colors(act_idx_save{1,curr_animal}(3,curr_neuron),:);
        plot([1:4],[act_idx_save{1,curr_animal}(3,curr_neuron),act_idx_save{1,curr_animal}(1,curr_neuron),act_idx_save{1,curr_animal}(2,curr_neuron),act_idx_save{1,curr_animal}(4,curr_neuron)]+curr_rand,'Color',curr_color)
        scatter([1:4],[act_idx_save{1,curr_animal}(3,curr_neuron),act_idx_save{1,curr_animal}(1,curr_neuron),act_idx_save{1,curr_animal}(2,curr_neuron),act_idx_save{1,curr_animal}(4,curr_neuron)]+curr_rand,'MarkerFaceColor',curr_color,'MarkerEdgeColor',curr_color)
    end
end
ax = gca;
ax.YTick = [1:3];
ax.YTickLabel = ['Active','Non-cons','Inactive'];
ax.YTickLabel = {'Active','Non-cons','Inactive'};
ax.XTick = [1:4];
ax.XTickLabel = {'100 Hz', '10 Hz', '10 Hz B', 'TBS'}; %{'10 Hz', '10 Hz B', '100 Hz', 'TBS'};