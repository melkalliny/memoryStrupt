function [] = plotGroupConsistencyCorrelations(correlations_testGroup, correlations_controlGroup, locationStrings, regionToDo)




for region=regionToDo
    temp_data = squeeze(correlations_testGroup(region,:,4));
    temp_behavior = squeeze(correlations_testGroup(region,:,5));
    valid = find(~cellfun(@isempty,temp_data));
    gcf = figure('Units', 'Normalized', 'OuterPosition', [0 0 0.8 0.8]);
    clf
    for i=1:length(valid)
        subplot(1,length(valid),i)
        scatter(temp_data{1,valid(i)},temp_behavior{1,valid(i)},'o');
        lsline
        [a,pval] = corr(temp_data{1,valid(i)},temp_behavior{1,valid(i)},'rows','complete','type','Spearman');
        xlim([-1 1]); ylim([-2000 2000]);
        xlabel('mean pairwise consistency')
        ylabel('mean reaction time change')
        title(sprintf('window size 500ms, corr %0.2f, pval %0.2f',a, pval));
        suptitle(sprintf('location %s',locationStrings{1,region}))
        set(gca, 'FontName', 'Helvetica');
        set(gca,'FontSize',18)
    end
    
end


test_data = squeeze(correlations_testGroup(region,:,1));
valid = find(~cellfun(@isempty,test_data));
test_data = cell2mat(test_data(valid));
control_data = squeeze(correlations_controlGroup(region,:,1));
valid = find(~cellfun(@isempty,control_data));
control_data = cell2mat(control_data(valid));


gcf = figure('Units', 'Normalized', 'OuterPosition', [0 0 0.8 0.8]);
clf
jitter = 0.2; markerSize = 8;
for ii=1:length(test_data)
    tmp=test_data(ii);
    x=1;
    x = x+(rand(size(x))-0.5)*jitter; %add a little random "jitter" to aid visibility
    plot(x,tmp,'.r','markersize',markerSize)
end
% lets overlay mean and ste
hold on
% databar = bar(session,[nanmean(data)],0.5,'FaceColor','none','EdgeColor','k','LineWidth',2);
hold on
ebar = errorbar(1,[nanmean(test_data)],[nanstd(test_data) / sqrt(length(find(~isnan(test_data))))],'Color','r','LineWidth',8,'LineStyle','none');



hold on
jitter = 0.2; markerSize = 8;
for ii=1:length(control_data)
    tmp=control_data(ii);
    x=2;
    x = x+(rand(size(x))-0.5)*jitter; %add a little random "jitter" to aid visibility
    plot(x,tmp,'.r','markersize',markerSize)
end
% lets overlay mean and ste
hold on
% databar = bar(session,[nanmean(data)],0.5,'FaceColor','none','EdgeColor','k','LineWidth',2);
hold on
ebar = errorbar(2,[nanmean(control_data)],[nanstd(control_data) / sqrt(length(find(~isnan(control_data))))],'Color','r','LineWidth',8,'LineStyle','none');

suptitle(sprintf('location %s',locationStrings{1,region}))
set(gca, 'FontName', 'Helvetica');
set(gca,'FontSize',18)

end

