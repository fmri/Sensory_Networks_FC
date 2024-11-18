function swarmchart_PSC(x, group1, group2, masks, xlab, title_text)
%SWARMCHART_PSC 

N = size(group1,1);
mask1 = masks(:,1);
mask2 = masks(:,2);
mask3 = masks(:,3);

figure; 
swarmchart(repmat(x,N,1), [mean(group1(:,mask1),2,'omitnan'), mean(group1(:,mask2),2,'omitnan'), mean(group1(:,mask3),2,'omitnan'),...
                           mean(group2(:,mask1),2,'omitnan'), mean(group2(:,mask2),2,'omitnan'), mean(group2(:,mask3),2,'omitnan')], 'filled');
hold on;
scatter(x, [mean(group1(:,mask1),'all','omitnan'), mean(group1(:,mask2),'all','omitnan'), mean(group1(:,mask3),'all','omitnan'),...
                           mean(group2(:,mask1),'all','omitnan'), mean(group2(:,mask2),'all','omitnan'), mean(group2(:,mask3),'all','omitnan')], 100, 'r', 'filled');
xticks(x);
xticklabels({'Vis ROIs', 'Aud ROIs', 'MD ROIs', 'Vis ROIs', 'Aud ROIs', 'MD ROIs'});
xtickangle(45);
xlabel(xlab);
ylabel('PSC');
title(title_text);
grid on;

end

