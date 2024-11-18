function swarmchart_PSC_perROI(x, group1, group2, ROIs, xlab, title_text)
%SWARMCHART_PSC 

N = size(group1,1);
N_ROIs = length(ROIs);

figure;
count = 0;
for rr = 1:N_ROIs
    count = count + 1;
    swarmchart(x(count), group1(:,rr), 'filled');
    hold on;
    scatter(x(count), mean(group1(:,rr),'omitnan'), 100, 'r', 'filled');
    count = count + 1;
    swarmchart(x(count), group2(:,rr), 'filled');
    scatter(x(count), mean(group2(:,rr),'omitnan'), 100, 'r', 'filled');
end

xticks(x);
xticklabels(repelem(ROIs,2));
xtickangle(45);
xlabel(xlab);
ylabel('PSC');
title(title_text);
grid on;
set(gcf,'Position',[100 100 1500 800])


end

