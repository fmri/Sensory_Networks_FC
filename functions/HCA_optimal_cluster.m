function [cluster_stats] = HCA_optimal_cluster(clusters, distance_measure, names)
%

assert(nargin==3);
[a,b] = size(distance_measure);
assert(a==b && mod(a,2)==0, 'distance_measure must be square and have even number of rows/cols');

cluster_wcs = nan(a, 1);
cluster_ROIs = cell(a,1);
num_clusts = nan(a,1);

clusts = clusters(end);

for cc = 1:a-2
    if cc == 1
        clusts = clusters(end);
        cluster_wcs(cc) = sum(triu(distance_measure), 'all');
        cluster_ROIs{cc} = clusts;
        node_ind = a-2;
        num_clusts(cc) = 1;
        continue
    elseif cc == 2 
        clusts = clusters(node_ind-1:node_ind);
        node_ind = node_ind-1;
    else
        clust_split = clusters{node_ind};
        clusts_prev = cellfun(@(x) x(~ismember(x, clust_split)), clusts, 'UniformOutput', false);
        clusts = [clusts_prev, {clust_split}];
        clusts = clusts(~cellfun('isempty', clusts));
    end
    
    % Find within cluster sum of distance
    wcs = 0;
    for ii = 1:length(clusts)
        clust = clusts{ii};
        ROI_inds = ismember(names, clust);
        clust_mat = distance_measure(ROI_inds, ROI_inds);
        wcs = wcs + sum(clust_mat, 'all');
    end
    cluster_wcs(cc) = wcs;
    cluster_ROIs{cc} = clusts;
    num_clusts(cc) = length(clusts);
    node_ind = node_ind-1;
end

cluster_stats = table(num_clusts, cluster_wcs, cluster_ROIs, 'VariableNames', {'Num_clusters', 'within_cluster_sum', 'cluster_ROIs'});

end
