function [clusters, clusters_txt] = linkage_output_extract(Z, names)
%LINKAGE_OUTPUT_EXTRACT Extract the clusters from the "linkage" function output

assert(nargin == 2);

N_clusters = size(Z,1);
N_leafs = N_clusters+1;
assert(N_leafs==length(names));

clusters = {};
clusters_txt = {};

for cc = 1:N_clusters
    cluster_IDs = Z(cc,[1,2]);
    flag = false;
    %cluster_txt = [];
    for ii = 1:2
        if any(cluster_IDs(ii)>N_leafs)
            cluster{ii} = clusters{cluster_IDs(ii)-N_leafs};
            %cluster_txt = [cluster_txt clusters_txt{cluster_IDs(ii)-N_leafs}];
            flag=true;
        else
            cluster{ii} = names{cluster_IDs(ii)};
            %cluster_txt = [cluster_txt cluster{ii}];
        end
    end
    if flag
        clusters{cc} = [cluster{1}, cluster{2}];
    else
        clusters{cc} = {cluster{1}, cluster{2}};
    end
    clusters_sorted = sort(clusters{cc});
    clusters_txt{cc} = [clusters_sorted{:}];

end

