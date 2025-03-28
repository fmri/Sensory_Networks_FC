function [prob_clusters_tbl] = HCA_bootstrap(conn_matrices, iters, clusters_orig, clusters_orig_txt, names)
%HCA_BOOTSTRAP run bootstrap analysis on hierarchical clustering analysis

assert(nargin>=4, 'at least 4 inputs required for function HCA_bootstrap');

N_leafs = length(names);

%% Extract clusters from original data
% mean_connmat_z_orig = mean(conn_matrices_z, 3, 'omitnan');
% mean_connmat_orig = (exp(2.*mean_connmat_z_orig) - 1) ./ (exp(2.*mean_connmat_z_orig) + 1); % inverse fisher transform
% distmat_orig = 1-abs(mean_connmat_orig);
% distmat_orig(1:1+N_leafs:end) = 0; % Make diagonal distance zero
% linkage_cluster_orig = linkage(distmat_orig, 'average');
%[clusters_orig, clusters_orig_txt] = linkage_output_extract(linkage_cluster_orig, names);

%% Get bootstrap sampling indices
N = size(conn_matrices,3);
sample_inds = randi(N, [N,iters]);

%% Loop through interations and calculate HCA
bs_clusters = cell(iters,1);
prob_clusters = NaN(iters,length(clusters_orig_txt));

for ii = 1:iters

    connmats_bs = conn_matrices(:,:,sample_inds(:,ii));
    mean_connmat_bs = mean(connmats_bs, 3, 'omitnan');
    distmat_bs = 1-abs(mean_connmat_bs);
    distmat_bs(1:1+N_leafs:end) = 0; % Make diagonal distance zero
    linkage_cluster_bs = linkage(distmat_bs, 'ward');
    [~, bs_clusters{ii}] = linkage_output_extract(linkage_cluster_bs, names);

    for cc = 1:length(clusters_orig_txt) % loop through each original cluster and see if it's in this bootstrapped data cluster
        cluster = clusters_orig_txt{cc};
        prob_clusters(ii,cc) = ismember(cluster, bs_clusters{ii});
    end

end

prob_clusters_tbl = table(clusters_orig',[mean(prob_clusters,1)]');
prob_clusters_tbl.Properties.VariableNames = {'cluster', 'bs_probability'};

end

