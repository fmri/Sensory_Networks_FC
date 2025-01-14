function exclusive_regions = find_exclusive_regions(regions)
% FIND_EXCLUSIVE_REGIONS
% Finds items in binary columns exclusive to each column

assert(size(regions,2)>1, 'input must have more than 1 column')

ncols = size(regions, 2);

exclusive_regions = nan(size(regions));
for cc = 1:ncols
    exclusive_regions(:,cc) = (regions(:,cc) - sum(regions(:,1:end~=cc),2) ) == 1;
end

end

