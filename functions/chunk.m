function [chunk_sequence, n_per_chunk, change_detect] = chunk(num_sequence, count_per_chunk)
%CHUNK

%% Check arguments
assert(nargin>=1, 'At least one argument required for chunk function');
if nargin < 2 || isempty(count_per_chunk)
    count_per_chunk = false;
end

%% Remove consecutive repeats but keep the values in the same order
change_detect = [true diff(num_sequence(:)')~=0];
chunk_sequence = num_sequence(change_detect);

%% Count number of occurences in each chunk
if count_per_chunk
    n_per_chunk = nan(1,length(chunk_sequence));
    count = 1;
    ind_count = 1;
    for ii = 2:length(num_sequence)
        if change_detect(ii)
            n_per_chunk(ind_count) = count;
            ind_count = ind_count + 1;
            count = 1;
        else
            count = count + 1;
        end
    end
    n_per_chunk(end) = count;
else
    n_per_chunk = 'use count_per_chunk=true argument to get n_per_chunk output';
end
end

