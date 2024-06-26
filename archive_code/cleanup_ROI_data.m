all_ROIs_use = {'tgPCS', 'cIFS/G', 'FO', 'CO', 'pAud', 'pVis', 'preSMA-V', 'SPCS', ...
                'IPCS', 'midIFS', 'cmSFG_mult', 'Ins_mult'};
basedir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';
dirs = {dir(basedir).name};
dirs = dirs(3:end);

for dd = 1:length(dirs)

    filecurr = dirs{dd};
    if ~contains(filecurr, all_ROIs_use) && ~contains(filecurr, '.func')
        unix(['mv ' basedir filecurr ' ' basedir 'unused_ROIs/']);
    end

end