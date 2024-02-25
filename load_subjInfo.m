function [subjDf] = load_subjInfo(filepath)
%LOAD_SUBJINFO loads csv with various subject information into a table

% If not filepath provided, used default location
if nargin == 0
 filepath = '/projectnb/somerslab/scripts/jupyter/subjectInfo.csv';
end

% force the date and t1runs columns to be read in as char instead of double (which would produce NaNs for those with multiple dates separated by '/')
opts = detectImportOptions('/projectnb/somerslab/scripts/jupyter/subjectInfo.csv');
date_columns = {'x1WayLocalizerDate', 'spacetimeDate', 'restDate', 'longDelayDate', 'retinoDate', 'speechDate', 'languageDate', ...
    't1Date', 't2Date', 't1Runs'};
column_names = opts.VariableNames;
date_column_inds = ismember(column_names, date_columns);
vartypes = opts.VariableTypes;
vartypes(date_column_inds) = {'char'};
opts.VariableTypes = vartypes;

subjDf = readtable(filepath, opts);

end

