function [ROI_names] = extract_borderfile_ROIs(filepath)
%% The purpose of this funciton is to extract the ROI information from a cwb border file
%
% Inputs:
%       filepath: str - filepath to border file
% Outputs:
%       ROI_names: cell array of strs - ROI names contained in border file
%
% Tom Possidente - May 2024

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["Var1", "ClosedTrue", "Var3", "Var4", "Var5", "Var6"];
opts.SelectedVariableNames = ["Var1", "ClosedTrue", "Var3", "Var4", "Var5", "Var6"];
opts.VariableTypes = ["string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "ClosedTrue", "Var3", "Var4", "Var5", "Var6"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "ClosedTrue", "Var3", "Var4", "Var5", "Var6"], "EmptyFieldRule", "auto");

% Import the data
borderfile_table = readtable(filepath, opts);
name_column = borderfile_table.ClosedTrue;

% Convert to strings of just ROI names
all_names =  cellstr(name_column(contains(name_column, 'Name="')));
ROI_names = {};

counter = 0;
for nn = 1:length(all_names)
    name = all_names{nn};
    stripped_name = name(7:end-1);
    if ~any(strcmpi(stripped_name, {'Auditory', 'Tactile', 'Visual', 'Multisensory'}))
        counter = counter + 1;
        ROI_names{counter} = stripped_name;
    end
end

%% Clear temporary variables
clear opts


end

