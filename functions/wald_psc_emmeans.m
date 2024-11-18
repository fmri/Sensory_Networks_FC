function [wald_res_tbl] = wald_psc_emmeans(lme,emm)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

wald_res_tbl = table();
N_conds = height(emm.table);
for ii = 1:N_conds
    for jj = 1:N_conds
        if jj>ii
            contrast = zeros(1,N_conds);
            contrast(ii) = 1;
            contrast(jj) = -1;
            res_table = contrasts_wald(lme, emm, contrast);
            wald_res_tbl = [wald_res_tbl; {emm.table.Row{ii}, emm.table.Row{jj}, res_table.Wald, res_table.pVal}];
        end
    end
end

wald_res_tbl.Properties.VariableNames = {'Condition1', 'Condition2', 'Wald_test_stat', 'pval'};