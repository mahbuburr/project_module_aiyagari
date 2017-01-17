function [ better ] = per_better( current,last,store )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
sim_e = sim_e - 1;
employed_n = sum(sim_e);
unemployed_n = size(sim_e,2) - employed_n;
aux = store(current).U - store(last).U;
aux = aux(sim_e==1);
aux_unem = aux(sim_e==0);
em_better = (sum(aux_em(aux_em>0))/employed_n)*100;
unem_better = (sum(aux_unem(aux_unem>0))/unemployed_n)*100;

end

