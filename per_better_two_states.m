function [ em_better, unem_better, total_better ] = per_better_two_states( current,last,store, sim_e )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
sim_e = sim_e - 1;
employed_n = sum(sim_e);
unemployed_n = size(sim_e,2) - employed_n;
aux = store(current).U - store(last).U;
aux_em = aux(sim_e==1);
aux_unem = aux(sim_e==0);
em_better = (sum(aux_em>0)/employed_n)*100;
unem_better = (sum(aux_unem>0)/unemployed_n)*100;
total_better = (sum(aux>0)/size(sim_e,2))*100;

end

