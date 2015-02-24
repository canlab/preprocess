function fitness = prospect_organism(wh,pop,truep,iter)


x = pop(wh,1);
y = pop(wh,2);
p = pop(wh,3);

% true underlying utility
[trueu] = prospect_utility(truep,x,y,p);
% true report
ureport = vinv(trueu,truep(3),truep(4));

[fitness,covmtx,ste_params,best_params] = prospect_ste_mc(x,y,p,ureport,iter);

return

