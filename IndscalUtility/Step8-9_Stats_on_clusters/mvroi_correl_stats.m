function c = mvroi_correl_stats(DATA,names)
% function c = mvroi_correl_stats(DATA,names)
%
% prints correlations for average and difference among conditions
%
% tor wager
%where = DATA.where;

comp = DATA.SPEC.comps(1,:);            % first contrast specified

% -------------------------------
% get data in right format
% -------------------------------

% Get AVERAGE within subject across conditions

[avgcor,states] = contrast3d(DATA.DATA.xc,ones(size(comp))./length(comp));  % Conditions are nested within subjects
%avgcor = avgcor(where,where,:);

% get CONTRAST across conditions, for each subject, specified in DATA.SPEC
[difcor,states] = contrast3d(DATA.DATA.xc,comp);              % Conditions are nested within subjects
%difcor = difcor(where,where,:);

c.avgcor = avgcor;
c.difcor = difcor;
c.states = states;


% get average for EACH SEPARATE state
for i = 1:length(comp)
    thiscond = zeros(size(comp)); thiscond(i) = 1;
    
    [statecor] = contrast3d(DATA.DATA.xc,thiscond);              % Conditions are nested within subjects
    c.statecorrs{i} = statecor;
end

% -------------------------------
% stats on data
% -------------------------------

fprintf(1,'\nAverage across conditions\n------------------------------\n)');
[c.avgcor_mean, dummy, c.AVGSTATS] = xc_stats(avgcor,[],names);

fprintf(1,['\nContrast across conditions [' repmat('%3.0f ',1,length(comp)) ']\n------------------------------\n)'],comp);
[c.difcor_mean, dummy, c.DIFSTATS] = xc_stats(difcor,[],names);


% get stats for EACH SEPARATE state
for i = 1:length(comp)
    fprintf(1,'\nSeparate analysis for state %3.0f\n------------------------------\n)',i);
    [c.statecor_mean{i}, dummy, c.STATESTATS(i)] = xc_stats(c.statecorrs{i},[],names);
end


return

