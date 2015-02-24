function [match2,manip,vec,vec2] = reorder(actionstring,loc,match2,manip)
% [match2,manip,vec,vec2] = reorder(actionstring,loc,manip)
% loc
%   string, string matrix, or cell array of strings to search for
%   '051505jk' or {'4444xx' '33333xx'}
%
% ACTION STRINGS:
%
% 'delete'  deletes entries in manip where match2 entries match loc entries
% --------------------------------------------------------------------------
% [match2,manip] = reorder('delete','0619rbc2',dd,dd);
% [match2,manip] = reorder('delete',{'0619rbc2' '0530aa'},dd,dd);
%
% 'findempty'  finds empty cells in manip, deletes them, and returns match2
%              values for empty manip cells
% --------------------------------------------------------------------------
% [avgwave_notempty,d_empty,vec] = reorder('findempty','',avgwave,d);
%
% 'reorder'  re-orders match2 and manip to match the order of loc
%            in the example, matches subj and beh to sorder
%            vec returns indices of subj/beh needed to get newsubj
%            so, subj(vec) = newsubj
%            vec2 returns indices of sorder with matching subj
%            so, sorder(vec2) = newsubj
% --------------------------------------------------------------------------
%[newsubj,newbeh,vec] = reorder('reorder',sorder,subj,beh);
% [newsubj,newbeh,vec] = reorder('reorder',cell_strings_in_desired_order,cell_strings_to_order,beh_corresponding to cell strings to order);
% just find orders
%[dummy,dummy,vec] = reorder('reorder',submtx(:,1:8),d,[]);
%
% example:
% xlorder = parse_char_to_cell('reapp1	reapp2	reapp3	reapp4	reapp5	rea1	rea2	rea3	rea4	rea5	rea6	rea7	rea8	rea9	rea10	rea11	rea13	rea14	rea15	rea16	rea17	rea18	rea19	rea20	rea21	rea22	rea23	rea24	rea25	rea26	rea27	rea28	rea29	rea30	rea31	rea32	rea33', 'tab')'
% xlgender = parse_char_to_cell('M	M	F	M	F	F	M	F	F	M	M	M	F	F	F	M	M	F	M	M	M	F	M	F	F	F	F	F	M	F	F	M	F	M	F	M	M', 'tab')'
% [newsubj,newbeh,vec] = reorder('reorder',EXPT.subjects', xlorder, xlgender);



if isempty(manip), manip = zeros(size(match2));,end

% MAKE VECTOR OF WHICH CELLS/INDICES MATCH

docell = 1;
if ~iscell(loc), tmp = loc; loc = {}; 
    for i = 1:size(tmp,1)
        loc{i} = tmp(i,:);
    end
end

if ~iscell(match2), docell = 0; tmp = [];
    for i = 1:size(match2,1), tmp{i} = deblank(match2(i,:)); end; match2 = tmp;
end

vec = []; vec2 = [];
for L = 1:length(loc)
    for i = 1:length(match2)
        wh = find(strcmp(loc{L},match2{i}));  % findstr(loc{L},match2{i});                    % matching L
        if ~isempty(wh), vec = [vec i];, vec2 = [vec2 L]; end
    end
end

% PERFORM ACTION

switch actionstring
    
case 'reorder'
    match2 = match2(vec);
    if iscell(manip), 
        manip = manip(vec);
    else
        manip = manip(vec,:);
    end
    
    
case 'delete'
    
    match2(vec) = [];
    
    if iscell(manip), 
        manip(vec) = [];,
    else
        manip(vec,:) = [];
    end
    
    
    
case 'findempty'
    
    vec = [];
    for i = 1:length(match2)
        if isempty(match2{i}), vec(end+1) = i;,end
    end
    
    match2(vec) = [];
    
    if iscell(manip), 
        manip = manip(vec);
    else
        manip = manip(vec,:);
    end
    
    
    
otherwise
    
    disp('Unknown action')
    return
    
end


if ~docell
    match2 = str2mat(match2);
end
    
    
return


