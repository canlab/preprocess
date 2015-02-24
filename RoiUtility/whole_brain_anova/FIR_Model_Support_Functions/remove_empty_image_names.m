function EXPT = remove_empty_image_names(EXPT,fldname)
% function EXPT = remove_empty_image_names(EXPT,fldname)
%
% removes any empty rows from a cell vector of string matrices containing
% names

str = (['nms = EXPT.' fldname ';']);
eval(str)

for ss = 1:length(nms)
    
    tmpim = nms{ss};

    wh = []; for i=1:size(tmpim,1),if isempty(deblank(tmpim(i,:))), wh(end+1) = i;,end,end

    tmpim(wh,:) = [];
    nms{ss} = tmpim;
    
end

str = (['EXPT.' fldname ' = nms;']);
eval(str)


return

