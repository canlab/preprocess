function [XX,DXX,classXX,Xlin,rtout] = c2rtpredictor(c,rt,cnames,Sess,L,TR,evtspersession,wh,varargin)
%[XX,DXX,classXX,out,Xlin,rtout] = c2rtpredictor(c,rt,cnames,Sess,L,TR,evtspersession,wh,[even/odd])
%
% This function takes a c input vector of onset times, and rt cell array, 
% from tor's spm modeling scripts, and makes design matrices
%
% c         cell array of input times, in TRs
% rt        cell array of reaction times on each trial
% cnames    cell array, names of event types
% Sess      number of sessions
% L         number of images in experiment
% TR        repetition time
% evtsper   in c, number of event types per session, for [c rt], 
% wh        which event types to extract (vector, e.g., [3:6])
% [even/odd]    optional, 1 or 2 codes get odd or even sessions
%
% calls:
% scale.m
% pad.m
%
% RETURNS three matrices:
% XX, containing predictors and linear/quadratic pred x RT
% DXX, containing deconvolution matrix
% classXX, containing categorized predictors for fast, med, and slow trials,
% respectively
%
% X has basic, then linear, then quad columns in order
% classXX has fast, med, slow RT predictors for each type, within event
% types
% DXX has timepoints within event types, then event types in order
%
% YOU should remove session means separately - it's NOT done for you here!
%
% by Tor Wager

if ~any(evtspersession - mean(evtspersession))
    % check!  this only works of all evtspersession are the same.
    for i = 1:length(c), sizec(i) = length(c{i});,end
    for i = 1:length(rt), sizert(i) = length(rt{i});,end
    if length(sizec) ~= length(sizert), warning('Num of conditions in c (onsets) and rt (reaction time) vectors do not match. Wrong # conditions in one.'),keyboard,end
    wh2 = find(sizec - sizert);
    if any(wh2), warning(['Discrepancies found in conditions: ' num2str(wh2)]);, keyboard, else, end
end

ind = 1; DXX = []; classXX = []; rtout = [];

for E = wh

    %tmp = cnames(E:evtspersession(1):end); 
    t = c(E:evtspersession(1):end);
    r = rt(E:evtspersession(2):end);

    for i = 2:length(t),                          % list of onsets in TRs
        t{i} = t{i} + (i-1)*(L./Sess); ,
    end    

    if length(varargin) > 0
        evodd = varargin{1};
        t = t(evodd:2:end);
        r = r(evodd:2:end);
    end
    
    % eliminate empties from both lists
    wh2 = [];
    for i = 1:length(t), if isempty(t{i}) | isempty(r{i}), wh2(end+1) = i; end, end
    t(wh2) = []; r(wh2) = [];
    
    % t2 is onsets in s, and r2 is rts 
    clear t2 r2
    t2{1} = cat(1,t{:}) .* TR; 
    r2{1} = cat(1,r{:});
    [X,d,out] = rt2delta(t2,r2,TR);                    % "TR" of 1 to preserve, as TRs are input; 2nd input is TR for conv
    
    %if E ~= wh(end), X(:,end) = [];, end        % not the last one
    XX(:,ind) = pad(X(:,1),L-size(X,1));
    XXl(:,ind) = pad(out.rtlinearX,L-size(X,1));
    XXq(:,ind) = pad(out.rtquadX,L-size(X,1));
    
    rtout(:,ind) = pad(out.rtlineard,L-size(X,1));
    
    classXX = [classXX pad(out.rtclassX,L-size(out.rtclassX,1))];
    
                                        % make deconv matrix
    [DX] = tor_make_deconv_mtx3(out.basicd{1},20./TR,1);
    DX = pad(DX(:,1:end-1),L-size(X,1));

    DXX = [DXX DX];
    ind = ind + 1;
end

Xlin = XXl;
XX = [XX XXl XXq];
fprintf(1,'Cond nos: Xlin = %3.2f, Xlin+quad = %3.2f, DX = %3.2f, classX = %3.2f\n',cond(Xlin),cond(XX),cond(DXX),cond(classXX))

XX(:,end+1) = 1;
DXX(:,end+1) = 1;
classXX(:,end+1) = 1;


return

