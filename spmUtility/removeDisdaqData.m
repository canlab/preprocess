function [trimmedData] = removeDisdaqData(ts, numVolsPerRun, numDisdaqVols, varargin)
% [trimmedData] = removeDisdaqData(ts, numVolsPerRun, numDisdaqVols, ['strict'])
%       ts - timeseries of data to shorten, where each column is a timeseries
%       numVolsPerRun - vector of volume counts per run, not including disdaq vols
%       numDisdaqVols - constant describing how many data points to remove from the beginning of each run
%       'strict' - if set, will error out unless data given to it is
%          exactly proper length - defaults to 1
%
% E.g., removing some physio measures recorded during the disdaq period,
% for an experiment
% numVolsPerRun = [124 140 109];  % NOT including disdaqs
% numDisdaqVols = 4;
% physioData = whatever();  
% trimmedRatingData = removeDisdaqData(physioData, numVolsPerRun, numDisdaqVols);

strict = 1;

if ~isempty(varargin)
    for i = 1:length(varargin)
        if strcmp(varargin{i},'strict'), strict = varargin{i+1};end
    end
end



numRuns = length(numVolsPerRun);
exptLength = sum(numVolsPerRun) + (numRuns * numDisdaqVols);

if(strict && (exptLength ~= size(ts,1)))
    error('The number of vols + disdaqs don''t add up to the length of the timeseries.');
else
    trimmedData = ts(1:exptLength,:);
    st = [1 cumsum(numVolsPerRun(1:end-1)) + 1];
    en = st + numDisdaqVols - 1;
    for i=1:length(st)
        trimmedData(st(i):en(i),:) = [];
    end
end