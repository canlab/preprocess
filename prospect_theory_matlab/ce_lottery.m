%% get stuff set up
%read eprime file in, to get trial choice information. file is exported
%from eprime's .edat file to .xls. this means you don't get the first
%column, but all others import (so should have 67 cols in data, and 68 in
%textdata)

newData1 = importdata(fileToRead1);

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

% data column n label is textdata(2,n+1)... do something with this,
% probably change structure of data; labels don't really matter though.
datalabels = textdata(2,2:68);


%%  select random trial to count

TOTAL = [];
for j = 1:10,
    %get random trial between 1 and 50
    trialnum = ceil(50*rand(1));

    %get extremes of gamble
    value1 = data(trialnum,66);
    value2 = data(trialnum,67);

    % get random number between extremes of gamble
    distance = abs(value1-value2);
    lottery = ceil(distance*rand(1));
    if value1 < value2,
        lottery = value1 + lottery;
    else
        lottery = value2 + lottery;
    end

    % compare lottery value with CE
    %cert_equiv_index is which selection they chose from newtext 1:11 (which map onto 35
    %as nt1, 38-45 for nt2-nt9, and 36-37 for nt10-11)

    index = [35,38:1:45,36,37];
    cert_equiv_index = data(trialnum,21);
    cert_equiv = data(trialnum,index(cert_equiv_index));
    % should confirm that this is midpoint. seems like displayed value?

    if cert_equiv >= lottery, %make sure this works conceptually with negative values.
        !echo Lottery number is less than CE. No gamble!
        lottery = lottery;
        cert_equiv = cert_equiv;
        TOTAL(j) = cert_equiv
    else
        !echo Lottery number is greater than CE. Time to gamble...
        p_value1 = data(trialnum,46); %p of value 1, 1-p is value 2
        lottery2 = rand(1);
        if lottery2 < p_value1,
            !echo Gamble done!
            TOTAL(j) = value1
        else
            !echo Gamble done!
            TOTAL(j) = value2
        end
    end
end