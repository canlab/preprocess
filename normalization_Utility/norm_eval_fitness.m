function [fitness,tempvec,normvec,yhat_sorted] = norm_eval_fitness(tempvec,normvec,doplot)
% [fitness,tempvec,normvec,yhat_sorted] = norm_eval_fitness(tempvec,normvec,doplot)
%
% tor wager, nov. 06
%
% fitness is % variance explained by monotonic fit

dosimpleprob = 0;


if dosimpleprob
    yhat_sorted = [];

    % treat these vectors as gray matter probability maps
    % normalize to same total gray matter volume/distribution for each of template and
    % image
    % scale to avoid very small numbers

    % max fitness is 100,000  , min is 0, if every template vox is 1 or 0,
    % and every image vox is 1 or 0 and matches template vox...
    tempvec = 100000 * tempvec ./ sum(tempvec);
    normvec = 100000 * normvec ./ sum(normvec);

    fitness = tempvec' * normvec;

else
    % whole-image version

    % get rid of extraneous voxels
    % this should now be done ahead of time by passing in vectors for tempdat
    % and normdat of only voxels of interest!
    % % % tempvec = tempdat(:);
    % % % normvec = normdat(:);
    % % % whomit = tempvec < .1 | normvec < 50 | isnan(tempvec) | isnan(normvec);
    % % % tempvec(whomit) = [];
    % % % normvec(whomit) = [];

    % round to make ranking more manageable
    normvec = round(100 .* normvec ./ mean(abs(normvec)));
    tempvec = round(100 .* tempvec ./ mean(abs(tempvec)));


    % monotonic regression
    [yhat,yhat_sorted,err,fitness] = monotonic_regression(tempvec,normvec,doplot);
    %fitness = npos ./ df;

    %fitness = corrcoef(tempvec,normvec);
    %fitness = fitness(1,2);

    %if isnan(fitness), keyboard, end
end


return

