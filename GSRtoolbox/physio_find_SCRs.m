
function SCRs = physio_find_SCRs(cf, gradperc)
% cf is the handle for the figure in which the data is located
% gradperc is the percentile of slopes which are IDed as markers for
%   potential responses
    
    
    
% SCR detection and quantification
    %
    
    
    plotdata = get(cf, 'UserData');
    
    nsub = length(plotdata.signal);
     nsess = zeros(nsub,1);
    
    for i=1:nsub
        nsess(i) = length(plotdata.signal{i});
    end;
    fs = plotdata.fs;
    trialdata = plotdata.trials;
    
    
    total_loops = 0;
    actionloc = cell(nsub,1);
    heights = cell(nsub,1);
    humpdat = cell(nsub,1);
    for i = 1:nsub
         actionloc{i} = cell(nsess(i),1);
         heights{i} = cell(nsess(i),1);
         humpdat{i} = cell(nsess(i),1);
         total_loops = total_loops + nsess(i);
    end;
    %
    
    
    if plotdata.parameters.detect == 4
        threshhold = pick_percentile(gradient(concatcell(trialdata),1/fs), gradperc);
    end;
    
    spm_progress_bar('init');
    total_item = 0;
    for i = 1:nsub
        
        if plotdata.parameters.detect == 3
               threshhold = pick_percentile(gradient(concatcell(trialdata{i}),1/fs{i}), gradperc);
        end;
        for j = 1:nsess(i)  
            
            if plotdata.parameters.detect == 1
         
                for k = 1:size(trialdata{i}{j},1)
                    
                    trialdat = plotdata.clean_signal{i}{j}(trialdata{i}{j}(k,1):trialdata{i}{j}(k,2));
                    threshhold = pick_percentile(gradient(trialdat,1/fs{i}{j}), gradperc);
                    
                    % Do the SCR detection
                    % smoothingw may need to be a separate field for
                    % derivative smoothing...
                    [temptime, tempheight, temphump] = SCRdetect5(trialdat, threshhold, fs{i}{j}, ...
                        plotdata.parameters.smoothingw, plotdata.parameters.minperiod/1000, 0, 'noplot');
                    
                    if ~isempty(temptime)
                        temptime = tempstart + trialdata{i}{j}(k,1);
                        if ~isempty(temphump)
                            temphump(:,1) = temphump(:,1) + trialdata{i}{j}(k,1);
                            humpdat{i}{j} = [humpdat{i}{j}; temphump];
                        end;
                        actionloc{i}{j} = [actionloc{i}{j}; temptime];
                        heights{i}{j} = [heights{i}{j}; tempheight];
                    
                    end;
                end;
            else
                
                if plotdata.parameters.detect == 2
                    threshhold = pick_percentile(gradient(plotdata.clean_signal{i}{j},1/fs{i}{j}), gradperc);
                end

                % Do the SCR detection
                % smoothingw may need to be a separate field for
                % derivative smoothing...
                [actionloc{i}{j}, heights{i}{j}, humpdat{i}{j}] = SCRdetect5(plotdata.clean_signal{i}{j}, threshhold, fs{i}{j}, ...
                    plotdata.parameters.smoothingw, plotdata.parameters.minperiod/1000, 0, 'noplot');

            end
            total_item = total_item + 1;
            spm_progress_bar('set', total_item/total_loops);
        end
        
    end
    
    spm_progress_bar('clear');
    
    SCRs.location = actionloc;
    SCRs.height = heights;
    SCRs.humps = humpdat;
end



function output = concatcell(incell)
    if isnumeric(incell)
        output = incell;
    elseif iscell(incell)
        if isnumeric(incell{1}) %is this the lowest level of the cell?
            
            output = [];
            for i = 1:length(incell)
                output = [output; incell{i}];
            end;
            
        else %not the lowest level, call self recursively
            output = [];
            for i = 1:length(incell)
                output = [output; concatcell(incell{i})];
            end;
        end;
    end
end


function value = pick_percentile(inarray, perc)
    temp = sort(inarray);
    pick_ind = perc * length(inarray);
    
    if pick_ind ~= floor(pick_ind)
        floorco = pick_ind-floor(pick_ind);
        value = temp(floor(pick_ind))*floorco + temp(ceil(pick_ind))*(1-floorco);
    else
        value = temp(pick_ind);
    end;
end
    