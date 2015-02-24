% [xY, Y] = cltoVOI(cls, SPM, ['save'], ['plot'], ['filter', 0|1 ])
%   Returns the xY structure and the Y data contained in VOI*.mat files
%
%   Optional inputs:
%       'save' - save as a VOI .mat file based on the shorttitle of the cluster
%       'plot' - display the cluster and the first computed eigenvector
%       'filter' - HP filter the data and whiten/weight data - default: 1
%
%   NB: This is not exact, as SPM VOIs are run/session-specific, and clusters are not.

function [xYs, Ys] = cltoVOI(cls, SPM, varargin)
    saving_VOI = 0;
    plotting_VOI = 0;
    filtering_data = 1;
    nruns = 1;
    data_field = 'timeseries';

    parse_inputs();

    for i=1:length(cls)
        cl = cls(i);

        xY.xyz = cl.mm_center(:);
        xY.Ic = 0;
        xY.def = 'cluster';
        xY.XYZmm = cl.XYZmm;

        if(~isfield(cl, 'shorttitle') || isempty(cl.shorttitle))
            xY.name = 'VOI';
        else
            xY.name = cl.shorttitle;
        end


        all_y = cl.(data_field);
        if(filtering_data)
            all_y = spm_filter(SPM.xX.K, SPM.xX.W * all_y);
        end

        all_X0 = SPM.xX.xKXs.X(:,[SPM.xX.iB SPM.xX.iG]);
        
        for j=1:length(SPM.Sess)
            xY.Sess = j;
            
            wh_sess = SPM.Sess(xY.Sess).row;
            y = all_y(wh_sess,:);
            xY.X0 = all_X0(wh_sess,:);
            
            if(isfield(SPM.xX.K(xY.Sess), 'X0'))
                xY.X0 = [xY.X0 SPM.xX.K(xY.Sess).X0];
            end
            if(isfield(SPM.xX.K(xY.Sess), 'KH'))
                xY.X0 = [xY.X0 SPM.xX.K(xY.Sess).KH];
            end

            [Y v s] = compute_regional_response(y);

            xY.y = y;
            xY.u = Y;
            xY.v = v;
            xY.s = s;

            xY = order_xY_fields(xY);
            
            if(saving_VOI)
                save_VOI(xY, Y);
            end
            if(plotting_VOI)
                plot_VOI(xY, Y, SPM);
            end

            xYs(i, j) = xY;
            Ys{i,j} = Y;
        end
    end

    function parse_inputs()
        for iv=1:length(varargin)
            if(ischar(varargin{iv}))
                switch(lower(varargin{iv}))
                    case 'save'
                        saving_VOI = 1; %#ok
                    case 'plot'
                        plotting_VOI = 1; %#ok
                    case {'filter' 'filter_data'}
                        filtering_data = varargin{iv + 1}; %#ok
                    case {'nruns' 'num_runs'}
                        nruns = varargin{iv + 1}; %#ok
                    case 'data_field'
                        data_field = varargin{iv + 1}; %#ok
                end
            end
        end
    end
end


function [Y v s] = compute_regional_response(y)
    [m n] = size(y);
    if m > n
        [v s v] = svd(spm_atranspa(y));
        s = diag(s);
        v = v(:,1);
        u = y * v / sqrt(s(1));
    else
        [u s u] = svd(spm_atranspa(y'));
        s = diag(s);
        u = u(:,1);
        v = y' * u / sqrt(s(1));
    end

    d = sign(sum(v));
    u = u * d;
    v = v * d;
    Y = u * sqrt(s(1) / n);
end

function plot_VOI(xY, Y, SPM)
    % show position
    %------------------------------------------------------------------------
    Fgraph = spm_figure('GetWin','Graphics');
    spm_results_ui('Clear', Fgraph);
    figure(Fgraph);
    subplot(1, 2, 1)
    spm_dcm_display(xY, [], [], [[1 0 0];[0 1 0]]', 64)


    % show dynamics
    %------------------------------------------------------------------------
    subplot(1, 2, 2);
    try
        plot(SPM.xY.RT * 1:length(xY.u), Y);
        str = 'time (seconds}';
    catch
        plot(Y);
        str = 'scan';
    end
    title(['1st eigenvariate: ' xY.name], 'FontSize', 10);
    str = {	str;
        ' ';
        sprintf('%d voxels in VOI at [%3.0f %3.0f %3.0f]', length(Y), xY.xyz);
        sprintf('Variance: %0.2f%%', xY.s(1) * 100/sum(xY.s))
        };
    xlabel(str);
    axis tight square
    drawnow
end

function xY = order_xY_fields(xY)
    xY_fieldnames = {'xyz', 'Ic', 'Sess', 'def', 'XYZmm', 'name', 'X0', 'y', 'u', 'v', 's'};
    xY = orderfields(xY, xY_fieldnames);
end

function save_VOI(xY, Y) %#ok
    VOI_mat_name = ['VOI_' xY.name];
    if isfield(xY, 'Sess') && length(xY.Sess) == 1
        VOI_mat_name = sprintf('VOI_%s_%i', xY.name, xY.Sess);
    end
    save(VOI_mat_name, 'Y', 'xY');
end