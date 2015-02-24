% Bold deconvolution to create physio- or psycho-physiologic interactions
% FORMAT PPI = peb_ppi(SPM, ppi_type, ppi_name, VOI_files);
%
% See spm_peb_ppi.m for more

function [PPI U] = peb_ppi(SPM, ppi_type, ppi_name, VOI_files, varargin)
    U = [];
    show_results = 1;

    for i=1:length(varargin)
        if(ischar(varargin{i}))
            switch(varargin{i})
                case {'psych inputs', 'psychological inputs', 'U'}
                    U = varargin{i+1};
                case {'show', 'display'}
                    show_results = varargin{i+1};
            end
        end
    end
    
    % set up the graphical interface
    %----------------------------------------------------------------------
    if(show_results || isempty(U))
        Finter = spm_figure('GetWin', 'Interactive');
        Fgraph = spm_figure;
        header = get(Finter, 'Name');
    end

    % check inputs and set up variables
    %----------------------------------------------------------------------
    if ~nargin
        error('Must supply an SPM variable');
    end
    RT     = SPM.xY.RT;
    dt     = SPM.xBF.dt;
    NT     = RT/dt;


    % Ask whether to perform physiophysiologic or psychophysiologic interactions
    %--------------------------------------------------------------------------
    if(show_results || isempty(U))
        set(Finter, 'name', 'PPI Setup');
    end

    switch ppi_type
        case  'simple deconvolution'
            P = cellstr(VOI_files);
            if(size(P, 1) ~= 1), error('Must input 1 VOI'); end
            p      = load(P{1}, 'xY');
            xY(1)  = p.xY;
            Sess   = SPM.Sess(xY(1).Sess);

        case  'physiophysiologic interaction' % interactions between 2 regions
            P = cellstr(VOI_files);
            if(size(P, 1) ~= 2), error('Must input 2 VOIs'); end
            for i=1:2
                p      = load(P{i}, 'xY');
                xY(i)  = p.xY;
            end
            Sess   = SPM.Sess(xY(1).Sess);

        case  'psychophysiologic interaction'  % get hemodynamic response
            P = cellstr(VOI_files);
            if(size(P, 1) ~= 1), error('Must input 1 VOI'); end
            p      = load(P{1}, 'xY');
            xY(1)  = p.xY;
            Sess   = SPM.Sess(xY(1).Sess);

            % get 'causes' or inputs U
            %----------------------------------------------------------------------
            U = setup_inputs(Sess, U);

        otherwise
            error('Unknown ppi_type: "%s". Must be one of "simple deconvolution", "physiophysiologic interaction", or "psychophysiologic interaction"', ppi_type);
    end


    % name of PPI file to be saved
    %-------------------------------------------------------------------------
    %     PPI.name    = spm_input('Name of PPI', 3, 's', 'PPI');
    PPI.name = ppi_name;


    % Setup variables
    %-------------------------------------------------------------------------
    N     = length(xY(1).u);
    k     = 1:NT:N*NT;  			% microtime to scan time indices


    % create basis functions and hrf in scan time and microtime
    %-------------------------------------------------------------------------
    %     spm('Pointer', 'watch')
    hrf   = spm_hrf(dt);


    % create convolved explanatory {Hxb} variables in scan time
    %-------------------------------------------------------------------------
    xb    = spm_dctmtx(N*NT + 128, N);
    Hxb   = zeros(N, N);
    for i = 1:N
        Hx       = conv(xb(:,i), hrf);
        Hxb(:,i) = Hx(k + 128);
    end
    xb    = xb(129:end,:);


    % get confounds (in scan time) and constant term
    %-------------------------------------------------------------------------
    X0    = xY(1).X0;
    M     = size(X0, 2);


    % get response variable,
    %-------------------------------------------------------------------------
    for i = 1:size(xY, 2)
        Y(:,i) = xY(i).u;
    end


    % remove confounds and save Y in ouput structure
    %-------------------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WOOOT! WOOOT! WOOOT!
    % The below lines are the *key* change for multi-session PPIs. A bug in SPM2 
    % and SPM5 has X0 include the zeros of intercept terms for other sessions.. a 
    % bug. However, when you take the inv, you get NaNs and gibberish. The key is 
    % to either remove the zero columns first (probably ideal) or use pinv, which 
    % I've chosen for speed's sake. Either produce nearly identical results it 
    % seems. Lastly, you can, of course, change your design to that of a single 
    % session, but you should not have to redesign around bad code. For more, see
    % these threads:
    %
    % http://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind05&L=SPM&D=0&I=-3&m=13917&P=232778
    % http://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind04&L=SPM&P=R351544&I=-3&m=13917
    % 
    %Yc    = Y - X0*inv(X0'*X0)*X0'*Y;
    Yc    = Y - X0*pinv(X0'*X0)*X0'*Y;
    PPI.Y = Yc(:,1);
    if size(Y, 2) == 2
        PPI.P  = Yc(:,2);
    end


    % specify covariance components; assume neuronal response is white
    % treating confounds as fixed effects
    %-------------------------------------------------------------------------
    Q      = speye(N, N)*N/trace(Hxb'*Hxb);
    Q      = blkdiag(Q, speye(M, M)*1e6  );

    % get whitening matrix (NB: confounds have already been whitened)
    %-------------------------------------------------------------------------
    W      = SPM.xX.W(Sess.row, Sess.row);

    % create structure for spm_PEB
    %-------------------------------------------------------------------------
    P{1}.X = [W*Hxb X0];		% Design matrix for lowest level
    P{1}.C = speye(N, N)/4;		% i.i.d assumptions
    P{2}.X = sparse(N + M, 1);	% Design matrix for parameters (0's)
    P{2}.C = Q;


    switch ppi_type
        case  'simple deconvolution'
            C       = spm_PEB(Y, P);
            xn      = xb*C{2}.E(1:N);
            xn      = spm_detrend(xn);

            % save variables
            %---------------------------------------------------------------------
            PPI.xn  = xn;

            % Plot so the user can see the results
            %---------------------------------------------------------------------
            if(show_results)
                figure(Fgraph);
                t       = RT*[1:N];
                T       = dt*[1:(N*NT)];

                subplot(2, 1, 1)
                plot(t, Yc, T, PPI.xn)
                title('hemodynamic and neuronal responses')
                xlabel('time (secs)')
                axis tight square
                grid on
                legend('BOLD', 'neuronal')
            end

        case  'physiophysiologic interaction' % PHYSIOPHYSIOLOGIC INTERACTIONS
            C       = spm_PEB(Y(:,1), P);
            xn1     = xb*C{2}.E(1:N);
            C       = spm_PEB(Y(:,2), P);
            xn2     = xb*C{2}.E(1:N);
            xn1     = spm_detrend(xn1);
            xn2     = spm_detrend(xn2);
            xnxn    = xn1.*xn2;

            % convolve and resample at each scan for bold signal
            %---------------------------------------------------------------------
            ppi     = conv(xnxn, hrf);
            ppi     = ppi(k);

            % save variables
            %---------------------------------------------------------------------
            PPI.xn  = [xn1 xn2];
            PPI.ppi = spm_detrend(ppi);


            % Plot so the user can see the results
            %---------------------------------------------------------------------
            if(show_results)
                figure(Fgraph);
                t       = RT*[1:N];
                T       = dt*[1:(N*NT)];

                subplot(2, 1, 1)
                plot(t, PPI.ppi)
                title('PPI')
                xlabel('time (secs)')
                axis tight square
                grid on

                subplot(2, 2, 3)
                plot(t, Yc(:,1), T, PPI.xn(:,1))
                title('hemodynamic and neuronal responses (1st)')
                xlabel('time (secs)')
                axis tight square
                grid on
                legend('BOLD', 'neuronal')


                subplot(2, 2, 4)
                plot(t, Yc(:,2), T, PPI.xn(:,2))
                title('hemodynamic and neuronal responses (2nd)')
                xlabel('time (secs)')
                axis tight square
                grid on
                legend('BOLD', 'neuronal')
            end


        case  'psychophysiologic interaction'
            % COMPUTE PSYCHOPHYSIOLOGIC INTERACTIONS
            % use basis set in microtime
            %---------------------------------------------------------------------
            % get parameter estimates and neural signal; beta (C) is in scan time
            % This clever trick allows us to compute the betas in scan time which is
            % much quicker than with the large microtime vectors. Then the betas
            % are applied to a microtime basis set generating the correct neural
            % activity to convolve with the psychological variable in mircrotime
            %---------------------------------------------------------------------
            C       = spm_PEB(Y, P);
            xn      = xb*C{2}.E(1:N);
            xn      = spm_detrend(xn);

            % setup psychological variable from inputs and contast weights
            %---------------------------------------------------------------------
            PSY     = zeros(N*NT, 1);
            for i = 1:size(U.u, 2)
                PSY = PSY + full(U.u(:,i)*U.w(:,i));
            end
            PSY     = spm_detrend(PSY);

            % multiply psychological variable by neural signal
            %---------------------------------------------------------------------
            PSYxn   = PSY.*xn;

            % convolve and resample at each scan for bold signal
            %---------------------------------------------------------------------
            ppi	    = conv(PSYxn, hrf);
            ppi     = ppi(k);

            % similarly for psychological effect
            %---------------------------------------------------------------------
            PSYHRF  = conv(PSY, hrf);
            PSYHRF  = PSYHRF(k);

            % save psychological variables
            %---------------------------------------------------------------------
            PPI.psy = U;
            PPI.P   = PSYHRF;
            PPI.xn  = xn;
            PPI.ppi = spm_detrend(ppi);


            % Plot so the user can see the results
            %---------------------------------------------------------------------
            if(show_results)
                figure(Fgraph);
                t       = RT*[1:N];
                T       = dt*[1:(N*NT)];

                subplot(2, 1, 1)
                plot(t, Yc(:,1), T, PPI.xn(:,1))
                title('hemodynamic and neuronal responses')
                xlabel('time (secs)')
                axis tight square
                grid on
                legend('BOLD', 'neuronal')

                subplot(2, 2, 3)
                plot(t, PPI.P, T, PSY, '--')
                title('[convolved] psych. variable')
                xlabel('time (secs)')
                axis tight square
                grid on

                subplot(2, 2, 4)
                plot(t, PPI.ppi)
                title('PPI')
                xlabel('time (secs)')
                axis tight square
                grid on
            end
    end % (switch)

    % setup other output variables
    %-------------------------------------------------------------------------
    PPI.xY  = xY;
    PPI.dt  = dt;
    str     = ['PPI_' PPI.name];
    save(fullfile(SPM.swd, str), 'PPI');

    % clean up
    %-------------------------------------------------------------------------
    %     spm('Pointer', 'arrow')
    %     spm('FigName', header);
    fprintf('\ndone\n')
end


function U = setup_inputs(Sess, U)
    if(isempty(U))
        spm_input('Psychological variable:...  ', 2, 'd');
        u      = length(Sess.U);
        U.name = {};
        U.u    = [];
        U.w    = [];
        for  i = 1:u
            for  j = 1:length(Sess.U(i).name)
                str   = ['include ' Sess.U(i).name{j} '?'];
                if spm_input(str, 3, 'y/n', [1 0])
                    U.u             = [U.u Sess.U(i).u(33:end, j)];
                    U.name{end + 1} = Sess.U(i).name{j};
                    str             = 'Contrast weight';
                    U.w             = [U.w spm_input(str, 4, 'e', [], 1)];
                end
            end
        end
    else
        u      = length(Sess.U);
        U.u    = [];
        for  i = 1:u
            for  j = 1:length(Sess.U(i).name)
                if any(strcmp(Sess.U(i).name{j}, U.name))
                    U.u = [U.u Sess.U(i).u(33:end, j)];
                end
            end
        end
    end
end