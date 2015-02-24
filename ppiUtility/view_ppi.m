function view_ppi(PPI)
    check_inputs();

    N = length(PPI.xY(1).u);
    RT = length(PPI.xn) / (length(PPI.Y) / PPI.dt);
    t = RT * (1:N);
    NT = RT / PPI.dt;
    T = PPI.dt * (1:(N*NT));

    PSY = gen_PSY();

    subplot(2, 1, 1);
    plot(t, PPI.Y, T, PPI.xn(:,1));
    title('hemodynamic and neuronal responses');
    xlabel('time (secs)');
    axis tight square
    grid on
    legend('BOLD', 'neuronal');

    subplot(2, 2, 3);
    plot(t, PPI.P, T, PSY, '--');
    title('[convolved] psych. variable');
    xlabel('time (secs)');
    axis tight square
    grid on
    legend('convolved', 'orig');

    subplot(2, 2, 4);
    plot(t, PPI.ppi);
    title('PPI');
    xlabel('time (secs)')
    axis tight square
    grid on


    function check_inputs()
        if(ischar(PPI))
            if(exist(PPI, 'file'))
                load(PPI);
            else
                error('Unable to locate PPI file: %s\n', PPI);
            end
        elseif(~isstruct(PPI))
            error('Unknown format for PPI variable. Not a filename or structure.');
        end
    end
    
    function PSY = gen_PSY()
        PSY = zeros(N*NT, 1);
        for i = 1:size(PPI.psy.u, 2)
            PSY = PSY + full(PPI.psy.u(:,i)*PPI.psy.w(:,i));
        end
        PSY = spm_detrend(PSY);
    end
end
