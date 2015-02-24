load('../analysis/rating');
load('../analysis/version_info');

diary Physio_correlations_output.txt

TR = 2;
scanlen = 540;
numTR = scanlen./TR;

%subphys = {'cei1' 'cei10' 'cei11' 'cei12' 'cei13' 'cei14' 'cei15' 'cei16' 'cei17' 'cei18' 'cei19' 'cei2' 'cei20' 'cei21' 'cei22' 'cei4' 'cei5' 'cei6' 'cei7' 'cei8'};

subphys = {'cei1' 'cei2' 'cei4' 'cei5' 'cei6' 'cei7' 'cei8' 'cei10' 'cei11' 'cei12' 'cei13' 'cei14' 'cei15' 'cei16' 'cei17' 'cei18' 'cei19'  'cei20' 'cei21' 'cei22' };


% subidxs = 11:length(subphys);

subidxs = 1:length(subphys);
%subidxs = [12];

if(length(subidxs) == length(subphys))
    physio = [];
    Xphys = [];
end

subjects = fieldnames(rating);

for i = subidxs

    %subj = ['rea' subjects{i}(2:end)];

    subj = subphys{i};      % make sure they match !  rename so they have same names!
    subj2 = subjects{i};

    fprintf(1,'THESE SHOULD MATCH: %s  %s',subj,subj2);

    if(~(exist([subj '.txt']) == 2))
        disp(['No file: ' subj]);
    else

        try

            physio = read_physio_data({[subj '.txt']},physio);

            % loop thru each run
            [scanner_onperiods] = physio_get_scanperiods(physio.(subj).scanner_pulse, scanlen);

            if isempty(scanner_onperiods)
                disp('Cannot find scan period of correct length!!');
                [scanner_onperiods,allon] = physio_get_scanperiods(physio.(subj).scanner_pulse, scanlen,'doplot',1);
                allon
                scanner_onperiods = input(['Enter scanner on and off times in []: ']);
            end
            
            %plot_physio(physio, subj,scanner_onperiods(1,1),scanner_onperiods(1,2));


            % limit all vars to scanner on only
            N = fieldnames(physio.(subj));
            N(strcmp(N,'scanner_pulse')) = [];
            for j = 1:length(N)
                physio.(subj).(N{j}) = physio.(subj).(N{j})(scanner_onperiods(1):scanner_onperiods(2));
            end


            % get pulse
            % ---------------
            physio.(subj) = physio_pulse2hr(physio.(subj),0);
            %physio.(subj).hr_tr = resample(physio.(subj).hrs,1,TR);  % sampled in TRs
            % resample produces big edge effects!!!
            disp('NOTE: This program requires work at this point for TRs that are fractional!!!');
            physio.(subj).hr_tr = physio.(subj).hrs(1:2:end);  % sampled in TRs
            physio.(subj).hr_tr = physio.(subj).hr_tr(1:numTR);

            % get gsr
            % ---------------
            gsr = physio.(subj).gsr;
            gsr = [gsr(1:(100*TR*10)); gsr; gsr(end-(100*TR*10):end)];
            gsr = resample(gsr,1,100*TR);
            gsr = gsr(11:end-10);
            gsr = gsr(1:numTR);
            physio.(subj).gsr_tr = gsr;

            % get rating
            % ---------------
            if (subjversion.(subj) == 1|| subjversion.(subj) == 3)
                r = [rating.(subj2).Baseline.rating rating.(subj2).Panelists.rating rating.(subj2).Computer.rating rating.(subj2).Recovery.rating];
                o = [rating.(subj2).Baseline.onset rating.(subj2).Panelists.onset rating.(subj2).Computer.onset rating.(subj2).Recovery.onset];
            elseif (subjversion.(subj) == 2|| subjversion.(subj) == 4)
                r = [rating.(subj2).Baseline.rating rating.(subj2).Computer.rating rating.(subj2).Panelists.rating rating.(subj2).Recovery.rating];
                o = [rating.(subj2).Baseline.onset rating.(subj2).Computer.onset rating.(subj2).Panelists.onset rating.(subj2).Recovery.onset];
            else
                disp('ERROR IN VERSION INFO!!!'); 
                error('bad!!!')
            end
            
            o = o./(1000*TR);

            % pad end points
            r = [r(1) r r(end)];
            o = [0 o numTR+1];

            % interpolate
            r2 = interp1(o',r,(1:numTR)');

            physio.(subj).rating_tr = r2;
            physio.(subj).rating = r;
            physio.(subj).ratetimes_in_trs = o;

            %figure;plot(r2)
            %hold on; plot(o,r,'ro','MarkerFaceColor','r')

            physio.(subj).covs = [physio.(subj).hr_tr  physio.(subj).gsr_tr physio.(subj).rating_tr];

            % FIGURE 1
            
            tor_fig(3,1);
            nms = {['Heart Rate ' subj ' ' subj2] 'GSR' 'Emotion ratings'};
            for j=1:3
                subplot(3,1,j)
                plot(physio.(subj).covs(:,j));
                title(nms)
            end
            
            set(gcf,'Position',[131   300   441   606])
            drawnow
            try
                saveas(gcf,[subj],'fig');
                saveas(gcf,[subj],'tif');
            catch
                disp('Error saving figure.')
            end

            Xphys.(subphys{i}) = linear_detrending(physio.(subj).covs')';

            % FIGURE 2
            tor_fig;
            nms = {'Heart Rate' 'GSR' 'Emotion ratings'};
            plot((1:numTR)*2,scale(Xphys.(subphys{i})));
            legend(nms,7); ylabel('Standardized values'), xlabel('Time (s)');
            try
                saveas(gcf,[subj '_2'],'fig');
                saveas(gcf,[subj '_2'],'tif');
            catch
                disp('Error saving figure.')
            end
            
            
            Xphys.correls(:,:,i) = corrcoef(Xphys.(subphys{i}));
            print_matrix(squeeze(Xphys.correls(:,:,i)));

            [nxc,nxl] = shift_correl_matrix(Xphys.(subphys{i}),6,0,0);
            Xphys.crosscorrels(:,:,i) =nxc;
            Xphys.crosslatency(:,:,i) =nxl;
            
        catch
            disp('Error processing!')
        end

    end % if phys file exists

end

[Xphys.meancor,Xphys.tcor,Xphys.sigcor] = ttest3d(Xphys.correls);
disp('Mean correlation:')
print_matrix(Xphys.meancor)
disp('t correlation:')
print_matrix(Xphys.tcor)
disp('sig correlation:')
print_matrix(Xphys.sigcor)

[Xphys.meanxcor,Xphys.txcor,Xphys.sigxcor] = ttest3d(Xphys.crosscorrels);
disp('Mean cross-correlation:')
print_matrix(Xphys.meanxcor)
disp('t correlation:')
print_matrix(Xphys.txcor)
disp('sig correlation:')
print_matrix(Xphys.sigxcor)

diary off



