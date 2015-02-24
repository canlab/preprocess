function c = spm_config_gsrpreproc(varargin)

addpath(fullfile(spm('dir'),'toolbox','GSRtoolbox'));



data.type = 'files';
data.name = 'Raw GSR signal';
data.tag  = 'gsrsig';
data.num  = [1 Inf];
data.filter = '.*';
data.help = {'Select the signals to preprocess'};

trialdat.type = 'files';
trialdat.name = 'ASCII files with trial onset data';
trialdat.tag  = 'trialdata';
trialdat.num  = [1 Inf];
trialdat.filter = '.*';
trialdat.help = {'Select file with trial information, or files if you have multiple sessions.', ...
    'If fewer trial files are specified than signal files, trial information will be duplicated to compensate.', ... 
    'Onsets should be entered in a column vector.', ...
    'The columns in the file should be:', ...
    'Onsets (sec) Offsets (sec) Integer Code for Trial Type', ...
    };

sigcol.type = 'entry';
sigcol.name = 'Signal column';
sigcol.tag  = 'sigcol';
sigcol.strtype = 'n';
sigcol.num  = [1 1];
sigcol.help = {'Enter that column number in the text file that corresponds to the GSR signal'};


trigcol.type = 'entry';
trigcol.name = 'Trigger column in signal data';
trigcol.tag  = 'trigcol';
trigcol.strtype = 'w';
trigcol.num  = [1 1];
trigcol.help = {'Enter that column number in the text file that corresponds to the event triggers', ...
    'Enter 0 if non-existent.'};

freq.type = 'entry';
freq.name = 'Acquisition frequency';
freq.tag  = 'fs';
freq.strtype = 'n';
freq.num  = [1 1];
freq.help = {'Frequency of acquisition', ...
    'If onsets are given in a unit other than seconds, adjust this value accordingly'};


choice.type = 'choice';
choice.name = 'Trial format';
choice.tag = 'choice';
choice.values = {trigcol, trialdat};
choice.help = {'Choose trial onset format'};


c.type = 'branch';
c.name = 'GSR preprocessing';
c.tag  = 'gsrpreproc';
c.val  = {data, sigcol, freq, choice};
c.prog = @spm_gsrpreproc;
c.help = {'Combine signal and trial trigger information.',...
'Requires ascii text input.'};


end
%_______________________________________________________________________

%_______________________________________________________________________
function spm_gsrpreproc(job)
    
    spm_progress_bar('init');
    
    if isfield(job.choice,'trialdata') 
        if length(job.choice.trialdata)<length(job.gsrsig)
            dups = ceil(length(job.choice.trialdata)/length(job.gsrsig));
            for i=1:dups
                evalstr = [evalstr 'job.choice.trialdata; '];
            end;
            eval(['job.choice.trialdata = [' evalstr '];']);
        end;
      
    end;
    
    
    
    
    for i=1:length(job.gsrsig)
        fname = job.gsrsig{i};
        disp('Loading data.')
        dat = load(fname);
        
        n = length(dat(:,job.sigcol));
        
        outsig = dat(:,job.sigcol);
        if ~isfield(job.choice,'trialdata')
            if job.choice.trigcol ~= 0
                %outdat = [outdat dat(:,job.choice.trigcol)];
                %convert
                
                tdat = find(dat((1:n-1),job.choice.trigcol)~=dat((2:n),job.choice.trigcol));
                if dat(1,job.choice.trigcol)
                    tdat = [0; tdat];
                end;
                if dat(n,job.choice.trigcol)
                    tdat = [tdat; n];
                end;
                tdat = [tdat(1:2:length(tdat))+1 tdat(2:2:length(tdat))];
                tdat = [tdat dat(tdat(:,1),job.choice.trigcol)];
                
            else
                tdat = [1 n 1];
                %outdat = [outdat ones(n,1)];
            end;
        else
            % we do have 'trialdata' field with filename for onsets
            % Load onsets from ASCII file
            
            tdat = load(job.choice.trialdata{i});
            if size(tdat,2) == 1
                % We have only onsets; Generate offsets and trial types
                dur = input('Cannot find durations in file. Enter duration (in sec) for all trial types: ');
                offsets = dur + tdat;
                trialtypes = ones(size(tdat, 1), 1);
                
                tdat = [tdat offsets trialtypes];
                
            elseif size(tdat,2)==2
                % We have  onsets and offsets; Generate trial types
                trialtypes = ones(size(tdat, 1), 1);
                tdat=[tdat trialtypes];
                
            elseif size(tdat,2)==4
                disp('WARNING!!!! THIS BIT OF CODE NOT UPDATED!!!');
                tdatn = size(tdat,1);
                tdatch = [1; find(tdat(1:tdatn-1,1)~=tdat(2:tdatn,1))];
                tdat = tdat(find(tdat(:,1)==tdat(tdatch(i),1)),2:4);
            
            end;    
                
            % Convert to samples, assuming input in sec    
            tdat(:,1:2) = tdat(:,1:2)*job.fs;
            
            % print your onsets
            fprintf('Your event data in samples (tab delimited):\n');
            print_matrix(tdat, {'Onsets' 'Offsets' 'Conditions'});
            
        end 
        
        dirdelim = find(fname==filesep);
        fdir = fname(1:dirdelim(length(dirdelim)));
        fname = fname(dirdelim(length(dirdelim))+1:length(fname)-3);
        
        output.signal = outsig;
        output.trialonsets = tdat;
        output.fs = job.fs;
        
        fprintf('Saving: %s\n', [fdir fname 'mat']);
        save([fdir fname 'mat'], 'output');
        
        
        spm_progress_bar('set', i/length(job.gsrsig));
    end 
    spm_progress_bar('clear');
end


