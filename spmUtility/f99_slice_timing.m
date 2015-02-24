function f99_slice_timing(DIR,EXP,Seq,sl_order,refslice,TR,TA)

% '========================================================================='
%    'SLICE TIMING CORRECTION'
% '========================================================================='
% '               --- The Script Parameters Explained ---                   '


%  ------------------------------------------------------------------------
%  DIR       - [matrix containing directories of experiments]
%  EXP       - [matrix containing directories of images wildcard]
%  Seq       - slice acquisition order (1,2,3 = asc, desc, interl,user
%                                       specified)
%  sl_order  - chronological order of slice acquisition (if Seq = 4)
%  refslice  - slice for time 0
%  TR        - repitition time
%  TA        - acquisition time (time between start of acq of first
%              slice and the end of acq of the last slice)
%  ------------------------------------------------------------------------


%function spm_slice_timing(P, Seq, refslice, timing)
% function spm_slice_timing(P, Seq,refslice,timing)
% INPUT:
%       P               nimages x ?     Matrix with filenames
%       Seq             slice acquisition order (1,2,3 = asc, desc, interl)
%       refslice        slice for time 0
%       timing          additional information for sequence timing
%                       timing(1) = time between slices
%                       timing(2) = time between last slices and next volume
%
%       If no input is specified the function serves as a GUI
%
% OUTPUT:
%       None
%
% NMH_ACQCORRECT  Correct differences in image acquisition time between slices
%   This routine is intended to correct for the staggered order of
%   slice acquisition that is used during echoplanar scanning. The
%   correction is necessary to make the data on each slice correspond
%   to the same point in time. Without correction, the data on one
%   slice will represent a point in time as far removed as 1/2 the TR
%   from an adjacent slice (in the case of an interleaved sequence).
%
%   This routine "shifts" a signal in time to provide an output
%   vector that represents the same (continuous) signal sampled
%   starting either later or earlier. This is accomplished by a simple
%   shift of the phase of the sines that make up the signal.
%
%   Recall that a Fourier transform allows for a representation of any
%   signal as the linear combination of sinusoids of different
%   frequencies and phases. Effectively, we will add a constant
%   to the phase of every frequency, shifting the data in time.
%
%    Shifter - This is the filter by which the signal will be convolved
%    to introduce the phase shift. It is constructed explicitly in
%    the Fourier domain. In the time domain, it may be described as
%    an impulse (delta function) that has been shifted in time the
%    amount described by TimeShift.
%
%   The correction works by lagging (shifting forward) the time-series
%     data on each slice using sinc-interpolation. This results in each
%     time series having the values that would have been obtained had
%     the slice been acquired at the beginning of each TR.
%
%   To make this clear, consider a neural event (and ensuing hemodynamic
%     response) that occurs simultaneously on two adjacent slices. Values
%     from slice "A" are acquired starting at time zero, simultaneous to
%     the neural event, while values from slice "B" are acquired one
%     second later. Without corection, the "B" values will describe a
%     hemodynamic response that will appear to have began one second
%     EARLIER on the "B" slice than on slice "A". To correct for this,
%     the "B" values need to be shifted towards the Right, i.e., towards
%     the last value.
%
%   This correction assumes that the data are band-limited (i.e. there
%     is no meaningful information present in the data at a frequency
%     higher than that of the Nyquist). This assumption is support by
%     the study of Josephs et al (1997, NeuroImage) that obtained
%     event-related data at an effective TR of 166 msecs. No physio-
%     logical signal change was present at frequencies higher than our
%     typical Nyquist (0.25 HZ).
%
%   NOTE WELL:  This correction should be the first performed (i.e.,
%     before orienting, motion correction, padding, smoothing, etc.).
%     Additionally, it should only be performed once!
%
% Written by Darren Gitelman at Northwestern U., 1998
%
% Based (in large part) on ACQCORRECT.PRO from Geof Aquirre and
% Eric Zarahn at U. Penn.
%
% v1.0  07/04/98        DRG
% v1.1  07/09/98        DRG     fixed code to reflect 1-based indices
%                               of matlab vs. 0-based of pvwave
%
% Modified by R Henson, C Buechel and J Ashburner, FIL, 1999, to
% handle different sequence acquisitions, analyze format, different
% reference slices and memory mapping.
%
% Modified by M Erb, at U. Tuebingen, 1999, to ask for non-continuous
% slice timing and number of sessions.
%
% Modified by D. Gitelman 00/06/26 to adjust timing based on 
% chronological slice order.
% This is not as prone to error in specification of the reference 
% slice particularly for
% interleaved sequences.
%_______________________________________________________________________
% @(#)spm_slice_timing.m        2.7.1 00/06/26

global fmriDIR fmriTEST

if f99_exist(fmriDIR,'LOG') == 0
   str = ['!mkdir ' fmriDIR filesep 'LOG'];
   eval(str);
end
if f99_exist([fmriDIR filesep 'LOG'],'SLICETIME.OK')
   str = ['!rm -f ' fmriDIR filesep 'LOG' filesep 'SLICETIME.OK'];
   eval(str);
end

% Start logging
%--------------
logfile = [fmriDIR filesep 'LOG' filesep 'slicetime.log'];
tijd=spm('time');
disp('');
disp('*****************************************************************');
f99_log(logfile,['Slice Timing Script started at ' mat2str(tijd)]);
str = ['Start logging in ' logfile];
f99_log(logfile,str);

%SPMid = spm('FnBanner',mfilename,'2.6');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Slice timing');
spm_help('!ContextHelp',mfilename);

nsubjects = 1;
%if nargin < 1,
        % Choose the images
        %P = spm_get(+Inf,'*.img','Select images to acquisition correct');
% Modified by M Erb
         % get number of subjects
%         nsubjects = spm_input('number of subjects/sessions',1, 'e', 1);
%         if (nsubjects < 1)
%                 spm_figure('Clear','Interactive');
%                 return;
%         end
%         for i = 1:nsubjects
%                % Choose the images
%                P = [];
%                P = spm_get(+Inf,'*.img',...
%                        ['Select images to acquisition correct for subject ' num2str(i)]);
%                 eval(['P'    num2str(i) ' = P;']);
%         end
% end of Modified by M Erb
%end;

D = [];
for i = 1:size(DIR,1)
   d = [deblank(fmriDIR) filesep deblank(DIR(i,:))];
   D = [D; d];
end
DIR = D;

for s=1:size(DIR,1)
 P   = f99_P(DIR(s,:),EXP(s,:));
 nf  = size(P,1);
 str = ['  SESSION ' spm_str_manip(DIR(s,:),'t') ' with number of timepoints = ' mat2str(nf)];
 f99_log(logfile,str);
 f99_checkP(DIR(s,:),EXP(s,:),1);
end

P   = f99_P(DIR,EXP);

% map image files into memory
Vin     = spm_vol(P(1,:));
%nimgo   = size(P,1);
%nimg    = 2^(floor(log2(nimgo))+1);
nslices = Vin(1).dim(3);

%if nargin < 2,
%        % select fMRI acquisition sequence type
%        Stype = str2mat(...
%                'ascending (first slice=bottom)',...
%                'descending (first slice=top)',...
%                'interleaved (first slice=top)',...
%                'user specified');
%        str   = 'Select sequence type';
%        Seq   = spm_input(str,'!+1','m',Stype,[1:size(Stype,1)]);
%end;


% DRG added chronorder variable, modified reslice selection by time bin, and modified
% refslice calculation.
%--------------------------------------------------------------------
if Seq==[1],
    sliceorder = [1:1:nslices];
    chronorder = sliceorder;
elseif Seq==[2],
    sliceorder = [nslices:-1:1];
    chronorder = sliceorder;
elseif Seq==[3],
    % Assumes interleaved sequences top-middle downwards
    for k=1:nslices,
        sliceorder(k) = round((nslices-k)/2 + (rem((nslices-k),2) * (nslices - 1)/2)) + 1;
    end;
    chronorder = [1:2:nslices 2:2:nslices];
elseif Seq==[4],
%       sliceorder = [];
%       while length(sliceorder)~=nslices | max(sliceorder)>nslices | ...
%           min(sliceorder)<1 | any(diff(sort(sliceorder))~=1),
%           sliceorder = spm_input('Anatomical order of slices (1=bottom)','!+0','e');
%           chronorder = spm_input('Time order of slices','!+0','e');
%       end;
        sliceorder = sl_order;
        chronorder = sliceorder;
end;

fprintf(['The slice order is :\n']);
fprintf(['     ' mat2str(chronorder) '\n']);

%if nargin < 3,
%    % Choose reference slice in time format
%    % FIXED BY DRG ->NOW CHECKS. Note: no checking that 1 < refslice < no.slices (default = middle slice)
%    c = sprintf('[%i..%i]',1,nslices);
%    refslice = 0;
%    while ~any(1:nslices==refslice)
%       refslice = spm_input(['Reference slice by time bin 'c],'!+0','e',floor(Vin(1).dim(3)/2));
%    end
%end;

refslice = find(sliceorder == chronorder(refslice));
% end of change by DRG
 
%if nargin < 4,
%
% changed by M Erb
%       factor = 1/nslices;
%        TR = spm_input('Interscan interval (TR) {secs}','!+1','e',3);
%        TA = spm_input('Acquisition Time (TA) {secs}','!+1','e',TR);
%        TA = spm_input('Acquisition Time (TA) {secs}','!+1','e',TR-TR/nslices);
        while TA > TR | TA <= 0,
                TA = spm_input('Acquisition Time (TA) {secs}','!+0','e',TA);
        end;
        timing(2) = TR - TA;
        timing(1) = TA / (nslices -1);
        factor = timing(1)/TR;
% end of changed by ME
%
%else,
        TR      = (nslices-1)*timing(1)+timing(2);
        fprintf('Your TR is %1.1f\n',TR);
        factor = timing(1)/TR;
%end;


spm('Pointer','Watch')
%spm('FigName','Slice timing: working',Finter,CmdLine);


if ~fmriTEST

% Modified by M Erb
for subj = 1:nsubjects
    task=['Slice timing: working on session ' num2str(subj)];
    spm('FigName',task,Finter,CmdLine);
    %eval(['P    =    P' num2str(subj) ';']);
    if nsubjects > 1, eval(['P    =    P' num2str(subj) ';']); end;
    if ~fmriTEST
      Vin     = spm_vol(P);
    end
    nimgo   = size(P,1);
    nimg    = 2^(floor(log2(nimgo))+1);
    nslices_t= Vin(1).dim(3);
    if subj==1 
      nslices=nslices_t;
    end
    if ( nslices_t ~= nslices )
        fprintf('Number of slices differ! %d %\n', nimg);
    else
        % end of Modified by M Erb
        % create new header files
        Vout    = Vin;
        for k=1:nimgo,
            [pth,nm,xt,vr] = fileparts(deblank(Vin(k).fname));
            Vout(k).fname  = fullfile(pth,['a' nm xt vr]);
            if isfield(Vout(k),'descrip'),
                desc = [Vout(k).descrip ' '];
            else,
                desc = '';
            end;
            Vout(k).descrip = [desc 'acq-fix ref-slice ' int2str(refslice)];
            Vout(k) = spm_create_image(Vout(k));
        end;

        % Set up large matrix for holding image info
        % Organization is time by voxels
        slices = zeros([Vout(1).dim(1:2) nimgo]);
        stack  = zeros([nimg Vout(1).dim(1)]);
        
        spm_progress_bar('Init',nslices,'Correcting acquisition delay','planes complete');

        % For loop to read data slice by slice do correction and write out
        % In analzye format, the first slice in is the first one in the volume.
        for k = 1:nslices,

            % Set up time acquired within slice order
            shiftamount  = (find(sliceorder==k) - find(sliceorder==refslice)) * factor;

            % Read in slice data
            B  = spm_matrix([0 0 k]);
            for m=1:nimgo,
                slices(:,:,m) = spm_slice_vol(Vin(m),B,Vin(1).dim(1:2),1);
            end;

            % set up shifting variables
            len     = size(stack,1);
            phi     = zeros(1,len);

            % Check if signal is odd or even -- impacts how Phi is reflected
            %  across the Nyquist frequency. Opposite to use in pvwave.
            OffSet  = 0;
            if rem(len,2) ~= 0, OffSet = 1; end;

            % Phi represents a range of phases up to the Nyquist frequency
            % Shifted phi 1 to right.
            for f = 1:len/2,
                phi(f+1) = -1*shiftamount*2*pi/(len/f);
            end;

            % Mirror phi about the center
            % 1 is added on both sides to reflect Matlab's 1 based indices
            % Offset is opposite to program in pvwave again because indices are 1 based
            phi(len/2+1+1-OffSet:len) = -fliplr(phi(1+1:len/2+OffSet));

            % Transform phi to the frequency domain and take the complex transpose
            shifter = [cos(phi) + sin(phi)*sqrt(-1)].';
            shifter = shifter(:,ones(size(stack,2),1)); % Tony's trick

            % Loop over columns
            for i=1:Vout(1).dim(2),

                % Extract columns from slices
                stack(1:nimgo,:) = reshape(slices(:,i,:),[Vout(1).dim(1) nimgo])';

                % fill in continous function to avoid edge effects
                for g=1:size(stack,2),
                        stack(nimgo+1:end,g) = linspace(stack(nimgo,g),stack(1,g),nimg-nimgo)';
                end;

                % shift the columns
                stack = real(ifft(fft(stack,[],1).*shifter,[],1));

                % Re-insert shifted columns
                slices(:,i,:) = reshape(stack(1:nimgo,:)',[Vout(1).dim(1) 1 nimgo]);
            end;

            % write out the slice for all volumes
            for p = 1:nimgo,
                Vout(p) = spm_write_plane(Vout(p),slices(:,:,p),k);
            end;
            spm_progress_bar('Set',k);
        end; % loop over slices

        spm_progress_bar('Clear');
    nslices = nslices_t;
        % Modified by M Erb
        % endelse of "if ( nslices_t ~= nslices )"
    end
      
end % endfor of "for subj = 1:nsubjects"

% end of Modified by M Erb

spm('FigName','Slice timing: done',Finter,CmdLine);
spm('Pointer');

end % if ~fmriTEST

statusfile = [fmriDIR filesep 'LOG' filesep 'SLICETIME.OK'];
f99_log(statusfile,'');

tijd=spm('time');
f99_log(logfile,['Slice Timing Script ended at ' mat2str(tijd)]);
disp('-----------------------------------------------------------------');
disp('');
