function outfiles = tor_mr2an(subjid,output_file_base)
% output_file_names = tor_mr2an(subjid,output_file_base)
%
% subjid:   ID code for subject, e.g., ow1
%   disp('You type:  /export/data/scanner/gmr1/e264/s2/E264S2I   ')
% output_file_base: dir to write output to
%   e.g., /export/data/analysed/tor/Class/
%   This program appends your subjid and 'anatomy' to the path:
%   e.g., /export/data/analysed/tor/Class/ow1/anatomy
%
% by Jack Grinband, modified 2 / 22 / 04 by Tor Wager to function-alize.
% Tor also added smoothing, brain extraction, and segmentation with FSL
% tools (only for structurals with more than 18 slices).
% Also fixed close file bug:  needed call to fclose
%
% examples:
% tor_mr2an('ow1','/export/data/analyzed/tor/Class');
%
% batch:
% for s = {'ow1' 'ow2' 'ow3' 'ow4' 'ow5' 'ow6'}            
%   out = tor_mr2an(s{1},'/export/data/analyzed/tor/Class')  
% end

% Input path stuff

ifb = fullfile('/export/data/scanner/', subjid);
tmp = dir([ifb filesep 'e*']);
eval(['cd ' fullfile(ifb,tmp(1).name)])
ifb = pwd;

% Output path stuff (later, add struct name)

ofb = fullfile(filesep,output_file_base,subjid);
if ~(exist(ofb) == 2), eval(['!mkdir ' ofb]);, end
ofb = fullfile(ofb,'anatomy');
if ~(exist(ofb) == 2), eval(['!mkdir ' ofb]);, end

diaryname = [ofb filesep 'stacklog.txt'];
disp(['tor_mr2an.m: Logging stack in ' diaryname])
diary(diaryname)

output_file_ext='.img';

d = dir('s*');  % find all s* subdirectories

for myimage = 1:length(d)
    
    dd = dir(d(myimage).name); [dd,ff,input_file_ext] = fileparts(dd(3).name);
    input_file_base = fullfile(ifb,deblank(d(myimage).name),ff);
    output_file_base = fullfile(ofb,['T1_' d(myimage).name]);
    
    fprintf(1,'mr2an.m: Working on %s\n',input_file_base); 
    
% Jack's code
% ------------------------

lastslash=strfind(output_file_base,'/');
short_output_file_base=output_file_base(1:lastslash(end));
outputfilename=sprintf('%s%s',output_file_base,output_file_ext);

s=sprintf('unix(''lximg %s.MR>%sheader.txt'');',input_file_base,short_output_file_base); 
eval(s);
fid = fopen([short_output_file_base 'header.txt'], 'r');
header=fscanf(fid,'%s');

temp='ImageSliceThickness:';
tempstart=findstr(header,temp);
zstep=sscanf(header(tempstart+length(temp):tempstart+length(temp)+7),'%f');

temp='Numberofslices:';
tempstart=findstr(header,temp);
topslice=sscanf(header(tempstart+length(temp):tempstart+length(temp)+4),'%f');

temp='ImageXDimension:';
tempstart=findstr(header,temp);
dimx=sscanf(header(tempstart+length(temp):tempstart+length(temp)+4),'%f');

temp='ImageYDimension:';
tempstart=findstr(header,temp);
dimy=sscanf(header(tempstart+length(temp):tempstart+length(temp)+4),'%f');

temp='DisplayFOV-X:';
tempstart=findstr(header,temp);
fov=sscanf(header(tempstart+length(temp):tempstart+length(temp)+15),'%f');

fprintf(1,'zstep=%3.2f topslice=%3.0f dimx=%3.0f dimy=%3.0f fov=%3.0f\n',zstep,topslice,dimx,dimy,fov);

botslice=1;

% Remove old output file


s=sprintf('unix(''rm %s*'');', output_file_base); 
%fprintf(1,'Removing old files: %s',s); eval(s);

% Start cat loop

%disp('Concatenating files...');

for i = botslice : topslice
    inputfilename = sprintf('%s%d%s',input_file_base(1:end-1),i,input_file_ext);
    fid = fopen(inputfilename);     % open input file
    
    if fid == -1
        disp('OOPS - tried to open image file and failed')
        disp(['File is: ' inputfilename])
        disp('If 2 below, the file exists.')
        exist(inputfilename,'file')
    else
        fseek(fid,8432,'bof');
        dat= fread(fid);%,'int16');%,[1 inf], 'single'); ''                              % read input file into F
    end
        
    fclose(fid);
    
    if i == botslice                                                              % MRdat is the array into which data is concatenated
        MRdat = dat;                                                              % Define MRdat
    else
        MRdat = cat(1,MRdat,dat);                                                 % Append dat to MRdat
    end       
end

fprintf(1,'Saving structural (%3.0f slices) in %s',i,outputfilename);

fid = fopen(outputfilename,'w');
fwrite(fid,MRdat,'uint8');
fclose(fid);

% Create Analyze header

%disp('Creating an Analyze header for T1 volume...');
s=sprintf('unix(''avwcreatehd %d %d %d %d %d %d %d %d %d %d %d %d %s'');',dimx,dimy,topslice,1,fov/dimx,fov/dimy,zstep,0,0,0,0,4,output_file_base); 
eval(s);

disp('Flipping along x and y i.e. making neurological coordinates, L=L and R=R');
s=sprintf('unix(''avwswapdim %s -x -y z %s'');',output_file_base,output_file_base);
eval(s);  

% end Jack's code

if myimage == 1, outfiles(1,:) = outputfilename;,
else, outfiles = str2mat(outfiles, outputfilename);
end
  
% Brain extraction and segmentation
%/Applications/fsl/bin/susan_smooth /Users/tor/Desktop/ow1_anatomy/T1_s3.hdr 52.2454 /Users/tor/Desktop/ow1_anatomy/T1_s3_susan.hdr 0 3D 1 0
%/Applications/fsl/bin/bet /Users/tor/Desktop/ow1_anatomy/T1_s5 /Users/tor/Desktop/ow1_anatomy/T1_s5_brain -f 0.5 -g 0
% /Applications/fsl/bin/fast -t1 -c 3 -n -os -od /Users/tor/Desktop/ow1_anatomy/T1_s5_susan_brain /Users/tor/Desktop/ow1_anatomy/T1_s5_susan_brain.hdr

if topslice > 18    % do this stuff if we have a full image, otherwise don't bother
    
    % 'Nonlinear' noise reduction
    
    [dd,ff,ee] = fileparts(outputfilename);
    V = spm_vol(outputfilename);v = spm_read_vols(V);
    nn = std(v(:)) ./ 3;
    sout = fullfile(dd,['s' ff ee]);
    s = ['!susan_smooth ' outputfilename ' ' num2str(nn) ' ' sout ' 0 3D 1 0'];
    disp(s)
    eval(s);
    
    % brain extraction
    
    [dd,ff,ee] = fileparts(sout);
    bout = fullfile(dd,['e' ff ee]);
    s = ['!bet ' sout ' ' bout ' -f 0.5 -g 0'];
    disp(s)
    eval(s);
    
    % segmentation
    [dd,ff,ee] = fileparts(sout);
    fout = fullfile(dd,['e' ff ee]);
    s = ['!fast -t1 -c 3 -n -os -od ' sout(1:end-4) ' ' sout(1:end-4) '.hdr'];
    disp(s)
    eval(s);
    
end
    
end % end loop through structurals

% create tar.gz file of structurals
cwd = pwd;
cd(ofb)
s = ['!gtar zcf ' subjid '_structurals.tar.tz *'];
eval(s);
cd(cwd)

disp('Voila!')
diary off


