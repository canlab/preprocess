function sdt2an(exptname,runs,outpath)
%
% sdt2an.m
% Convert a file from Stimulate to Analyze format.
% (Copy .sdt to .img, add a header.)
%
% by Jack Grinband
% modified to make function by Tor Wager
%
% Inputs: 
%       exptname = experiment subject code in /export/data/reconstruction
%       runs = string with names of runs r1 ... r*, e.g., 'r1 r2 r3'
%       outpath = output path to save data in

%-----------------
% Get input
%-----------------

stimfile   = 'stimulate.sdt';
%exptname   = from shell
numruns = length(find(runs == 'r'));  %input('How many expts/runs?  ');
inpath   = '/export/data/reconstruction/';

disp(['Saving 4-D images in: ' outpath])

%outpath   = mdir;   %input('Output directory? (Ex: /home/jgrinband/data/tst/):  ' ,'s');

%------------------
% Loop through runs
%------------------
for run_num  = 1 : numruns
    s=sprintf('Run number %d of %d', run_num,numruns); disp(s);
    
    reconfile=sprintf('%s%s/r%d/recon.log',inpath,exptname,run_num);
    disp('Reading paramaters...');
    fid = fopen(reconfile, 'r');
    header=fscanf(fid,'%s');
    
    temp='FinalReconstructionXres=';
    tempstart=findstr(header,temp);
    dimx=sscanf(header(tempstart+length(temp):tempstart+length(temp)+7),'%f')
    
    temp='Yres=';
    tempstart=findstr(header,temp);
    dimy=sscanf(header(tempstart(2)+length(temp):tempstart(2)+length(temp)+7),'%f')
    
    temp='SliceThickness=';
    tempstart=findstr(header,temp);
    dimz=sscanf(header(tempstart+length(temp):tempstart+length(temp)+7),'%f')
    
    temp='FOV(x)=';
    tempstart=findstr(header,temp);
    fovx=sscanf(header(tempstart+length(temp):tempstart+length(temp)+7),'%f')
    
    temp='FOV(y)=';
    tempstart=findstr(header,temp);
    fovy=sscanf(header(tempstart+length(temp):tempstart+length(temp)+7),'%f')
    
    temp='NumberEPIReps=';
    tempstart=findstr(header,temp);
    nvols=sscanf(header(tempstart+length(temp):tempstart+length(temp)+7),'%f')
    
    temp='NumberofSlicesperRep=';
    tempstart=findstr(header,temp);
    nslices=sscanf(header(tempstart+length(temp):tempstart+length(temp)+7),'%f')
    
    
    stim_file = [inpath,exptname,'/r',int2str(run_num),'/',stimfile];
    analyze_file  = [outpath,'/r',int2str(run_num),'.img'];
    
    disp('Copying sdt file...');
    s=sprintf('unix(''cp %s %s'');',stim_file,analyze_file); eval(s);
    header_file  = [outpath,'/r',int2str(run_num)];

    disp('Creating an Analyze header...');
    s=sprintf('unix(''avwcreatehd %d %d %d %d %d %d %d %d %d %d %d %d %s'');',dimx,dimy,nslices,nvols,fovx/dimx,fovy/dimy,dimz,0,0,0,0,4,header_file); 
    eval(s);
    
    disp('Flipping along x and y i.e. making neurological coordinates, L=L and R=R');
    s=sprintf('unix(''avwswapdim %s -x -y z %s'');',header_file,header_file);
    eval(s);  
 
end


disp('Voila!')



 
 
 
 
