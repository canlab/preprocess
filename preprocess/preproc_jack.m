exptName=input('Experiment name:  ','s');
disp('  ');
disp('Enter experiment path')
expPath=input('(e.g. /home/jgrinband/data/decide/decide6  ):  ','s');

numRuns=input('Enter number of runs to preprocess:  ');
inPath   = '/export/data/reconstruction/';

% deleteVol=input('Enter number of volumes to delete');
deleteVol=3;

% disp('Enter path for the high res image');
% hiRes=input('(e.g. /home/jgrinband/data/jg/spgr1/spgr_bet   )','s');
hiRes='/home/jgrinband/data/jg/spgr1/spgr_bet';
standard='/export/apps/fsl_3.1/etc/standard/avg152T1_brain.hdr';

for i=1:numRuns
    filename=[expPath '/r' num2str(i)];
    reconFile=sprintf('%s%s/r%d/recon.log',inPath,exptName,i)
    disp('Reading paramaters...');
    fid = fopen(reconFile, 'r');
    header=fscanf(fid,'%s');
    
    % search for TR
    temp='TR=';
    tempstart=findstr(header,temp);
    TR=sscanf(header(tempstart+length(temp):tempstart+length(temp)+7),'%f')
    
    % search for number of volumes 
    temp='NumberEPIReps=';
    tempstart=findstr(header,temp);
    totalVols=sscanf(header(tempstart+length(temp):tempstart+length(temp)+7),'%f')
    
    numVols=totalVols-deleteVol
    
    s=sprintf('unix(''mkdir %s/'');',filename);eval(s);
    s=sprintf('unix(''mkdir %s/mc'');',filename);eval(s);
    
    % slice timing correction 
    disp('slice timing correction - interleaved (1,3,5,...,2,4,6,...)')
    s=sprintf('unix(''slicetimer -i %s -o %s_stc -r %d --odd --verbose'');',filename,filename,TR);
    eval(s);  
    
    % motion correction 
    disp('motion correction - aligning to 1st high contrast image - using mutual information')
    s=sprintf('unix(''mcflirt -in %s_stc -o %s_mc -cost mutualinfo -refvol 0 -stats -plots -rmsrel -rmsabs -report'');',filename,filename);
    eval(s);  
    
    
    
    % s=sprintf('unix(''echo "set term pbm color ; set size 1,0.3 ; set title ''MCFLIRT estimated rotations (radians)'' ; plot ''prefiltered_func_data_mcf.par'' using 1 title ''x'' with lines , ''prefiltered_func_data_mcf.par'' using 2 title ''y'' with lines , ''prefiltered_func_data_mcf.par'' using 3 title ''z'' with lines" | gnuplot > rot.ppm'');');
    % eval(s);
    % s=sprintf('unix(''echo "set term pbm color ; set size 1,0.3 ; set title ''MCFLIRT estimated translations (mm)'' ; plot ''prefiltered_func_data_mcf.par'' using 4 title ''x'' with lines , ''prefiltered_func_data_mcf.par'' using 5 title ''y'' with lines , ''prefiltered_func_data_mcf.par'' using 6 title ''z'' with lines" | /export/apps/fsl_3.1/bin/gnuplot > trans.ppm'');');
    % eval(s);
    % s=sprintf('unix(''echo "set term pbm color ; set size 1,0.3 ; set title ''MCFLIRT estimated mean displacement (mm)'' ; plot ''prefiltered_func_data_mcf_abs.rms'' title ''absolute'' with lines , ''prefiltered_func_data_mcf_rel.rms'' title ''relative'' with lines" | /export/apps/fsl_3.1/bin/gnuplot > disp.ppm'');');
    % eval(s);
    % s=sprintf('unix(''convert rot.ppm rot.gif;convert trans.ppm trans.gif;convert disp.ppm disp.gif'')');
    % eval(s);
    %     
    
    %    s=sprintf('unix(''mv -f %s_mc* %s/mc'');',filename,filename);
    %    eval(s);  
    
    %;rm -f rot.ppm trans.ppm disp.ppm
    
    % save first volume 
    disp('generating first volume for later registration')
    s=sprintf('unix(''avwroi %s_mc %s_first 0 1'');',filename,filename);
    eval(s); 
    
    % extract brain from first volume and create mask
    disp('bet - extracting brain, f = 0.2')
    s=sprintf('unix(''bet %s_first %s_first -m -f .2 '');',filename,filename);
    eval(s); 
    
    % multiply all t2* images by brain mask to extract
    disp('multiplying t2* images by brain mask')
    s=sprintf('unix(''avwmaths %s_mc -mul %s_first_mask %s_mc_brain'');',filename,filename,filename);
    eval(s); 
    
    % copy 1st vol to example_func
    disp('copying first volume to example_func for later registration')
    s=sprintf('unix(''cp %s_first.img %s/example_func.img'');',filename,filename);
    eval(s); 
    s=sprintf('unix(''cp %s_first.hdr %s/example_func.hdr'');',filename,filename);
    eval(s); 
    
    % delete first 3 vols   
    disp('deleting first 3 volumes')
    s=sprintf('unix(''avwroi %s_mc_brain %sfil 3 %d'');',filename,filename,numvols);
    eval(s); 
    
    % phasemap unwarping of first vol
    
    % register first vol to hiRes SPGR
    disp('registering first vol')
    s=sprintf('unix(''flirt -in %s_first -ref %s -out %s_first2hires.hdr -omat %s_first2hires.mat -cost normmi -dof 7'');',filename,hiRes,filename,filename);
    eval(s);  
    
    % register hiRes SPGR to MNI152
    s=sprintf('unix(''flirt -in %s -ref %s -out %s_hires2std.hdr -omat %s_hires2std.mat -cost normmi -dof 12'');',hiRes,standard,filename,filename);
    eval(s);  
    
    % combine the registration matrices
    s=sprintf('unix(''convert_xfm -matonly -omat %s_first2std -concat %s_hiRes2std %s_first2hiRes '');',filename,filename,filename);
    
    
    s=sprintf('unix(''mv %s*.* %s/'');',filename,filename);
    eval(s);  
    
    s=sprintf('unix(''mv %s/*mc_* %s/mc/'');',filename,filename);
    eval(s);  
    
    
    %     s=sprintf('unix(''rm %s/r%d_* '');',filename,i);
    %     eval(s);  
    
end
