% sdt2analyze(stim_file, output_base, [recon_file])
%
% This function converts Stimulate (.sdt) files to Analyze format (*.hdr/*.img).
% 
% Optional inputs:
%   recon_file - if not specified, defaults to 'recon.log'
%
% Sample usage: 
% sdt2analyze('stimulate.sdt', 'r3') % to convert the stimulate.sdt file to r3.img and r3.hdr
%
% or
%
% subjdirs = {'rea24' 'rea25' 'rea26'};
% rundirs = {'r1' 'r2' 'r3' 'r4' 'r5' 'r6'};
% for i=1:length(subjdirs)
%     for j=1:length(rundirs)
%         sdt2analyze(fullfile(subjdirs{i}, rundirs{j}, 'stimulate.sdt'), rundirs{j});
%     end
% end

function sdt2analyze(stim_file, output_base, recon_file)
    global FSLDIR;

    scn_setup();

    if(~exist('recon_file', 'var') || isempty(recon_file))
        recon_file = 'recon.log';
    end

    [path, filename, ext] = fileparts(stim_file);

    recon_file = fullfile(path, recon_file);
    if(~exist(recon_file, 'file'))
        error(['Unable to locate ' recon_file]);
    end

    disp('Reading parameters...');
    fid = fopen(recon_file, 'r');
    header = fscanf(fid, '%s');

    param = 'FinalReconstructionXres=';
    paramstart = strfind(header, param);
    dimx = sscanf(header(paramstart+length(param):paramstart+length(param)+7), '%f');

    param = 'Yres=';
    paramstart = strfind(header, param);
    dimy = sscanf(header(paramstart(2)+length(param):paramstart(2)+length(param)+7), '%f');

    param = 'SliceThickness=';
    paramstart = strfind(header, param);
    dimz = sscanf(header(paramstart+length(param):paramstart+length(param)+7), '%f');

    param = 'FOV(x)=';
    paramstart = strfind(header, param);
    fovx = sscanf(header(paramstart+length(param):paramstart+length(param)+7), '%f');

    param = 'FOV(y)=';
    paramstart = strfind(header, param);
    fovy = sscanf(header(paramstart+length(param):paramstart+length(param)+7), '%f');

    param = 'NumberEPIReps=';
    paramstart = strfind(header, param);
    nvols = sscanf(header(paramstart+length(param):paramstart+length(param)+7), '%f');

    param = 'NumberofSlicesperRep=';
    paramstart = strfind(header, param);
    nslices = sscanf(header(paramstart+length(param):paramstart+length(param)+7), '%f');

    param = 'TR=';
    paramstart = strfind(header, param);
    tr = sscanf(header(paramstart+length(param):paramstart+length(param)+6), '%f') / 1000;


    header_file = fullfile(path, sprintf('%s.hdr', output_base));
    image_file = fullfile(path, sprintf('%s.img', output_base));

    disp('Creating an Analyze header...');
    header_cmd = sprintf('. %s/etc/fslconf/fsl.sh && %s/bin/avwcreatehd %d %d %d %d %g %g %g %d %d %d %d %d %s', ...
        FSLDIR, FSLDIR, dimx, dimy, nslices, nvols, fovx/dimx, fovy/dimy, dimz, tr, 0, 0, 0, 4, header_file);
    unix(header_cmd);

    % Must come after creating the header, as avwcreatehd (in v3.2)
    % generates an empty .img file when executed
    fprintf('Copying %s to %s...\n', stim_file, image_file);
    copyfile(stim_file, image_file, 'f');

    % The axis input to avwswapdim seems weird because it is. It makes
    % certain assumptions about how you want things if you're generating
    % Analyze images. Rest assured it works (as of 2005.8.5).
    % As insane as it sounds, rotating around y and not both x and y will spit
    % out an image in neurological space.
    disp('Reordering data to neurological orientation...');
    flip_cmd = sprintf('. %s/etc/fslconf/fsl.sh && %s/bin/avwswapdim %s x -y z %s', FSLDIR, FSLDIR, image_file, image_file);
    unix(flip_cmd);

    fclose(fid);
end