function [dat,cdat,ncond,cmat,names]=geteegtest


%%%%load EEG data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mea_path='\\gamma\chris\cross\beast\'
load([mea_path,'distance1']);
ncond=3;
nsub=19;

%%%%reduce to 20 electrode configuration & reshape
%load('elecmat');
load('elecmat20');
nchan=20;
%%%%%elecmat12=findelec('AF3FpzAFzNz T7 T9 CP5FT9PO3Oz POzPO4')
names={'Oz' 'Pz' 'FCz' 'Fpz' 'PO3' 'CP3' 'FC3' 'AF3' 'PO4' 'CP4' 'AF4' 'FC4' 'P7' 'T7' 'F7' 'P8' 'T8' 'F8' 'CB1' 'CB2'};
dat=d.dwg(:,:,:,:);
%names={'AF3' 'Fpz' 'AFz' 'Nz' 'T7' 'T9' 'CP5' 'FT9' 'PO3' 'Oz' 'POz' 'PO4'}
cdat=reshape(d.dwg,nsub*ncond,nchan,nchan);
cmat=[ones(1,nsub) ones(1,nsub)*2 ones(1,nsub)*3]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cdat=shiftdim(cdat,1);