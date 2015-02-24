% make sample fMRI data with various conditions
% simmat - a 3D matrix specifying which channels go together

function [dat,cdat]=makefdat(ncond,s,rois,dp,humps,simmat,mod);
close all;
dat=rand(s,ncond,rois,dp)/100;
%dat=rand(s,rois,dp);

if ndims(mod)<2
    error('need jitter values for all conditions and channels');
end
%if ndims(simmat)<3 | size(simmat,1)~=ncond
%    error('simmat wrong size');
%end

if dp<1000;
    disp('dp too short');
end

load('canon');

% basic random structure with correlated channels
for n=1:size(simmat,1);
   humper=zeros(s,ncond,size(simmat,1),dp);    
       for h=1:humps;
        hrand=fix(rand*dp);
        for subs=1:s;
            for chans=1:size(simmat,2);
                for c=1:ncond
                    hmod=fix(rand*mod(c,n));
                    humper(subs,c,chans,hrand+hmod:hrand+hmod+256)=canon;
                end
            end
        end
    end
    dat(:,:,simmat(n,:),:)=dat(:,:,simmat(n,:),:)+humper(:,:,:,1:dp);
end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%MAKE DATA TOTALLY RANDOM%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dat=rand(20,2,9,2000);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n=1:ncond;          %visualise
    dat1(n,:,:,:)=squeeze(dat(:,n,:,:));
end

viz4(dat1);colormap gray

dat=reshape(dat,ncond*s,rois,dp);



for subs=1:s*ncond;
    cdat(:,:,subs)=squareform1(pdist1(squeeze(dat(subs,:,:))));
end

%figure;viz(shiftdim(cdat,2));colormap jet

