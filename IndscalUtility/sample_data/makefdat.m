% make sample fMRI data
% simmat - a 3D matrix specifying which channels go together
function [dat,cdat]=makefdat(s,rois,dp,humps,simmat,mod);
close all;
dat=rand(s,rois,dp)/100;
%dat=rand(s,rois,dp);

if dp<1000;
    disp('dp too short');
end

load('canon');

% basic random structure with correlated channels

for n=1:size(simmat,1);
    humper=zeros(s,size(simmat,1),dp);
    for h=1:humps;
        hrand=fix(rand*dp);
        for subs=1:s;
            for chans=1:size(simmat,2);
            hmod=fix(rand*mod);
            humper(subs,chans,hrand+hmod:hrand+hmod+256)=canon;
            end
        end
    end
    dat(:,simmat(n,:),:)=dat(:,simmat(n,:),:)+humper(:,:,1:dp);
end


%dat=rand(s,rois,dp);
viz(dat);colormap gray

for subs=1:s;
    cdat(:,:,subs)=squareform1(pdist1(squeeze(dat(subs,:,:))));
end

figure;viz(shiftdim(cdat,2));colormap jet

