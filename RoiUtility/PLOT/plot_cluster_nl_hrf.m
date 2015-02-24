function plot_cluster_nl_hrf(myd,EXPT)
% function plot_cluster_nl_hrf(myd,EXPT)
%
% Tor Wager
%
% For Batch:
% D = dir; D=str2mat(D.name); wh=strmatch('snpm',D); E = D(wh,:);
% for i = 1:size(E,1),plot_cluster_nl_hrf(deblank(E(i,:)),EXPT);, end

% load file

eval(['cd ' EXPT.studydir filesep myd])
D = dir;
D = str2mat(D.name);
for i = 1:size(D,1), 
    a = findstr('clusters.mat',D(i,:));, 
    if ~isempty(a), mynm = deblank(D(i,:));, end
end
if ~exist('mynm') == 1, cd .., return, end
eval(['load ' mynm])

mys = findstr(myd,'000');
if isempty(mys), mys = findstr(myd,'00');,end
if isempty(mys), error('Dir name does not contain contrast number (xxx000#xxx or 00##).'),end 
mynum = str2num(myd(mys:mys+3));
wcon = zeros(1,size(EXPT.DX.contrasts,1));
wcon(mynum) = 1;

cd ..

% loop through clusters and plot each
for i = 1:length(clusters)

    nl_hrf_plot(clusters(i).XYZmm',EXPT,3,wcon);
    title(['Contrast Height: cl' num2str(i) ' ' myd])
    saveas(gcf,[myd '_cl_' num2str(i) '_nlcon_' num2str(mynum)],'jpg')
    close
    title(['Height: cl' num2str(i) ' ' myd])
    saveas(gcf,[myd '_cl_' num2str(i) '_nlfit_' num2str(mynum)],'jpg')
    close
end


return
