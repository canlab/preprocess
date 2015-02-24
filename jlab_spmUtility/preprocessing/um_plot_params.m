function out = um_plot_params(dirname,varargin)
% out = um_plot_params(dirname,varargin)
%
% show realignment parameters
% start in main task directory for individual subject
% images should be in scan* subdirectories
%
% outputs:  4 numbers
% 1) max range in displacement rpy along any single dimension
% 2) max range xyz along any single dimension
% 3) max range in rpy (roll, pitch, yaw) within a session
% 4) max range in xyz within a session
%
% varargin is string for subject code, to append to name
% varargin 2 is to use the first n runs (enter an integer)
%
% example: 
% to read from CD's in Luis Hernandez output format, 
% try this:
%
% indx = 1;
%d = dir('d:'); fcode = lower(d(1).name); d = ['d:\' fcode '\func\intext\run_*\ra_img'];
%out = um_plot_params(d,fcode,8)        % for the 1st 8 runs
%scode{indx} = fcode; params(indx,:) = out; indx = indx+1;
%
% boxplot(params,1)
%scode = {};indx = 1;
%for i = 1:length(d)
%	fcode = d(i).name;
%	dd = [fcode '/scan*'];
%	out = um_plot_params(dd,fcode,8) 
%	if ~isempty(out)
%		scode{indx} = fcode; params(indx,:) = out; indx = indx+1;
%	end
%end


%if ~(exist('movement') == 7), mkdir movement, end
%try, cd movement,catch,end

% old

%!grep PARAMS ../*/realign.output > real_params.txt
%[run,file,c,d,e,f,g,h] = textread('real_params.txt','%s%s%f%f%f%f%f%f');
%params = [c d e f g h];

% new
myp = [];

file = get_filename2(fullfile(dirname,'realign.output'));
file = sortruns(file);
if length(varargin) > 1, n = varargin{2};, else, n = size(file,1);, end

for i = 1:min(n,size(file,1))
    if isempty(myp), 
        [myp,indx] = read_file(deblank(file(i,:)));
        run = ones(1,indx) .* i;
    else
        [tmp,indx] = read_file(deblank(file(i,:)));
        myp = [myp; tmp];
        run = [run ones(1,indx) .* i];
    end
end

params = myp;


% max rpy, max xyz
p1 = params(:,1:3);
p2 = params(:,4:6);

out(1) = max(max(p1) - min(p1));    % max range rpy along any single dimension
out(2) = max(max(p2) - min(p2));    % max range xyz along any single dimension

for i = 1:max(run), 
    tmp2(i) = max(max(p1(run==i,:)) - min(p1(run==i,:))); 
    tmp3(i) = max(max(p2(run==i,:)) - min(p2(run==i,:))); 
end
out(3) = max(tmp2);
out(4) = max(tmp3);

if length(varargin) > 0
    eval(['save ' varargin{1} 'real_params run file params out'])
else
    save real_params run file params out
end

figure; subplot 211; plot(params(:,1:3))
ylabel('movement in degrees')
legend({'pitch' 'roll' 'yaw'})

set(gcf,'Color','w')

subplot 212; plot(params(:,4:6))
legend({'x' 'y' 'z'})
ylabel('mm of displacement')


subplot 211 

% make title with subject path
% ---------------------------------------------------------
if length(varargin) > 0
    mytitle = varargin{1};
else
    mytitle = which('real_params.txt');
end
title(mytitle,'FontSize',14)


% plot lines and text to mark start of runs
% ---------------------------------------------------------

plot_lines(run)
subplot 212
plot_lines(run)


% set figure position (and print, if not commented out)
% ---------------------------------------------------------
set(gcf,'Position',[59   101   872   543])
try,print -dps2 -Pjjon-print1,catch,end

if length(varargin) > 0
    saveas(gcf,[varargin{1} 'real_params'],'fig')
    saveas(gcf,[varargin{1} 'real_params'],'tif')
else
    saveas(gcf,'real_params','fig')
    saveas(gcf,'real_params','tif')
end

%cd ..

return


% ---------------------------------------------------------
% ---------------------------------------------------------



% subfunctions
% ---------------------------------------------------------
function plot_lines(run)
hold on
yLim = get(gca,'YLim');

i = find(diff(run))';
plot([i+.5 i+.5]',repmat(yLim,length(i),1)','k')
		
for j = 1:length(i)
        st = i(j) - sum(run==j) + 5;
		text(st+5,yLim(2)-.2*yLim(2),['run ' num2str(j)])
        text(st+5,yLim(1)+.3,[num2str(sum(run == j)) ' imgs'])
        
end


function [ps,indx] = read_file(f)

fid = fopen(f);

tmp = '       '; indx = 1;

while tmp ~= -1
    
    tmp = fgetl(fid);
    if isempty(tmp) 
        tmp = 'aaaaaaaa';,
        
    elseif tmp ~= -1
    
        if strcmp(tmp(1:6),'PARAMS')
            [a] = sscanf(tmp,'%s%s%f%f%f%f%f%f');
            ps(indx,:) = a(end-5:end)';
            indx = indx + 1;
        end
        
    end
    
end
indx = indx - 1;
%fprintf(1,'Read %3.0f lines\n',indx)

return


    

function file = sortruns(file)

if size(file,1) > 10
    nplus = size(file,1) - 9;
    tmp = file(2:1+nplus,:);
    file(2:1+nplus,:) = [];
    file = [file; tmp];
end

return

