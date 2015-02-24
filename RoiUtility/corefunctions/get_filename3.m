function pp = get_filename3(wcard)
% pp = get_filename3(wcard)
% wcard is like ls
%
% lists files; returns string matrix, one file per string
% operates in OSX.
% Returns wrong order!

p = ls(wcard);
wh = find(p == sprintf('\t') | p == sprintf('\n'));
x = diff(wh);
x = find(x > 1);
start = [0 wh(x)] + 1;
ends = [wh(x) length(p)] - 1;

pp = []; for i = 1:length(start), pp = str2mat(pp,p(start(i):ends(i)));, end
pp = pp(2:end,:);
