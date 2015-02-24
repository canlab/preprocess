function html = canlab_html_creator(filename, filescript, title, subject)

%% initialize html code
html{1,1} = '<html>';
html{2,1} = '<head>';
html{3,1} = sprintf('<title>%s</title>',title);
html{4,1} = '</head>';
html{5,1} = '<body>';
html{6,1} = sprintf('<h2>%s</h2>', subject);

%% write images and descriptions
imagenum = size(filename,1);
for n = 1:imagenum
    html{7+(n-1)*2,1} = sprintf('<div><img src="%s"></div>',filename{n});  %image
    html{8+(n-1)*2,1} = sprintf('<div>%s</div>',filescript{n});              %related text
    html{7+(n-1)*2,1} = strrep(html{7+(n-1)*2,1},'\','\\')      %for PCs, make file writeable using fprintf
end
%% close html code
html{7+imagenum*2,1} = '</body>';
html{8+imagenum*2,1} = '</html>';

