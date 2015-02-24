%name_networks_plugin

donames = dointeractive;
if ~isfield(DATA,'APPLY_CLUSTER')
    DATA.APPLY_CLUSTER = [];
elseif isfield(DATA.APPLY_CLUSTER,'names') && dointeractive
    donames = input('Names found in DATA.APPLY_CLUSTER.names.  Re-name (1/0)? ');
end


if donames


    classes = DATA.CLUSTER.classes;

    % name 'networks'
    for i = 1:max(classes)

        disp(['Network ' num2str(i)])
        tmp = find(classes == i);
        for j = 1:length(tmp)
            fprintf(1,'%s . ',DATA.CLUSTER.names{tmp(j)});
            if j > 4, fprintf(1,'\n'),end
        end
        fprintf(1,'\n')
        DATA.APPLY_CLUSTER.names{i} = input('Name this network: ','s');

        if length(DATA.APPLY_CLUSTER.names{i}) < 2, DATA.APPLY_CLUSTER.names{i} = [DATA.APPLY_CLUSTER.names{i} ' '];,end
    end

    disp('Saved network names in DATA.APPLY_CLUSTER.names')


end     % end donames

