% sorts a bunch of strings based on the numeric value of a field

function sorted_str = str_numsort(str, field_nums, field_separator)
    str = implode(cellstr(str), ' ');
    
    field_options = '';
    for i=1:length(field_nums)
        field_options = [field_options ' -k ' field_nums(i)];
    end
    command = sprintf('echo ''%s'' | sort -t ''%s'' -n %s', str, field_separator, field_options);
    disp(command);
    [status, sorted_str] = system(command);
end