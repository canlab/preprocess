function [line_nbr, byte_offset] = locate_line(filename, str)
% locate_line - returns the line nbr and byte offset of the first line containing the
% string

fid = fopen(filename);
if(fid == -1)
    error('Could not open data file %s.', filename);
    return;
end
tline = '';
line_nbr = 0;
try
    while(1)
        tline = fgetl(fid);
        if(~ischar(tline))
            fclose(fid);
            disp(sprintf('Did not find raw data in eyetracking data file %s.', filename));
            return;
        end
        line_nbr = line_nbr + 1;
        if(~isempty(strfind(tline, str)))
            byte_offset = ftell(fid);
            break;
        end
    end
    fclose(fid);
catch
    fclose(fid);
end

end