function make_img(name,data,hdr)
% function make_img(name,data,hdr)
%
% name: basename of image file to write
% data: 3-D array of data to write as file.
% hdr: hdr struct, as read in with read_hdr.
%

fid= fopen([name '.img'],'wb');
switch hdr.datatype     
   	case 0
      fmt = 'int8';
   	case 2
      fmt = 'uint8';
   	case 4
      fmt = 'short';
   	case 8
      fmt = 'int';
   	case 16
      fmt = 'float';
   	case 32
      fmt = 'float';
      xdim = hdr.xdim * 2;
      ydim = hdr.ydim * 2;

   	otherwise
         error(['Data Type ' num2str(hdr.datatype) 'Unsupported. Aborting.']);
	end
fwrite(fid, data, fmt);
fclose(fid);

write_hdr([name '.hdr'],hdr);
return


