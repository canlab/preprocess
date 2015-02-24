function  write_img_data(name, data, hdr)
% Luis hernandez
% last edit 4-6-98
%
% function  write_img_data(name, data, hdr)
%
% Writes the data to an analyze format file 'name' containing mutislice image data 
% Use basenames to leave off extensions
% Also writes header
% Modified 3/08/01 by Tor Wager

basename = name;
name = [name '.img'];
disp(['output name is ' name])

   [pFile,messg] = fopen(name, 'wb');
   if pFile == -1
      errormesg(messg);   
      return;
   end
   
      
   switch hdr.datatype     
   case 2
      fmt = 'char';
   case 4
      fmt = 'int16'; %'short';
   case 8
      fmt = 'int';
   case 16
      fmt = 'float';
   case 32
      fmt = 'float';
           
   otherwise
      errormesg(sprintf('Data Type %d Unsupported. Aborting',hdr.datatype));
      return
      
   end

   
   fwrite(pFile, data, fmt);
   fclose(pFile);
   
   disp('Writing header')
   name = [basename '.hdr'];
   write_hdr(name,hdr);
   
 return

