img = '/data2/inhib/011018cg/anatomy/homocor_sn3d.mat';
meanadir = '/data2/inhib/011018cg/';


f99_sn3d(2,...																			% option
      img, ...																			% determine params from
      0,'',... 																			% m	asking
      [meanadir 'flanker/scan5';meanadir 'flanker/scan6'],... 							% directories
      ['ravol*img';'ravol*img'],... 													% EXP: file name wildcards
         [], ...																		% template
         0,[],... 																		% masking of canonical				
         'affine3',[],13,16,0.01,... 													% type and bounding box - neurological
         1,[],[],'') 																	% other params & voxel size - bilinear interp