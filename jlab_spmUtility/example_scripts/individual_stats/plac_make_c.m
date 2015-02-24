function [mout1,mout2] = plac_make_c(msub)

	instr = repmat('%f ',1,16);

	try
		[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p] = ...
		textread(['behavior/' msub '.txt'],instr);
	catch
		msub
		error('Problem reading text file')
		
	end

 	mout1 = {a b c d e f g h};
	mout2 = {i j k l m n o p};

return
