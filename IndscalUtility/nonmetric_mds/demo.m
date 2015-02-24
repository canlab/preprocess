function Config=demo(A,ndims);
% ----------------------------------------------------------------------------------
%    MULTIDIMENSIONAL SCALING in matlab by Mark Steyvers 1999 
%      needs optimization toolbox      
% 
%
% ----------------------------------------------------------------------------------

no 			= size(A,1);						% size of input matrix
%A  			= rand( no , no );	% dissimilarities

plotdim1 	= 1;		% which dimension on x-axis
plotdim2 	= 2;		% which dimension on y-axis
R			= 2.0;	% R values in minkowski metric
maxiter		= 100
;		% maximum number of iterations
conv     	= 0.001;	% convergence criterion
R1       	= 1;     % 1=Torgeson Young scaling for initial configuration 0 = intial random configuration
seed     	= 1;		% seed for random number generator
minoption 	= 1;  	% 1 = minimize stress1    2 = minimize stress2

% provide some text labels for the stimulus points here
for i=1:no
   labels{ i } = sprintf( '%d' , i );
end

userinput{2} = R;
userinput{3} = ndims;
userinput{4} = maxiter;
userinput{5} = conv;
userinput{6} = 1;         % 0=no 1=yes, printed comments
userinput{7} = 0;
userinput{8} = 0;
userinput{9} = R1;
userinput{10}= 0;
userinput{11}= seed;
userinput{13}= 1;         % 0=do not symmetrize input matrix % 1=do symmetrize
userinput{14}= 0;
userinput{15}= minoption;




% ---------------------------------------------------------------------------------
%    CALL THE MDS ROUTINE
[ Config,DHS,DS,DeltaS,Stress1,StressT1,Stress2,StressT2,Rs,RsT] = ...
              nmds( userinput,A );
% ---------------------------------------------------------------------------------



