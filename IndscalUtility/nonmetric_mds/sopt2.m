function f = sopt2( x );
% calculate the stress based on current position
global l Config no nt ns r DHS DO DS referindex;
global effected effecteds oldsum1 oldsum2 minoption;

% call the C FUNCTION
f = sopt2C( Config , l , x , no , nt , ns , referindex , r , DO , ...
   				DS , DHS , effected , effecteds , oldsum1 , oldsum2 , minoption );





