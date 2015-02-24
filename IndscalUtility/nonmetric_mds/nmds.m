function [ Config,DHS,DS,DeltaS,Stress1,StressT1,Stress2,StressT2,Rs,RsT] = ...
              nmds( userinput,TrainMatrix,TestMatrix,NTests )

global l Config no nt ns r DHS DS DO sortindex DeltaS referindex;
global effected effecteds oldsum1 oldsum2 noteffecteds;
global minoption;

r           = userinput{2};
nt          = userinput{3};
maxiter     = userinput{4};
convcrit    = userinput{5};
doprint     = userinput{6};
dograph     = userinput{7};
dospecific  = userinput{8};
startoption = userinput{9};
plot3dok    = userinput{10};
seed        = userinput{11};
dosym       = userinput{13};
dotest      = userinput{14};
minoption   = userinput{15};

DistsSub    = TrainMatrix; 

if (nt<1) 
   error( 'Number of dimensions must be bigger than 1' );
end

no = size( DistsSub , 1 ); % no is the Number of Objects 

if (no < 2)
   error('Loaded a similarity matrix with less than two objects' );
end

if dotest~=1
   StressT1 = 1.0;
   StressT2 = 1.0;
   RsT      = 0.0;
end

% do we symmetrize the matrix?
if dosym==1
   DistsSub = symmetrize( DistsSub );   
end

% manipulate the subjects' rating so that the triangular inequality and
% positivity holds

c1 = 0;
c2 = 0;
if (startoption==0)
   % make a random start configuration
   % changed from randn to rand [DJ]:
   rand( 'state',seed ); 
   Config = rand( no , nt ); % draw from normal distribution [NO LONGER NORMAL];
   
else % (startoption==1)
   % make an initial configuration with Torgeson-Young scaling
   %old code [DistsSub,c1,c2]    = manipulate( no , DistsSub  );
   [DistsSub,c1,c2] = ManipulateC( DistsSub,no  );
   Config = calcratconfig( nt , no , DistsSub );
end

% convert these two dimensional arrays to single dimension arrays
ns = no * (no - 1) / 2; % ns is #entries in single dimension array

Delta      = zeros( ns , 1 ); % the subjects' data 
DO         = zeros( ns , 1 ); % the distances of the start configuration
referindex = zeros( no , no );

ind = 0; % ind is the index
for i=1:no-1
   for j=i+1:no
      ind = ind + 1;
      Delta( ind )      = DistsSub( i,j );  
      DO( ind )         = Mink( Config( i,: ) , Config( j,: ) , r );
      
      referindex( i,j ) = ind;
   end
end

% figure out what indices are effected when object l is referred to
% this will speed up several processes
effected  = zeros( no , no-1 );
effecteds = zeros( no , no-1 );

for l=1:no
  count = 0; 
  if (l<no)
    i=l;
    for j=i+1:no
      ind = referindex( i,j );
      count = count + 1;
      effected( l , count ) = ind; 
    end
  end  

  if (l>1)
    j=l;
    for i=1:j-1;
      ind = referindex( i,j );
      
      count = count + 1;
      effected( l , count ) = ind;
    end
  end   
end

% sort the subjects' dissimilarity array
[DeltaS,sortindex]=sort( Delta );

% update the effected array to incorporate this order
[dummy,newindex] = sort( sortindex );
wholelist = 1:ns;
noteffecteds = zeros( no , ns - no + 1 );
for l=1:no
   for i=1:no-1
      effecteds( l,i ) = newindex( effected( l,i ));
   end
   
   notinlist = setdiff( wholelist , effecteds( l,:) );
   noteffecteds( l,: ) = notinlist;
end

% and sort the distances accordingly
DS = DO(sortindex);

% calculate the d-hats from the distances with Kruskal's upward-downward algorithm
[DHS,inblock,ave,rsum,nmem] = calcdhatsc( ns , DS );

% calculate stress
[Stress1,Stress2] = calcstress( DS , DHS );

% calculate the spearman rankorder correlation coefficient
Rs = compute_rank( ns , DS , DeltaS );

if dograph==1
   showplots( 0 , Stress1 , Stress2 , plot3dok );
end;   
   
% do the optimization of the coordinates

if doprint==1
   % print some statements about the initial status of the MDS program
   fprintf( 'Values added to satisfy distance axioms %g %g\n' , c1 , c2 );
   fprintf( 'Minkowski metric r value                %g\n'   , r );
   fprintf( 'Number of objects in matrix             %d\n'   , no );
   fprintf( 'Number of useful similarity ratings     %d\n'   , ns );
   fprintf( 'Number of dimensions:                   %d\n\n' , nt ); 
   fprintf( 'Maximum number of iterations            %d\n' , maxiter  );
   fprintf( 'Convergence Criterion:                  %g\n' , convcrit );
   fprintf( 'Minimize function:                      stress%d\n' , minoption );
   fprintf( 'Starting configuration is:              ' );
   
   if (startoption==0)
      fprintf( 'random  (seed=%d)\n' , seed );
   elseif (startoption==1)
      fprintf( 'Torgeson-Young' );
   else 
      fprintf( 'User-supplied' );
   end
   
   fprintf( '\n\n' );
   fprintf( 'The following results apply to training file only:\n' );
   fprintf( '(Rs is the rank order correlation coefficient\nbetween observed and predicted dissimilarities\n\n' );
   fprintf( 't=0 Stress1=%1.5f Stress2=%1.5f Rs=%1.5f\n' , Stress1 , Stress2 , Rs );
end

iter     = 0;
stopiter = 0;

if (iter==maxiter) 
   stopiter=1;
end

while (stopiter==0);
   iter=iter+1;
   
   if (iter==maxiter)
      stopiter=1;
   end
   
   % on each iteration, move all points around
   optimize( dospecific );
   
   % calculate the new distances for the updated configuration
   ind = 0; % ind is the index
   for i=1:no-1
      for j=i+1:no
         ind = ind + 1;
         DO( ind ) = Mink( Config( i,: ) , Config( j,: ) , r );
      end
   end
 
   % and sort the distances according to subjects sorted array
   DS = DO(sortindex);
   
   % calculate the new d-hats from the distances with Kruskal's upward-downward algorithm
   [DHS,inblock,ave,rsum,nmem] = calcdhatsc( ns , DS );
   
   % calculate new stress
   if iter>1
      if minoption==1
         LastStress = Stress1;
      else
         LastStress = Stress2;
      end
   end;   
   
   [Stress1,Stress2] = calcstress( DS,DHS );
   
   if isnan( Stress1 )
      error( 'A NaN occured in calculations (1)' );
   else
     % What is the change with respect to last Stress Value
     if iter>1
        if minoption==1
           Change = ( LastStress - Stress1 );
        else
           Change = ( LastStress - Stress2 );
        end
        
        if Change < convcrit 
           stopiter = 2;
        end
     end;
        
     % calculate the spearman rankorder correlation coefficient
     Rs = compute_rank( ns , DS , DeltaS );
          
     if (doprint==1)
       fprintf( 't=%d Stress1=%1.5f Stress2=%1.5f Rs=%1.5f\n' , iter , Stress1 , Stress2 , Rs );
     end
   
     if dograph==1
       showplots( iter , Stress1 , Stress2 , plot3dok );  
     end;
   end 
end

if (dotest==1) & (stopiter~=3)
   
   if (doprint==1)
     fprintf( '\n\n' );
   end
   
   no_t = size( TestMatrix , 2 );
      
   if (no_t ~= no) 
      error( '# objects in validation matrix unequal to # objects in calibration matrx' );
   end   
      
   StressT1 = zeros( 1,NTests );
   StressT2 = zeros( 1,NTests );
   RsT      = zeros( 1,NTests );
   
   for s=1:NTests
      
      if (NTests==1)
         DistsSubT = squeeze( TestMatrix( 1  , : , : ));
      else
         DistsSubT = squeeze( TestMatrix( s  , : , : ));
      end
      
      if dosym==1
         DistsSubT = symmetrize( DistsSubT );
      end
      
      DeltaT = zeros( ns , 1 );
      
      ind = 0; % ind is the index
      for i=1:no-1
         for j=i+1:no
            ind = ind + 1;
            
            DeltaT( ind ) = DistsSubT( i,j );
         end
      end
      
      [DeltaST,sortindexT]=sort( DeltaT );
      
      % calculate the d-hats for the test data
      DS_T = DO(sortindexT);
      [DHS_T,inblock,ave,rsum,nmem] = calcdhatsc( ns , DS_T );
      
      % calculate new stress
      [ StressT1(s) , StressT2(s) ] = calcstress( DS_T,DHS_T );
      
      if isnan( StressT1(s) )
         error( 'A Nan occurred here' );
      end
            
      % calculate the spearman rankorder correlation coefficient
      RsT(s) = compute_rank( ns , DS_T , DeltaST );
      
      if (doprint==1)
         fprintf( 'S=%d Stress1=%1.5f Stress2=%1.5f Rs=%1.5f\n' , s , StressT1(s) , StressT2(s) , RsT(s) );
      end   
   end
end

if doprint==1
   fprintf( '\nEnd of simulation - reason: ' );
   if (stopiter==1) 
      fprintf( 'maximum number of iterations reached\n\n' );
   elseif (stopiter==2)
      fprintf( 'convergence criterion reached\n\n' );
   elseif (stopiter==3)
      fprintf( 'a NaN occured in stress computation\n\n' );
   
   end   
end

clear DistSub
%-----------------end of mds function----------------------------------



%----------------------------------------------------------------------
function  optimize( dospecific );
% optimize Stress for each point separately but jointly over each dimension

global l Config no nt ns r DHS DO DS sortindex oldsum1 oldsum2;
global noteffecteds;

options(1)=-1;

for l=1:no
   % optimize position of point l
   
   [oldsum1,oldsum2 ] = NotEffectedC( noteffecteds , DS , DHS , l , ns , no );
      
   x0 = Config(l,:); 
   warning off;
   % using the matlab program sopt2 that calls a C function
   x  = fminu1( 'sopt2',x0,options  );
   warning on;
   Config(l,:)=x;
  
   % center all dimensions again
   for e=1:nt
     Config(:,e) = Config(:,e) - mean( Config(:,e) );
   end
   
   % and rescale
   Config = Config / std( Config(:) );
end



%--------------------------------------------------------------------
function showplots( iter , Stress1 , Stress2 , plot3dok );

global Config DeltaS DHS DS no nt;

% plot the d's and d-hats vs. the data
figure( 1 );
subplot( 1,2,1 );

if (nt==1)
   plot( Config(:,1) , Config(:,1) , 'r*' );
   grid on;
   axis square;
   xlabel( 'dimension 1' );
   ylabel( 'dimension 1' );
   
   for i=1:no
      text( Config(i,1)+0.1 , Config(i,1) , [sprintf( '%d' , i )] );
   end;   
elseif (nt==2)
   plot( Config(:,1) , Config(:,2) , 'r*' );
   grid on;
   axis square;
   xlabel( 'dimension 1' );
   ylabel( 'dimension 2' );
   
   for i=1:no
      text( Config(i,1)+0.1 , Config(i,2) , [sprintf( '%d' , i )] );
   end;   
elseif (nt>2) & (plot3dok==1)
   plot3( Config(:,1) , Config(:,2) , Config(:,3) , 'r*' );
   grid on;
   axis square;
   hold on;
   
   for i=1:no
      plot3( [Config(i,1) Config(i,1)] , [Config(i,2) Config(i,2)] , [ 0 Config(i,3)] , 'k' );
   end
   
   hold off;
   
   xlabel( 'dimension 1' );
   ylabel( 'dimension 2' );
   zlabel( 'dimension 3' );
   
   for i=1:no
      text( Config(i,1)+0.1 , Config(i,2) , Config(i,3) , [sprintf( '%d' , i )] );
   end;  
end

title( sprintf( 'Configuration after %d iterations' , iter ));

subplot( 1,2,2 );
plot( DS , DeltaS , 'r*' , DHS , DeltaS , '-gs' );
grid on;
axis square;
axis tight;
ylabel( 'dissimilarities data' );
xlabel( 'd' );
title( [ sprintf( 'Stress1=%1.4f Stress2=%1.4f ' , Stress1 , Stress2 ) ] );

drawnow;


%--------------------------------------------------------------------
function Config = calcratconfig( nt , no , DistsSub );
I = eye( no );
U = ones( no );
Z = I - (1/no) * U;
D = DistsSub.^2;
% convert the distances to inner products assuming centroid origin
B = -0.5 * Z * D * Z;

% do a eigenvector decomposition of B ( B = X * X' )
[V,D]=eig(B);
DD = diag( D );
[DDS,index]=sort(-DD);
DDS=-DDS;

DS=D(:,index);
VS=V(:,index);

% if we want n dimensions, there should be n positive eigenvalues
if DS(nt) < 0
   warning( 'negative eigenvalue in matrix: decrease number of dimensions' );
end

DS = DS(:,1:nt);
DS = DS.^0.5;

Config=VS*DS;

% normalize the variance of this starting configuration
Config = Config / std( Config(:) );


%--------------------------------------------------------------------
function [KSS1,KSS2,SStress] = calcstress( DOS , DHS );
% Kruskal's Stress formula 1
S1a  = sum( (DOS - DHS).^2 );
T1a  = sum( DOS.^2 );
KSS1 = sqrt( S1a / T1a ); 

% Kruskal's Stress formula 2
avDOS = mean( DOS );
T1a   = sum( ( DOS - avDOS ) .^ 2 );
KSS2  = sqrt( S1a / T1a );

% SStress
S1a     = sum( (DOS.^2 - DHS.^2).^2 );
T1a     = sum( DOS.^4 );
SStress = sqrt( S1a / T1a ); 

%--------------------------------------------------------------------
function d=Mink( x , y , r )
d1 = abs(x - y);
d = sum( d1.^r ) ^ (1/r);

%--------------------------------------------------------------------
function D = symmetrize( D );
global no;

for i=1:no-1
   for j=i+1:no
      rating1 = D( i,j );
      rating2 = D( j,i );
         
      newrating = (rating1 + rating2) / 2;
         
      D(i,j) = newrating;
      D(j,i) = newrating;
   end
end
   
for i=1:no
  D(i,i) = 0;
end





