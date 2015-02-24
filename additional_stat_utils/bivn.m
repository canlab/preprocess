function [Z]=bivn(X,Y,muX,muY,sigmaX,sigmaY,rho)

% 
% Usage:
% 
% [Z]=bivn(X,Y,muX,muY,sigmaX,sigmaY,rho)
% 
% Returns a matrix of equal size to matrices X and Y (which must be of
% identical size) of densities under the bivariate normal distribution with
% parameters definated by the remaining inputs.


Z=zeros(size(X,1),size(X,2));

for k=1:size(X(:),1)
    Z(k)=(1/(2*pi*sigmaX*sigmaY*sqrt(1-rho^2)))*...
        exp(-1/(2*(1-rho^2))*(((X(k)-muX)/sigmaX)^2-2*rho*((X(k)-muX)/sigmaX)*((Y(k)-muY)/sigmaY)+((Y(k)-muY)/sigmaY)^2));
end