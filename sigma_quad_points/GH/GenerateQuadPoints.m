function [xint,wint] = GenerateQuadPoints(Pcov,mcent,Np)

%
% Author: Puneet Singla.
% Last Modified: October 1st, 2010.
%
% This function computes the orthogonal polynomials and their inner
% products in n-dim parameter space.
%
% Input Variables:
% P is the covariance matrix
% mcent is the center of the Gaussian component
% Np number of integration points along each direction
%
% Output Variables:
% 
% xint is Ninteg x ND matrix of quadrature points
% wint Ninteg x 1 vector of quadrature weights.
%


[U,S,V] = svd(Pcov); S = diag(S);

if det(U-V) >= 1e-12
    disp('Covariance Matrix should be Symmteric');
    return
end

ND  =size(Pcov,1); % state dimension

for ct = 1: ND

        [xintg,wintg]=HermiteQuad(Np,0,S(ct)); xintg = xintg(:); wintg = wintg(:);
        xind(:,ct) = xintg; wind(:,ct) = wintg;
        
    
end


 

%%
% ND Integration points and corresponding weights
%%

index = GenerateIndex(ND,Np*ones(1,ND));
xint = xind(index(:,1),1); % ND quadrature points
wint = wind(index(:,1),1); % ND pdf

for j = 2:ND %good loop! - run through the number of dimensions!
    xint = [xint, xind(index(:,j),j)];
    wint = wint.*wind(index(:,j),j);


end
xint = xint'; xint = U*xint+repmat(mcent,1,size(xint,2));xint = xint';


    