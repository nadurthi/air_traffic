function [phi,phip,InnerProd,weigh,xind,wind] = GramSch(n,m,Np)

%
% This function implements the GramSchmidt procedure to generate orthogonal
% polynomials to a given weight function and compute corresponding inner
% product
%
% Input Variables:
% n is the desired order of polynomials
% m is the order of weight function (pdf), m=-2 coreesponds to std.
% Gaussian, m=-1 corresponds to uniform and m > -1 corresponds to GLOMAP
% weight functions
% Np is number of integration points
%
% Output Variables:
%  phip is n x N matrix of polynomial values at prescribed points, xp
%  InnerProd is n x n x n tensor of basis function inner product.
% weigh is N x 1 vector of pdf evaluations.


syms x real


%%
% Weight function Expression
%%
if m == -2
    weig = 1/(2*pi)*exp(-0.5*x^2); 
    lb = -inf; ub = inf;
    [xind,wind]=hermitequad(Np,0);
else if m == -1
        weig = 1/2;
        [xind,wind]=lgwt(Np,-1,1); lb = -1; ub = 1;
    else
        weig = GLOMAP(m,x); lb = -1; ub = 1;
    end
end

%%
% Compute Orthogonal Polynomials
%%
phi(2,1) = x; phi(1,1)=1;

for ct = 1:n+1
    psi(ct) = x.^(ct-1);
end

for ct = 3:n+1
    sump = 0;
    for k = 1: ct-1
        sump = sump + int(psi(ct)*phi(k)*weig,'x',lb,ub)/int(phi(k)*phi(k)*weig,'x',lb,ub)*phi(k);
    end
    phi(ct,1) = psi(ct)-sump;
end

% F  = [x,phi];
% matlabFunction(F,'file','BasisFunc.m')

for ct = 1:n+1
    InnerProd(:,:,ct) = (double(int(phi*phi'*phi(ct,1)*weig,'x',lb,ub)));
end

xp = xind;


for ct = 1:length(xp)
    x = xp(ct);
    phip(:,ct) = eval(phi);
    if m==-1
        weigh(ct,1) = 1/2;
    else
    weigh(ct,1) = eval(weig);
    end
end


function weig = GLOMAP(m,x)
sumw = 0;
for k = 0:m
    facm = factorial(m); facmk = factorial(m-k);  fack = factorial(k);
    sumw = sumw+(-1)^k/(2*m-k+1)*facm/(fack*facmk)*abs(x).^(m-k);
end
weig = 1-abs(x).^(m+1)*sumw*factorial(2*m+1)/factorial(m)^2*(-1)^m;