function [xint,wint,pw] = GenerateQuadPoints(m,ND,Np)

%
% Author: Puneet Singla.
% Last Modified: June 19, 2009.
%
% This function computes the orthogonal polynomials and their inner
% products in n-dim parameter space.
%
% Input Variables:
% N is the maximum order of polynomials along each direction and m detrmines the pdf.
% m = -2 corresponds to Gaussin pdf of Zero mean and variance 1
% m = -1 corresponds to uniform pdf over [-1,1]
% m >=0 corresponds to GLOMAP pdfs
% ND is number of uncertain parameters, i.e., our basis function will lie in
% ND dim space.
% Np number of integration points along each direction
%
% Output Variables:
% 
% xint is Ninteg x ND matrix of quadrature points
% wint,pw are Ninteg x 1 vectors of quadrature weights and pdf evaluations, respectively.
%


for ct = 1: ND
    if m(ct) == -2
        [xintg,wintg]=hermitequad(Np(ct),0); xintg = xintg(:); wintg = wintg(:);
        xind(:,ct) = xintg; wind(:,ct) = wintg; 
        weigh(:,ct) = 1/(2*pi)*exp(-0.5*xintg.^2);
        
        
    else if m(ct) == -1

            %[xintu,wintu]=lgwt(Np(ct),-1,1); xintu = xintu(:); wintu = wintu(:);
            [xintu,wintu]=ClenshawCurtis(Np(ct),-1,1); xintu = xintu(:); wintu = wintu(:);
            xind(:,ct) = xintu; wind(:,ct) = 1/2*wintu;
            weigh(:,ct) = 1/2*ones(size(xintu));
            
            
        else
             lb = -1; ub = 1;
            [xintgl,wintgl]=ClenshawCurtis(Np(ct),lb,ub);

            weig = GLOMAP(m(ct),xintgl); 
            xind(:,ct) = xintgl; wind(:,ct) = wintgl.*weig(:);
            weigh(:,ct) = weig(:);
        end
    end
end



%%
% ND Integration points and corresponding weights
%%

index = GenerateIndex(ND,Np);
xint = xind(index(:,1),1); % ND quadrature points
pw = weigh(index(:,1),1); % ND quadrature weights
wint = wind(index(:,1),1); % ND pdf

for j = 2:ND %good loop! - run through the number of dimensions!
    xint = [xint, xind(index(:,j),j)];
    wint = [wint, wind(index(:,j),j)];
    pw = pw.*weigh(index(:,j),j);

end


function weig = GLOMAP(m,x)
sumw = 0;
for k = 0:m
    facm = factorial(m); facmk = factorial(m-k);  fack = factorial(k);
    sumw = sumw+(-1)^k/(2*m-k+1)*facm/(fack*facmk)*abs(x).^(m-k);
end
weig = 1-abs(x).^(m+1).*sumw*factorial(2*m+1)/factorial(m)^2*(-1)^m;
    