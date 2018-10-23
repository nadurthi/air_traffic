function [Xnew,wnew]=Reweightpts(Nmoms,Xprior,wprior,hx,R,ymeas,Xpost,wpost)
% [y,M] give the moments of the posterior to be captured by weights
% pts position X stays the same.
% postfunX is posterior pdf evaluations at pts X
dim=size(Xpost,2);
N=length(wpost);


%% posterior moments
% y=[];
% M=[];
py=0;
ss=zeros(size(Xpost,1),1);
invR=inv(R);
y=zeros(1,size(Xpost,2));
M=1;
for i=1:1:N
    hh=hx(Xprior(i,:));
    ss(i)=-0.5*(ymeas(:)'-hh(:)')*invR*(ymeas(:)'-hh(:)')'; 
    py=py+wprior(i)*1/sqrtm(det(2*pi*R))*exp(-0.5*(ymeas(:)'-hh(:)')*invR*(ymeas(:)'-hh(:)')'); 
end

wmod=((wprior/py)*1/sqrt(det(2*pi*R))).*exp(ss);

n1=0;
n2=0;
for i=1:1:Nmoms
    [yy,MM]=Cal_moments_samples(Xprior,wmod,i,'raw');
    y=vertcat(y,yy);
    M=vertcat(M,MM);
    if i==1
        n1=length(MM);
    end
    if i==2
        n2=length(MM);
    end
    
    
end
[mu,P]=MeanCovPts(wmod,Xprior);
[Xnew,wnew]=conjugate_dir_gausspts_till_8moment(mu(:),P);

% [ypost,Mpost]=Cal_moments_samples(X,w,2,'raw');
% keyboard
% %%
% [mupost,Ppost]=MeanCovPts(wmod,Xprior);
% [muprior,Pprior]=MeanCovPts(wprior,Xprior);
% 
% isqPprior=sqrtm(inv(Pprior));
% sqPpost=sqrtm(Ppost);
% 
% 
% 
% % XX=zeros(size(Xprior));
% % for i=1:1:length(wprior)
% %     XX(i,:)=isqPprior*(Xprior(i,:)'-muprior);
% %     XX(i,:)=sqPpost*XX(i,:)'+mupost;
% % end
% % X=XX;
% % w=wprior;
% 
% X=Xpost;
% w=wpost;
% 
% A=zeros(size(y,1),N);
% B=zeros(size(y,1),1);
% 
% 
% for i=1:1:size(y,1)
% A(i,:)=prod(X.^repmat(y(i,:),N,1),2);
% B(i)=M(i);
% end
% 
% %% quad prog
% % Aeq=A(1:(1+n1+n2),:);
% % Beq=B(1:(1+n1+n2));
% % A(1:(1+n1+n2),:)=[];
% % B(1:(1+n1+n2))=[];
% % 
% % A=A(1:n1+n2,:);
% % B=B(1:n1+n2);
% 
% % opts = optimset('Algorithm','interior-point-convex','Display','on');
% % wnew=quadprog(A'*A,-A'*B,[],[],Aeq,Beq,zeros(1,N),ones(1,N));
% wnew=quadprog(eye(N),zeros(N,1),[],[],A,B,zeros(1,N),ones(1,N));



