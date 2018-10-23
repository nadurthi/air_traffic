function [T,W]=tens_prod_vec(u,v,wu,wv)
% 1 enitity is one row of any matrix u or v
% the rows of u are tensors producted with rows of v

if isempty(u)
    T=v;
    W=wv;
    return
end
if isempty(v)
    T=u;
    W=wu;
    return
end

n=size(u,1);
m=size(v,1);
T=[];
W=[];
for i=1:1:n
    T=vertcat(T,horzcat(repmat(u(i,:),m,1),v));
    W=vertcat(W,horzcat(repmat(wu(i),m,1),wv));
end
W=prod(W,2);
% W=W/sum(W);
% uc=size(u,2);
% vc=size(v,2);
% T=zeros(n*m,size(u,2)+size(v,2));
% k=1;
% for i=1:1:n
%     T(k,1:uc)=u(i,:);
%     for j=1:1:m
%         T(k,uc+1:uc+vc)=v(j,:);
%     end
%     k=k+1;
% end