function d=unifrom_moments_iid(dim,N)
%% m is the vector of moments
% N is the moment(always keep it even)
% n is dimension of system

mom=@(n)((-1-1)^n+(1+1)^n)/(2^(n+1)*(n+1));

n=dim;
% combos = combntns(1:N,2);
% ind=combntns(1:length(combos),N/2);
% A=[];
% a=[1:1:N];
% for i=1:1:length(ind)
%     r=[];
%     for j=1:1:N/2
%         r=horzcat(r,combos(ind(i,j),:));
%     end
%     if sum(abs(sort(r)-a))==0
%     A=vertcat(A,r) ;
%     end   
% end

% combos = combntns(0:N,n)
combos = GenerateIndex(n,(N+1)*ones(1,n));
combos(find(combos==(N+1)))=0;

x=[];
for i=1:1:length(combos)
    if sum(combos(i,:))==N
     x=vertcat(x,combos(i,:));
%      x=vertcat(x,wrev(combos(i,:)));
    end
end
size(x);
% x=vertcat(x,N/2*ones(1,n));
x=sortrows(x,-1);
nn=size(x);
% m=zeros(nn(1),1);
g=ones(nn(1),1);
h=[];
p=[];
for i=1:1:nn(1)
    for j=1:1:nn(2)
    g(i)=g(i)*mom(x(i,j));
    end
if g(i)~=0
   p=vertcat(p,g(i));
   h=vertcat(h,x(i,:));
end
end
d=[h,p];
end