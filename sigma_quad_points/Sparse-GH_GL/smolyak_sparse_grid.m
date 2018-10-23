function [x,w]=smolyak_sparse_grid(d,l,type)
% d is the dimension of the system
% type is eith 'GH' for gauss hermite
%              'GLgn'for gauss Legendre
% l is the number of points of the method type in 1 dimension.
% therefore the smolyak sparse grid scheme produces a 
% quadrature rule that can integrate all polynomials of total degree
% 2l-1.

switch lower(type)
    case 'gh'
        QUADD=@(m)GH_points(0,1,m);
    case 'glgn'
        QUADD=@(m)GLeg_pts(m, -1, 1);
    otherwise
        disp('Gods must be crazy; They made you!!!')
        return;
end

if d==1
    [x,w]=QUADD(l);
    return
end

x=[];
w=[];
for i=1:1:l
    [x1,w1]=QUADD(i);
    if i-1>0
    [x2,w2]=QUADD(i-1);
    x2=-x2;
    w2=-w2;
    [xd1,wd1]=smolyak_sparse_grid(d-1,l-i+1,type);
    [xd,wd]=tens_prod_vec(x1,xd1,w1,wd1);
    [xdm,wdm]=tens_prod_vec(x2,xd1,w2,wd1);
    x=vertcat(x,xd);
    w=vertcat(w,wd);
    x=vertcat(x,xdm);
    w=vertcat(w,wdm);

    else
     [xd1,wd1]=smolyak_sparse_grid(d-1,l-i+1,type);   
     [xd,wd]=tens_prod_vec(x1,xd1,w1,wd1);
     x=vertcat(x,xd);
     w=vertcat(w,wd);
    end

end

%% Now finding the duplicate points and adding their weights 
ep=1e-12;
i=1;
while i<=length(w)
%     if i>length(w)
%         break
%     end
    I=find(sum((x-repmat(x(i,:),size(x,1),1)).^2,2)<ep);
    if length(I)>1
    w(i)=sum(w(I));
    x(I(2:end),:)=[];
    w(I(2:end))=[];
    end
i=i+1;    
end

