function [xu,Pu]=EKF_disc(model,xu,Pu,ym)
%this is discrete dynamics

if model.propagate==1
    if strcmp(model.dynamics,'discrete')
        
        
    else
        
        [tt,xx]=ode45(@(t,x)lorenz_cont_ekf_prop(t,x,model),[0,model.dt],[xu(:);reshape(Pu,model.fn^2,1)]  );
        xf=xx(end,1:model.fn)';
        Pf=reshape(xx(end,model.fn+1:end),model.fn,model.fn);
        Pf=Pf+model.Q;
        
        %     [tt,xx]=ode45(@model.fx,[0,model.dt],xu(:));
        %     xf=xx(end,:)';
        %
        %     F=model.fx_jac_cont;
        %     X=[x;P]
        %     [tt,xx]=ode45(@(t,X)reshape( F(t,X(1:model.fn)) , model.fn^2,1),[0,model.dt],reshape(Pu,model.fn^2,1) );
        %     Pf=BB*Pu*BB'+model.Q
        
    end

else
    xf=xu;
    Pf=Pu;
    
    
end


%update
if ym==-1234
    xu=xf;
    Pu=Pf;
else
    
    AA=model.hx_jac(xf,model.para_dt);
    try
        K=Pf*AA'*inv(AA*Pf*AA'+model.R);
    catch
        keyboard
    end
    
    xu=xf+K*(ym-model.hx(xf,model.para_dt));
    Pu=(eye(model.fn)-K*AA)*Pf;
end
end

