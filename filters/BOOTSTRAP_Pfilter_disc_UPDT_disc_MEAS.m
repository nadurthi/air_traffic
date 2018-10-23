function [X_pf,w_pf]=BOOTSTRAP_Pfilter_disc_UPDT_disc_MEAS(model,pf,X_pf,w_pf,ym)
%Modified code of GABRIEL to suit my purpose

%% Propagate Particles
 X_pf = PF_time_update1(X_pf, model);
 w_pf=w_pf;
 %% CHECK MEasurement update if available
        if ym==-1234
% [PF.x{k},PF.P{k},PF.minB{k},PF.maxB{k},tmp_X_pf] = getPFdata(X_pf, w_pf);
         return
        end
        
 %% Meas UPDT
 [X_pf, w_pf] = PF_meas_update(X_pf, w_pf, model, pf, ym);
 
%  %% GET estimates
%  [PF.x{k},PF.P{k},PF.minB{k},PF.maxB{k},tmp_X_pf] = getPFdata(X_pf, w_pf);
end