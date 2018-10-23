function [mk1,Pk1]=KalmanUpdate(mk,Pk,mz,Pz,Pcc,ym)
             K=Pcc/Pz;
%              keyboard
             if max(isnan(ym))==1 || max(isempty(ym))==1
                 mk1=NaN;
                 error('no measurement here')
             else
             mk1=mk+K*(ym(:)-mz(:));
%              [mk1,mk,K*(ym(:)-mz(:))]
             end
             %update
             Pk1=Pk-K*Pz*K';

end
