function [Estimate] = FilterEstimate(Estimate,C,D,Obs)

Deltax = CalcDelta(D.nR,D.nt,D.L);
DeltaA = CalcADelta(D.nR,D.nt);
B = CalcB(D.nR,D.nt);

H=-Deltax;
y=(DeltaA*Obs.hv) ./ D.dt .* (B*Obs.wv);

for i=1:C.N,
    Qchain(:,i)=reshape(C.thetaQ(:,:,i)',D.nR*D.nt,1);     
end

P=diag(var(Qchain,[],2));
R=Obs.CA;

K=P*H'/(H*P*H'+R);

for i=1:C.N,
    xminus(:,i)=Qchain(:,i);     
    xplus(:,i)=xminus(:,i)+K*(y-H*xminus(:,i));
end

Qhatfv=mean(xplus(:,C.Nburn:C.N),2);
Estimate.QhatPostf=reshape(Qhatfv,D.nt,D.nR)';

return