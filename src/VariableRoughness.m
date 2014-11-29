function [nv,nvar]=VariableRoughness(Obs,D,c2,n)

deltah=.1;  %temporary -- should eventually estimate this from A0,W0
hmin=min(Obs.h,[],2);
Depth=Obs.h-hmin*ones(1,D.nt)+deltah;
Dbar=mean(Depth,2)*ones(1,D.nt);
Wbar=mean(Obs.w,2)*ones(1,D.nt);
Avar=(Obs.w.*Depth)./Wbar./Dbar;
c1=.85;
for i=1:D.nR
    nvar(i,:)=c1.*Avar(i,:).^c2(i).*n(i);
end

nv=reshape(nvar',D.nR*D.nt,1);

return