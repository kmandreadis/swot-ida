function [Estimate,C] = CalculateEstimates (C,D,Obs,Prior)

% 1) estimates on the chain A0, n, q: means & covariances
Estimate.A0hat=mean(C.thetaA0(:,C.Nburn+1:end),2);
Estimate.CA0=cov(C.thetaA0(:,C.Nburn+1:end)');
Estimate.stdA0Post=sqrt(diag(Estimate.CA0));

Estimate.nhat=mean(C.thetan(:,C.Nburn+1:end),2);
Estimate.Cn=cov(C.thetan(:,C.Nburn+1:end)');
Estimate.stdnPost=sqrt(diag(Estimate.Cn));

Estimate.qhat=mean(C.thetaq(:,C.Nburn+1:end),2);
Estimate.Cq=cov(C.thetaq(:,C.Nburn+1:end)');
Estimate.stdqpost=sqrt(diag(Estimate.Cq));

%2) calculate the Q chain, and estimate mean and std
for i=1:C.N,
    C.thetaQ(:,:,i) = 1./(C.thetan(:,i)*ones(1,D.nt)) .* ...
        (C.thetaA0(:,i)*ones(1,D.nt)+Obs.dA).^(5/3).*Obs.w.^(-2/3).*sqrt(Obs.S) ;
end

Estimate.QhatPost=mean(C.thetaQ(:,:,C.Nburn+1:end),3);
Estimate.QstdPost=std(C.thetaQ(:,:,C.Nburn+1:end),[],3);

%3) Calculate Q prior estimate
Estimate.QhatPrior=1./(Prior.meann*ones(1,D.nt)) .* ...
    (Prior.meanA0*ones(1,D.nt)+Obs.dA).^(5/3).*Obs.w.^(-2/3).*sqrt(Obs.S);

%4)   Discharge error budget: all done for Q(nr x nt)
%4.1) Uncertainty estimate of the dA term
Obs.sigdAv=sqrt(diag(Obs.CdA));
Obs.sigdA=reshape(Obs.sigdAv,D.nt,D.nR)';
%4.2) Uncertainty (variance) of the Manning terms
Estimate.QhatUnc.n=(Estimate.stdnPost./Estimate.nhat).^2;
Estimate.QhatUnc.w=(2/3*Obs.sigw./Obs.w).^2;
Estimate.QhatUnc.A0=(5/3.*Estimate.stdA0Post*ones(1,D.nt)./(Estimate.A0hat*ones(1,D.nt)+Obs.dA)).^2;
Estimate.QhatUnc.dA=(5/3.*Obs.sigdA./(Estimate.A0hat*ones(1,D.nt)+Obs.dA)).^2;
Estimate.QhatUnc.S=(1/2.*Obs.sigS./Obs.S).^2;

%4.3) Estimate correlation coefficient between A0 & n
for i=1:D.nR,
    R=corrcoef([C.thetaA0(i,:); C.thetan(i,:);]');
    Estimate.rho_A0n(i,1)=R(1,2);
end

%4.4) Estimate uncertainty of Manning's Q
%4.4.1) Estimate uncertainty in Q due to cross-correlation of A0 & n
Estimate.QhatUnc.A0n=(-2*(5/3).*Estimate.rho_A0n*ones(1,D.nt) .* ...
    (Estimate.stdA0Post*ones(1,D.nt)) .* (Estimate.stdnPost*ones(1,D.nt)) ./ ...
    (Estimate.nhat*ones(1,D.nt)) ./ (Estimate.A0hat*ones(1,D.nt)+Obs.dA)  );

%4.4.2) Estimate total Q uncertinty
Estimate.QhatUnc.Hat=sqrt( Estimate.QhatUnc.n+mean(Estimate.QhatUnc.w,2)+...
    mean(Estimate.QhatUnc.A0,2)+mean(Estimate.QhatUnc.dA,2)+mean(Estimate.QhatUnc.S,2)+...
    mean(Estimate.QhatUnc.A0n,2)  );

%4.4.2) Estimate total Q uncertinty wrt time
Estimate.QhatUnc.HatAll=sqrt( Estimate.QhatUnc.n*ones(1,D.nt) + Estimate.QhatUnc.w+...
    Estimate.QhatUnc.A0+Estimate.QhatUnc.dA+Estimate.QhatUnc.S+...
    Estimate.QhatUnc.A0n  );

%4.4.3) Discharge error budget
Estimate.QerrVarSum(1:D.nR,1)=Estimate.QhatUnc.n;
Estimate.QerrVarSum(1:D.nR,2)=mean(Estimate.QhatUnc.A0,2);
Estimate.QerrVarSum(1:D.nR,3)=mean(Estimate.QhatUnc.dA,2);
Estimate.QerrVarSum(1:D.nR,4)=mean(Estimate.QhatUnc.w,2);
Estimate.QerrVarSum(1:D.nR,5)=mean(Estimate.QhatUnc.S,2);

return