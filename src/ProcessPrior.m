function [Prior,jmp,AllObs]=ProcessPrior(Prior,Obs,D,AllObs,jmp,DAll)

% calculate minimum A0 values -- this is for "All"
for i=1:D.nR,
    if min(AllObs.dA(i,:))>=0
        allA0min(i,1)=0;
    else
        allA0min(i,1)=-min(AllObs.dA(i,:));
    end
end

% calculate minimum values for A0
for i=1:D.nR,
    if min(Obs.dA(i,:))>=0
        estA0min(i,1)=0;
    else
        estA0min(i,1)=-min(Obs.dA(i,:));
    end
end

% shift that A0 into the estimate window
i1=find(DAll.t==D.t(1));
AllObs.A0Shift=AllObs.dA(:,i1); %this is area at first estimation time > than time 0

%for future reference save the more restrictive
jmp.A0min=max(allA0min+AllObs.A0Shift,estA0min);

% calculate mean A0, from the AllObs perspective
N=1E4;

v=(Prior.covQbar*Prior.meanQbar)^2;
[mu,sigma] = logninvstat(Prior.meanQbar,v);
Qbar=lognrnd(mu,sigma,1,N);

dAbar=mean(AllObs.dA,2);
Sbar=mean(AllObs.S,2);
wbar=mean(AllObs.w,2);

for i=1:N,
    A0bar(:,i)= ( Qbar(i).*Prior.meann.*wbar.^(2/3).*Sbar.^(-.5) ).^(3/5)-dAbar;
end

j=all(A0bar>allA0min*ones(1,N));

%switch everything to estimation window perspective
A0bar=A0bar+AllObs.A0Shift*ones(1,N);

Prior.meanA0=mean(A0bar(:,j),2);
Prior.stdA0=std(A0bar(:,j),[],2);

for i=1:N,
    Q(:,:,i) = 1./(Prior.meann*ones(1,D.nt)) .* ...
        (A0bar(:,i)*ones(1,D.nt)+Obs.dA).^(5/3).*Obs.w.^(-2/3).*sqrt(Obs.S) ;
    Qavg=squeeze(mean(Q(:,:,i),1));
    Qbase(i)=min(Qavg);
end

Prior.meanQbase=mean(Qbase(j));
Prior.stdQbase=std(Qbase(j));



return