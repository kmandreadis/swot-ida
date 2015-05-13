function [Prior,jmp,AllObs]=ProcessPriorOld(Prior,Obs,D,AllObs,jmp,DAll)

%some of this is less than ideal...

jmp.A0min=ceil(-min(AllObs.dA,[],2));
if max(AllObs.dAv) < max(jmp.A0min),
    refA0max=max(5*jmp.A0min);
else
    refA0max=max(AllObs.dAv);
end
refA0=linspace(max(jmp.A0min),refA0max,100);

nbar=mean(Prior.meann);

for i=1:length(refA0),
    Qhat(i)=mean(1./nbar.*(refA0(i)+AllObs.dAv).^(5/3).*AllObs.wv.^(-2/3).*sqrt(AllObs.Sv));
end

p=polyfit(Qhat,refA0,2);

N=1E4;

v=(Prior.covQbar*Prior.meanQbar)^2;
[mu,sigma] = logninvstat(Prior.meanQbar,v);
Qbar=lognrnd(mu,sigma,1,N);

A0=polyval(p,Qbar);

i1=find(DAll.t==D.t(1));
AllObs.A0Shift=AllObs.dA(:,i1); %this is area at first estimation time > than time 0
A0e=AllObs.A0Shift*ones(1,N)+ones(D.nR,1)*A0;  %A0 with reference to estimation

for i=1:N,
    A(:,:,i)=A0e(:,i)*ones(1,D.nt)+Obs.dA;
    Q(:,:,i) = 1./(Prior.meann*ones(1,D.nt)) .* ...
        A(:,:,i).^(5/3).*Obs.w.^(-2/3).*sqrt(Obs.S) ;
    Qavg=squeeze(mean(Q(:,:,i),1));
    Qbase(i)=min(Qavg);
    QbarHat(i)=mean(Qavg);
end

for i=1:N,
    A0v(:,i)= ( Qbase(i).*Prior.meann.*Obs.w(:,1).^(2/3).*Obs.S(:,1).^(-.5) ).^(3/5);
end

j=squeeze(all(all(A>0,1),2));
Prior.meanQbase=mean(Qbase(j));
Prior.stdQbase=std(Qbase(j));
 
Prior.meanA0=mean(A0v(:,j),2);
Prior.minA0=min(Obs.dA,[],2);
Prior.stdA0=std(A0v(:,j),[],2);

%% transfer the minimum from all to estimate
jmp.A0min=ceil(jmp.A0min+AllObs.A0Shift);

return