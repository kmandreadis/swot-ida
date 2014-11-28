function Prior=ProcessPrior(Prior,Obs,D)

w0=Obs.w(:,1); %width at baseflow
S0=Obs.S(:,1); %slope at baseflow

%this should come out very close to the calculation of the mean, below...
Prior.meanA0=(Prior.meanQbase.*Prior.meann.*w0.^(2/3).*S0.^(-.5)).^(3/5);

%crude monte-carlo estimate for the prior standard deviatin of A0

N=1E4;

v=(Prior.covQbase*Prior.meanQbase)^2;
[mu,sigma] = logninvstat(Prior.meanQbase,v);
Qb=lognrnd(mu,sigma,1,N);

n=randn(D.nR,N).*(Prior.stdn*ones(1,N))+Prior.meann*ones(1,N);
n(n<0)=.001;

for i=1:D.nR,
    A0(i,:)=( Qb.*n(i,:).*w0(i).^(2/3).*S0(i).^(-.5) ).^(3/5);
end

Prior.meanA0=mean(A0,2);
Prior.stdA0=std(A0,[],2);


return