function [Obs,Prior] = GetCovMats(D,Obs,Prior)

%matrix sizes
M=D.nR*D.nt;
N=D.nR*(D.nt-1);

%matrices for operating on h & w
DeltaA = CalcADelta(D.nR,D.nt);
B = CalcB(D.nR,D.nt);

%Define covariance matrices of the raw observations: assume no
%  cross-correlation, for now
Obs.Ch=Obs.sigh^2.*eye(M);
Obs.Cw=Obs.sigw^2.*eye(M);
Obs.Cs=Obs.sigS^2.*eye(M);

%This is the Jacobian of the dAdt term w.r.t. h 
Obs.JAh=( B*Obs.wv*ones(1,M) ).*DeltaA ./ (D.dt*ones(1,M));
%This is the Jacobian of the dAdt term w.r.t. w
Obs.JAw=(DeltaA*Obs.hv*ones(1,M)) .* B ./ (D.dt*ones(1,M));
%combine terms to get the total covariance matrix of the dAdt term
Obs.JA=[Obs.JAh Obs.JAw]; 

%now calculate the covariance of the dAdt term w.r.t w&h
Chw=[Obs.Ch zeros(M);
     zeros(M) Obs.Cw;];
Obs.CA=Obs.JA*Chw*Obs.JA';        

%calculate covariance matrix of the dA term based on width & height errors
U = CalcU (D);
JdAh=U* ((ones(D.nR*(D.nt-1),1)*Obs.wv').*DeltaA);
JdAw=U* ((DeltaA*Obs.hv*ones(1,M)) .* B);
JdA=[JdAh JdAw];
Obs.CdA=JdA*Chw*JdA';

%Additional matrices used in the likelihood function: this for the dQ/dx term
Obs.CSdAw=[Obs.Cs zeros(D.nR*D.nt)  zeros(D.nR*D.nt);
           zeros(D.nR*D.nt) Obs.CdA zeros(D.nR*D.nt);
           zeros(D.nR*D.nt) zeros(D.nR*D.nt) Obs.Cw ;];

%This is the covariance of the lateral inflows
Prior.Cqf=eye(N).*Prior.stdq.^2;

return