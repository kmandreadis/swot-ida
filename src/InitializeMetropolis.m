function [Delta,DeltaA,B,C,thetauA0,thetaun,thetauq,thetauQb,R] = ...
               InitializeMetropolis (D,C,P,R)

Delta = CalcDelta(D.nR,D.nt,D.L);
DeltaA = CalcADelta(D.nR,D.nt);
B = CalcB(D.nR,D.nt);
              
%3.1) allocations, and initial state set
C.thetaA0=nan(D.nR,C.N); %this is just A0
C.thetaA0(:,1)=P.meanA0;
thetauA0=C.thetaA0(:,1);

C.thetan=nan(D.nR,C.N);
C.thetan(:,1)=P.meann;
thetaun=C.thetan(:,1);

C.thetaq=nan(D.nR*(D.nt-1),C.N);
C.thetaq(:,1)=P.meanq;
thetauq=C.thetaq(:,1);

C.thetaQb=nan(D.nR,C.N);
C.thetaQb(:,1)=P.meanQbase*ones(D.nR,1);
thetauQb=C.thetaQb(:,1);

%Get Random numbers
rng(R.Seed)
R.z1=randn(D.nR,C.N); 
R.z2=randn(D.nR,C.N);
R.z3=randn(D.nR*(D.nt-1),C.N);

R.u1=rand(C.N,1); %used for acceptance of A0
R.u2=rand(C.N,1); %used for acceptance of n
R.u3=rand(C.N,1); %used for acceptance of q

return