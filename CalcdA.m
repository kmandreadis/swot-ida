%calculate dA from h & w observations

function [Obs] = CalcdA(D,Obs)


%Compare DeltaA
for r=1:D.nR,
    for t=1:D.nt-1,
        DeltaAHat(r,t)=(Obs.w(r,t+1)+Obs.w(r,t))/2*(Obs.h(r,t+1)-Obs.h(r,t));
    end
end

%now calculate a dAHat from the DeltaAHat..
Obs.dA=[zeros(D.nR,1) DeltaAHat*triu(ones(D.nt-1))];

Obs.dAv=reshape(Obs.dA',D.nR*D.nt,1);

return