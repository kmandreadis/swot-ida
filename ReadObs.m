function [Domain,Obs] = ReadObs(fname)

fid=fopen(fname,'r');
fgetl(fid); Domain.nR=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Domain.xkm=fscanf(fid,'%f',Domain.nR); fscanf(fid,'\n');
fgetl(fid); Domain.L=fscanf(fid,'%f',Domain.nR)'; fscanf(fid,'\n'); %reach lengths
fgetl(fid); Domain.nt=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Domain.t=fscanf(fid,'%f',Domain.nt)'; fscanf(fid,'\n');
fgetl(fid); 

Domain.dt=reshape( (diff(Domain.t)'.*86400*ones(1,Domain.nR)),Domain.nR*(Domain.nt-1),1);

for i=1:Domain.nR, %read heights
    Obs.h(i,:)=fscanf(fid,'%f',Domain.nt); fscanf(fid,'\n');
end

fgetl(fid); Obs.h0=fscanf(fid,'%f',Domain.nR)'; fscanf(fid,'\n');

fgetl(fid); 
for i=1:Domain.nR, %read slopes
    Obs.S(i,:)=fscanf(fid,'%f',Domain.nt)./1E5; fscanf(fid,'\n');    %cm/km -> m/m
end
% fgetl(fid); 
% for i=1:Domain.nR, %read change in cross-sectional areas
%     Obs.dA(i,:)=fscanf(fid,'%f',Domain.nt); fscanf(fid,'\n');
%     
% end
fgetl(fid); 
for i=1:Domain.nR, 
    Obs.w(i,:)=fscanf(fid,'%f',Domain.nt); fscanf(fid,'\n');
end

fgetl(fid); Obs.sigS=fscanf(fid,'%f',1)/1E5; fscanf(fid,'\n');  %cm/km -> m/m
fgetl(fid); Obs.sigh=fscanf(fid,'%f',1)/1E2; fscanf(fid,'\n');  %cm -> m
fgetl(fid); Obs.sigw=fscanf(fid,'%f',1); fscanf(fid,'\n');  %,
fgetl(fid); Obs.sigdA=fscanf(fid,'%f',1); fscanf(fid,'\n');  %,

fclose(fid);


Obs.hv=reshape(Obs.h',Domain.nR*Domain.nt,1);
% Obs.dAv=reshape(Obs.dA',Domain.nR*Domain.nt,1);
Obs.wv=reshape(Obs.w',Domain.nR*Domain.nt,1);
Obs.Sv=reshape(Obs.S',Domain.nR*Domain.nt,1);

return