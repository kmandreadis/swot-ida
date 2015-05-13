function Truth = ReadTruth (fname,DAll)

fid=fopen(fname,'r');
fgetl(fid); Truth.A0=fscanf(fid,'%f',DAll.nR)'; fscanf(fid,'\n');
fgetl(fid); Truth.q=fscanf(fid,'%f',DAll.nR*(DAll.nt-1)); fscanf(fid,'\n');
fgetl(fid); Truth.n=fscanf(fid,'%f',DAll.nR); fscanf(fid,'\n');
fgetl(fid); QtAll=fscanf(fid,'%f',DAll.nR*DAll.nt); fscanf(fid,'\n');

for i=1:DAll.nR,
    Truth.Q(i,:)=QtAll( (i-1)*DAll.nt+1 : i*DAll.nt,:);
end

if feof(fid),
    fclose(fid);
    return
end

fgetl(fid); dAAll=fscanf(fid,'%f',DAll.nR*DAll.nt); fscanf(fid,'\n');
for i=1:DAll.nR,
    Truth.dA(i,:)=dAAll( (i-1)*DAll.nt+1 : i*DAll.nt,:);
end
Truth.dAv=reshape(Truth.dA',DAll.nR*DAll.nt,1);

if feof(fid),
    fclose(fid);
    return
end

fgetl(fid); hAll=fscanf(fid,'%f',DAll.nR*DAll.nt); fscanf(fid,'\n');
fgetl(fid); WAll=fscanf(fid,'%f',DAll.nR*DAll.nt); fscanf(fid,'\n');
for i=1:DAll.nR,
    Truth.h(i,:)=hAll( (i-1)*DAll.nt+1 : i*DAll.nt,:);
end

Truth.hv=reshape(Truth.h',DAll.nR*DAll.nt,1);

for i=1:DAll.nR,
    Truth.W(i,:)=WAll( (i-1)*DAll.nt+1 : i*DAll.nt,:);
end

Truth.Wv=reshape(Truth.W',DAll.nR*DAll.nt,1);

fclose(fid);

return