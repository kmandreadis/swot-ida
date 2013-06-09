function Truth = ReadTruth (fname,D)

fid=fopen(fname,'r');
fgetl(fid); Truth.A0=fscanf(fid,'%f',D.nR)'; fscanf(fid,'\n');
fgetl(fid); Truth.q=fscanf(fid,'%f',D.nR*(D.nt-1)); fscanf(fid,'\n');
fgetl(fid); Truth.n=fscanf(fid,'%f',D.nR); fscanf(fid,'\n');
fgetl(fid); QtAll=fscanf(fid,'%f',D.nR*D.nt); fscanf(fid,'\n');

for i=1:D.nR,
    Truth.Q(i,:)=QtAll( (i-1)*D.nt+1 : i*D.nt,:);
end

if feof(fid),
    fclose(fid);
    return
end

fgetl(fid); dAAll=fscanf(fid,'%f',D.nR*D.nt); fscanf(fid,'\n');
for i=1:D.nR,
    Truth.dA(i,:)=dAAll( (i-1)*D.nt+1 : i*D.nt,:);
end
Truth.dAv=reshape(Truth.dA',D.nR*D.nt,1);

if feof(fid),
    fclose(fid);
    return
end

fgetl(fid); hAll=fscanf(fid,'%f',D.nR*D.nt); fscanf(fid,'\n');
fgetl(fid); WAll=fscanf(fid,'%f',D.nR*D.nt); fscanf(fid,'\n');
for i=1:D.nR,
    Truth.h(i,:)=hAll( (i-1)*D.nt+1 : i*D.nt,:);
end

Truth.hv=reshape(Truth.h',D.nR*D.nt,1);

for i=1:D.nR,
    Truth.W(i,:)=WAll( (i-1)*D.nt+1 : i*D.nt,:);
end

Truth.Wv=reshape(Truth.W',D.nR*D.nt,1);

fclose(fid);

return