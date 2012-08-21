from numpy import empty
from numpy import array

def ReadObsData():
    fid=open("tst/JuneData.txt","r")
    infile=fid.readlines()

    nR=eval(infile[1])
    buf=infile[3]; buf=buf.split(); xkm=array(buf,float)
    buf=infile[5]; buf=buf.split(); L=array(buf,float)
    nt=eval(infile[7]);
    buf=infile[9]; buf=buf.split(); t=array([buf],float)
    hin=empty([nR,nt])
    for i in range(0,nR):
        buf=infile[i+11]; buf=buf.split(); hin[i,:]=array(buf,float)
    buf=infile[15]; buf=buf.split(); h0bar=array([buf],float)
    Sin=empty([nR,nt])
    for i in range(0,nR):
        buf=infile[i+17]; buf=buf.split(); Sin[i,:]=array(buf,float);
    Sin=Sin/1e5
    dA=empty([nR,nt])
    for i in range(0,nR):
        buf=infile[i+21]; buf=buf.split(); dA[i,:]=array(buf,float)
    wobs=empty([nR,nt])
    for i in range(0,nR):
        buf=infile[i+25]; buf=buf.split(); wobs[i,:]=array(buf,float)
    fid.close()   

    return nR,xkm,L,nt,t,hin,h0bar,Sin,dA,wobs
