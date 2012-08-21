from numpy import array

def ReadTruth():
    # read Truth
    fid=open("tst/Truth.txt","r")
    infile=fid.readlines()
    buf=infile[1]; buf=buf.split(); A0true=array([buf],float)
    buf=infile[3]; buf=buf.split(); qtrue=array(buf,float)
    buf=infile[5]; buf=buf.split(); xQ=array(buf,float)
    buf=infile[7]; buf=buf.split(); tQ=array(buf,float)
    buf=infile[9]; buf=buf.split(); Qtrue=array(buf,float)
    fid.close() 

    return A0true,qtrue,xQ,tQ,Qtrue
