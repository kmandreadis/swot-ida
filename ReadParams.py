from numpy import array

def ReadParams():
    # read June parameters
    fid=open("tst/JuneParams.txt","r")
    infile=fid.readlines()

    N=eval(infile[1])
    Nburn=eval(infile[3])
    buf=infile[5]; buf=buf.split(); meanA0=array([buf],float); 
    buf=infile[7]; buf=buf.split(); stdA0=array(buf,float)
    meann=eval(infile[9])
    stdn=eval(infile[11])
    buf=infile[13]; buf=buf.split(); meanq=array([buf],float)
    stdq=eval(infile[15])
    sigS=eval(infile[17])/1e5
    sigh=eval(infile[19])/1e2

    class Jump:
        def __init__(self,stdA0burn,stdA0sim,A0min,stdn,nmin,stdq,qmin):
            self.stdA0burn=stdA0burn
            self.stdA0sim=stdA0sim
            self.A0min=A0min
            self.stdn=stdn
            self.nmin=nmin
            self.stdq=stdq
            self.qmin=qmin
        def setStdA0(self,stdA0):
            self.stdA0=stdA0   

    stdA0burn=eval(infile[21]); stdA0sim=eval(infile[23]); 
    A0min=eval(infile[25]); stdn1=eval(infile[27]);
    nmin=eval(infile[29]); stdq1=eval(infile[31]); qmin=eval(infile[33]);

    jmp=Jump(stdA0burn,stdA0sim,A0min,stdn1,nmin,stdq1,qmin)
    del stdA0burn,stdA0sim,A0min,stdn1,nmin,stdq1,qmin
    fid.close()

    return N,Nburn,meanA0,stdA0,meann,stdn,meanq,stdq,sigS,sigh,jmp
