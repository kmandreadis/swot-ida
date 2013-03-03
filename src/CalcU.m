function U = CalcU (D)

M=D.nR*D.nt;
N=D.nR*(D.nt-1);

U=[  zeros(1,N);
     tril(ones(D.nt-1)) zeros(D.nt-1,2*(D.nt-1)) ;
     zeros(1,N);
     zeros(D.nt-1,D.nt-1) tril(ones(D.nt-1))   zeros(D.nt-1,D.nt-1)
     zeros(1,N);
     zeros(D.nt-1,2*(D.nt-1))  tril(ones(D.nt-1));    ];

 return