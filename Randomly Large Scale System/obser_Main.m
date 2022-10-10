clear;clc;

n=6;ny=10;nu=4;nz=3;p=0;nd=4 ;

[A,B,C,D]=sysgen(n,ny,nu);

Dd = randn(ny,nd);

Cz=randn(nz,n);

Bd = randn(n,nd);

e0= ones(ny ,1);
e = kron(ones(ny,1),[1;0]);
[L, R] = Permutations(ny);

sigma = 0.3;
theta =pi/2.15;
theta0 = pi/5;
r = 30;

tic
[g0,cvx_status0]=obser_Loop(A,B,C,Dd,Cz,p,nd,Bd,sigma,theta,theta0,r,e0);
toc

tic
[g1,cvx_status1]=obser_BlockHadamard_new(A,B,C,Dd,Cz,p,nd,Bd,sigma,theta,theta0,r,e0, L, R);
toc
[g0;g1]

