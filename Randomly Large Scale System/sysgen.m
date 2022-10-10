function [A,B,C,D]=sysgen(n,p,m)
%      [A,B,C,D]=sysgen(n,p,m)
%      Generates a random stable pxm transfer matrix of order n.
%      A is nxn  and stable
%      B is nxm
%      C is pxn
%      D is pxm

M1=randn(n);
M2=randn(n);
M3=randn(n);
A=M3*(-M1'*M1+(M2-M2'))/M3;
B=randn(n,m);
C=randn(p,n);
D=rand(p,m);