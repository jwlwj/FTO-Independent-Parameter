function [g,cvx_status,L]=NFT(A,B,C,Dd,Cz,p,nd,Bd,sigma,theta,theta0,r,e0)

[n,nu]=size(B);
ny=size(C,1);nz=size(Cz,1);

cvx_begin sdp quiet
cvx_precision best
variable P(n,n) symmetric ;
variable Y(n,ny);
variable g(1);

minimize g
subject to

P>=0;

%Fault-free
[P*A+A'*P+Y*C+C'*Y', P*Bd+Y*Dd , Cz';
(P*Bd+Y*Dd)', -g*eye(nd), zeros(nd,nz);
Cz, zeros(nz,nd),-g*eye(nz)] <=0;

%sigma
(P*A+A'*P+Y*C+C'*Y')+2*sigma*P <=0;
%theta
[sin(theta0)*(P*A+A'*P+Y*C+C'*Y'),cos(theta0)*(-P*A-Y*C+A'*P+C'*Y');
(cos( theta0)*(-P*A-Y*C+A'*P+C'*Y'))' ,sin(theta0)*(P*A+A'*P+Y*C+C'*Y') ] <=0;
%r
[-r*P, (P*A+Y*C);
(P*A+Y*C)', -r*P] <=0;

cvx_end
L=P\Y;
end
