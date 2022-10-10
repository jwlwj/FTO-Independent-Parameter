function [g,cvx_status,L]=CongRPP(A,B,C,Dd,Cz,p,nd,Bd,sigma,theta,theta0,r,e0,H,H1)

[n,nu]=size(B);
ny=size(C,1);nz=size(Cz,1);

cvx_begin sdp quiet
cvx_precision best
variable P(n,n) symmetric ;
variable Y(n,ny);
variable muh(1);
variable S1(ny ,ny) diagonal ;
variable S2(2*ny ,2*ny) diagonal;
variable S3(ny ,ny) diagonal ;
variable M1(ny,ny) diagonal;
variable M2(2*ny ,2*ny) diagonal;
variable M3(ny,ny) diagonal;
variable g(1);

minimize g
subject to

P>=0;
muh>=0;

%Fault-free
[P*A+A'*P+Y*C+C'*Y', P*Bd+Y*Dd , Cz';
(P*Bd+Y*Dd)', -g*eye(nd), zeros(nd,nz);
Cz, zeros(nz,nd),-g*eye(nz)] <=0;

%sigma
(P*A+A'*P+Y*C+C'*Y')+2*sigma*P <=0;
%theta
[sin(theta0)*(P*A+A'*P+Y*C+C'*Y'),cos(theta0)*(P*A+Y*C-A'*P-C'*Y');
(cos( theta0)*(P*A+Y*C-A'*P-C'*Y'))' ,sin(theta0)*(P*A+A'*P+Y*C+C'*Y') ] <=0;
%r
[-r*P, (P*A+Y*C);
(P*A+Y*C)', -r*P] <=0;

% --Constraints --

%sigma
T11 =P*A+(P*A)'+2*sigma*P;
T12 =C';
T13 =Y';

[T11+p*T12*(muh*eye(ny)-(M1+M1'))*T12', T13'+T12*S1', T12*M1*H, zeros(n,ny);
    (T13'+T12*S1')', -(S1+S1'),-M1*H,M1*sqrt(p);
    (T12*M1*H)',(-M1*H)',-muh*eye(ny-1),zeros(ny-1,ny);
    zeros(ny,n),(M1*sqrt(p))',zeros(ny,ny-1),-muh*eye(ny)]<=0;

% theta
T21 =[sin(theta)*(P*A+A'*P) ,cos(theta)*(P*A-A'*P);
    cos(theta)*(P*A-A'*P)', sin(theta)*(P*A+A'*P)];

T22 = [C', zeros(n,ny); zeros(n,ny),C'];

T23=[sin(theta)*Y',-cos(theta)*Y';
    cos(theta)*Y',sin(theta)*Y'];

[T21+p*T22*(muh*eye(2*ny)-(M2+M2'))*T22', T23'+T22*S2', T22*M2*H1, zeros(2*n,2*ny);
    (T23'+T22*S2')', -(S2+S2'),-M2*H1,M2*sqrt(p);
    (T22*M2*H1)',(-M2*H1)',-muh*eye(2*ny-1),zeros(2*ny-1,2*ny);
    zeros(2*ny,2*n),(M2*sqrt(p))',zeros(2*ny,2*ny-1),-muh*eye(2*ny)]<=0;

% r
T31 =[-r*P, P*A ; A'*P, -r*P];

T32 =[ zeros(n,ny);C'];

T33 =[Y' zeros(ny ,n)];

[T31+p*T32*(muh*eye(ny)-(M3+M3'))*T32', T33'+T32*S3', T32*M3*H, zeros(2*n,ny);
    (T33'+T32*S3')', -(S3+S3'),-M3*H,M3*sqrt(p);
    (T32*M3*H)',(-M3*H)',-muh*eye(ny-1),zeros(ny-1,ny);
    zeros(ny,2*n),(M3*sqrt(p))',zeros(ny,ny-1),-muh*eye(ny)]<=0;

cvx_end
L=P\Y;
end

