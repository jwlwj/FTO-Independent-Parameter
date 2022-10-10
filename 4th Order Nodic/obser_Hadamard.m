function [g,cvx_status,L]=obser_Hadamard(A,B,C,Dd,Cz,p,nd,Bd,sigma,theta,theta0,r,e0)

[n,nu]=size(B);
ny=size(C,1);nz=size(Cz,1);
e1=ones(2*ny,1);

cvx_begin sdp quiet
cvx_precision best
variable P(n,n) symmetric ;
variable Y(n,ny);

variable S1(ny ,ny) diagonal ;
variable S2(2*ny ,2*ny) diagonal ;
variable S3(ny ,ny) diagonal ;

variable M1(ny ,ny) symmetric;
variable M2(2*ny ,2*ny) symmetric;
variable M3(ny ,ny) symmetric;

expression M1e(ny ,ny);
expression M1I(ny ,ny);
expression M2e(2*ny ,2*ny);
expression M2I(2*ny ,2*ny);
expression M3e(ny ,ny);
expression M3I(ny ,ny);

variable g(1);

minimize g
subject to

M1>=0;
M2>=0;
M3>=0;
P>=0;

M1I=M1.*eye(ny);
M1e=M1.*(e0*e0');
M2I=M2.*eye(2*ny);
M2e=M2.*(e1*e1');
M3I=M3.*eye(ny);
M3e=M3.*(e0*e0');

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

% --Constraints --
% sigma
T11 =P*A+(P*A)'+2*sigma*P;
T12 =C';
T13 =Y';

[T11+T12*M1I*(ny-p)*T12'-T12*M1e*T12', (T13+S1*T12'+ M1e'*T12')';
 T13+S1*T12'+ M1e'*T12', -(S1+S1')-M1e ] <=0;

% theta
T21 =[sin(theta)*(P*A+A'*P) ,cos(theta)*(P*A-A'*P);
    cos(theta)*(P*A-A'*P)', sin(theta)*(P*A+A'*P)];

T22 = [C', zeros(n,ny); zeros(n,ny),C'];

T23=[sin(theta)*Y',-cos(theta)*Y';
    cos(theta)*Y',sin(theta)*Y'];
 
 [T21+T22*M2I*(ny-p)*T22'-T22*(M2e)*T22', (T23+S2*T22'+M2e'*T22')';
 T23+S2*T22'+M2e'*T22', -(S2+S2')-M2e] <=0;

% r
T31 =[-r*P, P*A ; A'*P, -r*P];

T32 =[ zeros(n,ny);C'];

T33 =[Y' zeros(ny ,n)];

[T31+T32*M3I*(ny-p)*T32'-T32*M3e*T32', (T33+S3*T32'+ M3e'*T32')';
T33+S3*T32'+ M3e'*T32', -(S3+S3')-M3e ] <=0;

cvx_end
L=P\Y;
end
