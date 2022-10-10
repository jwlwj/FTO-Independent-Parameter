function [g,cvx_status,L]=cong(A,B,C,Dd,Cz,p,nd,Bd,e0,H)

[n,nu]=size(B);
ny=size(C,1);nz=size(Cz,1);

cvx_begin sdp quiet
cvx_precision best
variable P(n,n) symmetric ;
variable Y(n,ny);
variable muh(1);
variable S(ny,ny) diagonal;
variable M(ny,ny) diagonal;
variable g(1);

minimize g
subject to

P>=0;
muh>=0;

%Fault-free
[P*A+A'*P+Y*C+C'*Y', P*Bd+Y*Dd , Cz';
(P*Bd+Y*Dd)', -g*eye(nd), zeros(nd,nz);
Cz, zeros(nz,nd),-g*eye(nz)] <=0;

T1=[P*A+A'*P, P*Bd , Cz';
(P*Bd)', -1.01*g*eye(nd), zeros(nd,nz);
Cz, zeros(nz,nd),-1.01*g*eye(nz)];
% T1=[P*A+A'*P, P*Bd , Cz';
% (P*Bd)', -g*eye(nd), zeros(nd,nz);
% Cz, zeros(nz,nd),-g*eye(nz)];
T2=[C';Dd';zeros(nz,ny)];
T3=[Y',zeros(ny,nd),zeros(ny,nz)];

[T1+p*T2*(muh*eye(ny)-(M+M'))*T2', T3'+T2*S', T2*M*H, zeros(n+nd+nz,ny);
    (T3'+T2*S')', -(S+S'),-M*H,M*sqrt(p);
    (T2*M*H)',(-M*H)',-muh*eye(ny-1),zeros(ny-1,ny);
    zeros(ny,n+nd+nz),(M*sqrt(p))',zeros(ny,ny-1),-muh*eye(ny)]<=0;


% d = (0:(2^ ny) -1)';
% b = de2bi(d,'left-msb');
% 
% for i =1:(2^ny)
% Delta = diag(b(i,:));
% if e0'*Delta*e0 >=p
%  [P*A+A'*P+Y*Delta*C+C'*Delta'*Y', P*Bd+Y*Delta*Dd , Cz';
% (P*Bd+Y*Delta*Dd)', -g*eye(nd), zeros(nd,nz);
% Cz, zeros(nz,nd),-g*eye(nz)] <=0;
% 
% end
% end

cvx_end
L=P\Y;
end

