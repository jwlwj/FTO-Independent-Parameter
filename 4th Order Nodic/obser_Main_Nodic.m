clear;clc;close all;
%% system Initialisation
load('Nord4stat_1in_6out');
B=Bu;

n = size (A, 2);
nu = size (B, 2);
ny = size (C, 1);
nd = size (Bd ,2) ;
nz = size (Cz , 1);
p=1;

e0= ones(ny ,1);
e = kron(ones(ny,1),[1;0]);
[L, R] = Permutations(ny);

sigma = 0.01;
theta = pi/2.05;
theta0 = pi/5;
r = 20;

%% Simulink Parameter
y2=[0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,1,0;0,0,0,0,0,0];
y1=eye(ny);
state3=[0,0,1,0];

%% Parameter needed by CongRPP
HH=ny*eye(ny)-e0*e0';
[T,Lambda]=eig(HH);
[y,i]=sort(diag(Lambda));
T=T(:,i);
Lambda=Lambda(i,i);
H=T(:,2:end)*diag(sqrt(diag((Lambda(2:end,2:end)))));

e1= ones(2*ny ,1);
HH1=ny*eye(2*ny)-e1*e1';
[T,Lambda]=eig(HH1);
[y,i]=sort(diag(Lambda));
T=T(:,i);
Lambda=Lambda(i,i);
H1=T(:,2:end)*diag(sqrt(diag((Lambda(2:end,2:end)))));

%% Result
[g0,cvx_status0,L0]=obser_Loop(A,B,C,Dd,Cz,p,nd,Bd,sigma,theta,theta0,r,e0);
[g1,cvx_status1,L1]=NFT(A,B,C,Dd,Cz,p,nd,Bd,sigma,theta,theta0,r,e0);
[g2,cvx_status2,L2]=obser_Hadamard(A,B,C,Dd,Cz,p,nd,Bd,sigma,theta,theta0,r,e0);
[g3,cvx_status3,L3]=obser_BlockHadamard_new(A,B,C,Dd,Cz,p,nd,Bd,sigma,theta,theta0,r,e0, L, R);
[gcr,cvx_status,L]=CongRPP(A,B,C,Dd,Cz,p,nd,Bd,sigma,theta,theta0,r,e0,H,H1);
[gc,cvx_statusc,Lc]=cong(A,B,C,Dd,Cz,p,nd,Bd,e0,H);

[g0;g1;g2;g3;gcr;gc]

