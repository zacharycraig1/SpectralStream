clc, clear, close,
n = 64;
delta = 20 / n;

x = linspace(-10,10,n+1);
y = linspace(-10,10,n+1);
tspan = 0:.5:4;

load('Afull1.mat');
load('Bfull1.mat');
load('Cfull1.mat');

A = (1/delta^2).*Afull; A(1,1) = 2;
B = (1/(delta*2)).*Bfull;
C = (1/(delta*2)).*Cfull;
A = sparse(A);
B = sparse(B);
C = sparse(C);

[X,Y] = meshgrid(x(1:n),y(1:n));
w0 = exp(-X.^2 - (Y.^2/20));

tic;[~,w] = ode45(@(t,w) system_rhs(t,w,[],A,B,C),tspan,w0);toc;


[L,U,P] = lu(A);
tic;[~,omega] = ode45(@(t,y) system_rhs_lu(t,y,[],A,B,C,L,U,P),tspan,w0);toc;

tic;[t,omg] = ode45(@(t,w) system_rhs_spectral(t,w,[],A,B,C),tspan,w0);toc;



tic;[t,omggmres] = ode45(@(t,w) sys_rhs_gmres(t,w,[],A,B,C), tspan,w0);toc;
tic;[t,omgbic] = ode45(@(t,w) sys_rhs_bic(t,w,[],A,B,C), tspan,w0);toc;
A1 = w;
A2 = omega;
A3 = omgbic;
A4 = omggmres;
A5 = omg;