clear all
close all
clc
format longG 

function b = TRUSTb(xi,D,D2,alpha)
eps = 0.01;
xbar = 0.2*sqrt(0.02*(4-alpha)/(2-alpha));
Abar = (2-alpha)^2.5/(4-alpha)^1.5/0.02^1.5/0.2/2;
b = 0.1*D + 0.5*D2*(0.15^2+2*Abar*eps^(2-alpha)/(2-alpha)*xbar^(1+alpha))...
    + Abar*(integral(@(x) (xi(x)-D*x).*(abs(x/xbar)).^(-alpha-1),-xbar,-eps)...
    + integral(@(x) (xi(x)-D*x).*(abs(x/xbar)).^(-alpha-1),eps,xbar));
end