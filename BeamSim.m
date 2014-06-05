% Euler Beam eqn simulations
ccc

P = 5;
[phi, D, x, w,gamma,lval] = GLL_Basis(P);
D1=D;
D2=D1*D1;
D3=D2*D1;
D4=D2*D2;

a = 0;
b = 2*pi;
J = (b-a)/2*x+a

u = sin(x);
plot(x,u)
plot(x,D4*u)
