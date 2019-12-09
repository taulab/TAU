%{
   analytical part of calculations is displayed in comments
%}
syms s;

Pl = s + 2;
Pr = 1;
Rl = s + 1;
Rr = s - 3;
num1 = sym2poly(Pr*Pl);
den1 = sym2poly(Rl*Rr);
Wo = tf(num1, den1);

%{

Wo =
 
      s + 2
  -------------
  s^2 - 2 s - 3


tp = 1
Tp = 7.2

alpha = tp / Tp = 0.1389

G(q) = q^3 + 5.1*q^2 + 6.35*q + 1

a0 = alpha^3 = 0.0027
a1 = alpha^2 * 5.1 = 0.0984
a2 = alpha * 6.35 = 0.8819
a3 = 1

G(s) = a0*s^3 + a1*s^2 + a2*s + a3
 
G(s) =
 
(125*s^3)/46656 + (85*s^2)/864 + (127*s)/144 + 1

npl = 1
npr = 0
nrl = 1
nrr = 1
ng = 3

by solving a system of inequalities we get 
    nn = 0
    nm = 2

as a result, 
the functions N(s) and M(s) will take the following form:

N(s) = a0
M(s) = b0*s^2 + b1*s + b2

to find the unknown coefficients we need to solve this equiation:

    Pr*M(s) + Rr(s)*N*s^r == G(s)

b0*s^2 + b1*s + b2 + c0*s^2*(s - 3)=a0*s^3 + a1*s^2 + a2*s + a3

c0*s^3 + s^2(b0 - 3*c0) + b1*s + b2

=>
%}
r=2;
N = 0.0027;
M = (1655*s^2)/15552 + (127*s)/144 + 1;

%{
to find the transfer function we must use this formula:
    
    Wp(s)=(Rl(s)*M(s))/(Pl(s)*N*s^r)

%}

num2 = sym2poly(Rl*M);
den2 = sym2poly(Pl * N * s^r);
Wp = tf(num2, den2);

%{
Wp =
 
  0.1064 s^3 + 0.9884 s^2 + 1.882 s + 1
  -------------------------------------
       0.002679 s^3 + 0.005358 s^2


we need to convert both functions
from continuous time to discrete time:

%}
T = 0.01;

Woz = c2d(Wo, T);
Wpz = c2d(Wp, T);

%{
Wo(z) =
 
    0.0102 z - 0.01
  --------------------
  z^2 - 2.021 z + 1.02


Wp(z) =
 
  39.72 z^3 - 115.5 z^2 + 111.9 z - 36.1
  --------------------------------------
     z^3 - 2.98 z^2 + 2.96 z - 0.9802

%}
Weg = feedback(1, Wpz*Woz);
Wyg = feedback(Wpz*Woz, 1);

info = stepinfo(Wyg,'SettlingTimeThreshold',0.05);
settling_time = roundn(info.SettlingTime, -4);
display(settling_time);

subplot(3,1,1);
step(Wyg, 100*T);
t = 0:T:1;
g1 = ones(1, 1/T + 1);
g2 = 2*t + 1;
subplot(3,1,2);
lsim(Weg,g1,t);
subplot(3,1,3);
lsim(Weg,g2,t);

%{
            output:

settling_time =

    0.2400

%}