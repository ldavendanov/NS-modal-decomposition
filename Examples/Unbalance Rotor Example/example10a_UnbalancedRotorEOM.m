clear
close all
clc

syms t
syms l1 l2 m1 m2 k1x k1y k2
syms z1x(t) z1y(t) z2(t) theta(t)
d1 = l1+z1y;
d2 = l2+z2;


p1 = [z1x; l1+z1y];
p2 = p1 + [d2*cos(theta); d2*sin(theta)];

p1_dot = diff(p1);
p2_dot = diff(p2);

T = (1/2)*m1*(p1_dot.'*p1_dot) + (1/2)*m2*(p2_dot.'*p2_dot);
V = (1/2)*k1x*z1x.^2 + (1/2)*k1y*z1y.^2 + (1/2)*k2*z2.^2; 
L = T - V;

D0 = diff( diff( L, diff(z1x,t) ), t );
D1 = diff( diff( L, diff(z1y,t) ), t );
D2 = diff( diff( L, diff(z2,t) ), t );
D3 = diff( diff( L, diff(theta,t) ), t );

E0 = diff( L, z1x );
E1 = diff( L, z1y );
E2 = diff( L, z2 );
E3 = diff( L, theta );

eq0 = simplify( D0 - E0 );
eq1 = simplify( D1 - E1 );
eq2 = simplify( D2 - E2 );
eq3 = simplify( D3 - E3 );

M = [diff(eq0,diff(z1x,t,t)) diff(eq0,diff(z1y,t,t)) diff(eq0,diff(z2,t,t)) diff(eq0,diff(theta,t,t))
     diff(eq1,diff(z1x,t,t)) diff(eq1,diff(z1y,t,t)) diff(eq1,diff(z2,t,t)) diff(eq1,diff(theta,t,t))
     diff(eq2,diff(z1x,t,t)) diff(eq2,diff(z1y,t,t)) diff(eq2,diff(z2,t,t)) diff(eq2,diff(theta,t,t))
     diff(eq3,diff(z1x,t,t)) diff(eq3,diff(z1y,t,t)) diff(eq3,diff(z2,t,t)) diff(eq3,diff(theta,t,t))];

f = simplify( M*[diff(z1x,t,t); diff(z1y,t,t); diff(z2,t,t); diff(theta,t,t)] - [eq0; eq1; eq2; eq3 ] );

% M = [diff(eq1,diff(z1y,t,t)) diff(eq1,diff(z2,t,t)) diff(eq1,diff(theta,t,t))
%      diff(eq2,diff(z1y,t,t)) diff(eq2,diff(z2,t,t)) diff(eq2,diff(theta,t,t))
%      diff(eq3,diff(z1y,t,t)) diff(eq3,diff(z2,t,t)) diff(eq3,diff(theta,t,t))];
% 
% f = simplify( M*[diff(z1y,t,t); diff(z2,t,t); diff(theta,t,t)] - [eq1; eq2; eq3 ] );