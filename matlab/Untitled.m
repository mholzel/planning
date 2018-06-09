clear;clc;close all

% Create the polynomials 
syms t T s_0 d_0
n = 5;
s = t.^(0:n) * [ s_0 ; sym( 's_%d', [n,1] ) ];
d = t.^(0:n) * [ d_0 ; sym( 'd_%d', [n,1] ) ];

% The velocity, acceleration, and jerk costs  
syms v_star a_star j_star a_v a_a a_j 
v_cost = exp( a_v * (             diff(s)^2 +             diff(d)^2 - v_star^2 ) );
a_cost = exp( a_a * (       diff(diff(s))^2 +       diff(diff(d))^2 - a_star^2 ) );
j_cost = exp( a_j * ( diff(diff(diff(s)))^2 + diff(diff(diff(d)))^2 - j_star^2 ) );

% syms p 
% v_cost = ((             diff(s)^2 +             diff(d)^2 ) / v_star^2)^p;
% a_cost = ((       diff(diff(s))^2 +       diff(diff(d))^2 ) / a_star^2)^p;
% j_cost = (( diff(diff(diff(s)))^2 + diff(diff(diff(d)))^2 ) / j_star^2)^p;
% int(v_cost,t,0,T)

% Collision avoidance cost
syms s_star d_star a_c
collision_cost = exp( a_c * (s-s_star) * (d-d_star));

% Lane keeping cost
beta = 2;
lane_cost = cos( pi * d / 4 )^2;

% Reference velocity cost 
syms v_ref sigma
reference_v_cost = 1 - exp( - (diff(s) - v_ref)^2 / sigma^2 );

% Total cost
total_cost =  v_cost ...
            + a_cost ...
            + j_cost ...
            + collision_cost ...
            + lane_cost ...
            + reference_v_cost;

total = simplify( total_cost );

% Next, integrate that cost over some interval of interest
syms T 
J = int(total, t, 0, T);
J