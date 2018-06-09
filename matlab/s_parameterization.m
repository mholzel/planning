clear;clc;close all

syms t v x0 v0 a0 j0 alpha beta b gamma g eta e

% s should have a parameterization that increases approximately linearly in
% time since we want ds/dt (that is, the velocity) to be constant in the 
% limit. 

% The following achieves that so long as the free variables b and g are positive. 

% First, note that we can either make the parameterization 
% in terms of a generic delta parameter, or we can rewrite it in terms of
% the final desired velocity, v. If you want the generic parameterization,
% simply comment this line out. However, it is typically easier to specify
% a good initial guess in terms of the velocity, so I would recommend that
% parameterization
% delta = -(a0 + v0*p_beta + v0*p_gamma - x*p_beta*p_gamma + x0*p_beta*p_gamma)/((p_beta - 1)*(p_gamma - 1));

% Next, we want to make sure that the initial conditions x0, v0, and a0 are 
% met. The following choices do that 
gamma = a0 / g^2;
alpha = x0 - gamma;
eta = v;
beta = v0 - eta + g * gamma;
% gamma = (a0 + delta - delta*p_beta + v0*p_beta)/(p_gamma*(p_beta - p_gamma));
% gamma = -(a0 - p_beta*v + p_beta*v0)/(p_gamma*(p_beta - p_gamma));
% beta = -a0 / (2*b);
f = alpha + t * ( eta + beta * exp( -b * t^2 ) ) + gamma * exp( -g * t );

f = simplify(f)
df = simplify(diff(f,t))
ddf = simplify(diff(df,t))
% dddf = simplify(diff(ddf,t))

f0 = simplify(subs(f,t,0))
df0 = simplify(subs(df,t,0))
ddf0 = simplify(subs(ddf,t,0))
% dddf0 = simplify(subs(dddf,t,0))

i = 1;
pre = 2;
fi = vpa(subs(f,t,i), pre)
dfi = vpa(subs(df,t,i), pre)
ddfi = vpa(subs(ddf,t,i), pre)
% dddfi = vpa(subs(dddf,t,i), pre)
