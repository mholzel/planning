clear;clc;close all

syms t x x0 v0 a0 j0 alpha beta b gamma g eta e 

% d should have a bounded parameterization, f. 
% The following achieves that so long as the free variables b and g

% First, note that we can either make the parameterization 
% in terms of a generic delta parameter, or we can rewrite it in terms of
% the final desired x-coordinate. If you want the generic parameterization,
% simply comment this line out. However, it is typically easier to specify
% a good initial guess in terms of x, so I would recommend that
% parameterization
% eta = -(a0 + v0*p_beta + v0*p_gamma - x*p_beta*p_gamma + x0*p_beta*p_gamma)/((p_beta - 1)*(p_gamma - 1));

% Next, we want to make sure that the initial conditions x0, v0, and a0 are 
% met. The following choices do that 
alpha = x;
beta = -v0 / b;
gamma = -(a0 + b*v0)/(2*g);
eta = x0 - alpha - beta - gamma;
% beta = -(a0 + delta - delta*p_gamma + v0*p_gamma)/(p_beta*(p_beta - p_gamma));
% gamma = (a0 + delta - delta*p_beta + v0*p_beta)/(p_gamma*(p_beta - p_gamma));
f = alpha + beta * exp( -b * t ) + gamma * exp( -g * t^2 ) + eta * exp( - e * t^3 );

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
