% Some examples to illustrate the use of this code.

% Most of these examples have been coded for clarity rather than speed.
% There are many places where this can be made faster, especially in 
% calculating the early growth rate for a distribution of household sizes.


% Basic Parameters
alpha = 1;
beta = 3;
gamma = 0.5;
sigma = 0.5;


% 1. Calculate R* and r for the SEEIIR household model, assuming no
% antivirals and N=5.

[R,r] = calculate_no_antivirals(5,alpha,beta,gamma,sigma)


% 2. Calculate same as in 1., but for the model with exponential delay, and specified values
% of intervention efficacy and mean delay.

tau=0.3;
rho = 0.5;
delay = 1;

ra = 1/delay;

[R_exp,r_exp] = calculate_exp_delay(5,alpha,beta,gamma,sigma,tau,rho,ra)


% 3. Calculate same as in 2., but for the model with constant delay.

[R_const,r_const] = calculate_const_delay(5,alpha,beta,gamma,sigma,tau,rho,delay)


% 4. Calculate R* and r for a distribution of households,
%    - assuming no antivirals.

% household size distribution in the UK 2001.
hh_dist = [0.3028,0.3407,0.1551,0.1332,0.0488,0.0141,0.0053];

[R_dist1,r_dist1] = dist_no_antivirals(hh_dist,alpha,beta,gamma,sigma)


% 5. Calculate R* and r for a distribution of households,
%    - assuming antivirals with an exponential delay. 

[R_dist2,r_dist2] = dist_exp_delay(hh_dist,alpha,beta,gamma,sigma,tau,rho,ra)


% 6. Calculate R* and r for a distribution of households,
%    - assuming antivirals with a constant delay.

[R_dist3,r_dist3] = dist_const_delay(hh_dist,alpha,beta,gamma,sigma,tau,rho,delay)