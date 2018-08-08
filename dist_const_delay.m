function [RR,rr] = dist_const_delay(hh_dist,alpha,beta,gamma,sigma,tau,rho,delay)
% Calculates R* (RR) and the early growth rate (rr) for the SEEIIR model,
%
% assuming a constant delay till antivirals, and
% assuming a distribution of household sizes.


Nmax = length(hh_dist);

% calculate the sized biased distribution.
sb_dist = (1:Nmax).*hh_dist/sum((1:Nmax).*hh_dist);


for ii=1:Nmax
    m(ii) = SEEIIR(ii);
end


% normalise beta correctly
beta_vec = ones(Nmax,1)*beta;
for nn=2:Nmax
    beta_vec(nn)=beta_vec(nn)/(nn-1);
end

rstar_hh = zeros(1,Nmax);

% calculate R*

for ii=1:Nmax
    
    v = zeros(length(m(ii).Cind),1);
    u = m(ii).totI(m(ii).Cind)*alpha;
    id = speye(length(m(ii).Cind));
    
    Q_pre = m(ii).input_params(beta_vec(ii),sigma,gamma,0);
    Q_post = m(ii).input_params(beta_vec(ii)*(1-tau)*(1-rho),sigma,gamma,0);
    
    part1 = phiv(delay,Q_pre(m(ii).Cind,m(ii).Cind),u,v);
    
    init_second = mexpv(delay,Q_pre',m(ii).init_cond);
    init_second =  init_second(m(ii).Cind);
    
    f = -m(ii).totI(m(ii).Cind)*alpha*(1-tau); % post admin, alpha is also reduced
    part2 = Q_post(m(ii).Cind,m(ii).Cind)\f;
    
    rstar_hh(ii) = part1(m(ii).Iind) + init_second'*part2;
    
end

RR = rstar_hh*sb_dist';


% If R* > 0 then calculate the early growth rate.
if RR > 0
    rr = fzero( @(r)eg_fun(r), [0, 6]);
else
    rr = 0;
end


    function [result] = eg_fun(r)
        % the function we want to minimize to find r
        
        part4 = zeros(1,Nmax);
        
        for jj=1:Nmax
            
            v = zeros(length(m(jj).Cind),1);
            u = m(jj).totI(m(jj).Cind)*alpha;
            id = speye(length(m(jj).Cind));
            
            Q_pre = m(jj).input_params(beta_vec(jj),sigma,gamma,0);
            Q_post = m(jj).input_params(beta_vec(jj)*(1-tau)*(1-rho),sigma,gamma,0);
            
            part1_mat = Q_pre(m(jj).Cind,m(jj).Cind)-id*r;
            
            part1 = phiv(delay,part1_mat,u,v);
            
            init_second = mexpv(delay,Q_pre',m(jj).init_cond);
            init_second =  init_second(m(jj).Cind);
            
            part2_mat = Q_post(m(jj).Cind,m(jj).Cind)-id*r;
            
            f = -m(jj).totI(m(jj).Cind)*alpha*(1-tau); % post admin, alpha is also reduced
            part2 = part2_mat\f;
            
            part4(jj) = part1(m(jj).Iind) + (init_second'*part2)*exp(-r*delay);
            
        end
        
        result = part4*sb_dist' - ((r+2*sigma)/(2*sigma))^2;
        
    end


end
