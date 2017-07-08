function [RR,rr] = dist_exp_delay(hh_dist,alpha,beta,gamma,sigma,tau,rho,ra)
% Calculates R* (RR) and the early growth rate (rr) for the SEEIIR model,
%
% assuming an exponentially distributed delay till intervention, and
% assuming a distribution of household sizes.


Nmax = length(hh_dist);

% calculate the sized biased distribution.
sb_dist = (1:Nmax).*hh_dist/sum((1:Nmax).*hh_dist);


for ii=1:Nmax
    m(ii) = SEEIIR_exp(ii);
end


% normalise beta correctly
beta_vec = ones(Nmax,1)*beta;
for nn=2:Nmax
    beta_vec(nn)=beta_vec(nn)/(nn-1);
end

rstar_hh = zeros(1,Nmax);

% calculate R*

for ii=1:Nmax
    
    QQ1 = m(ii).input_params(beta_vec(ii),sigma,gamma,0,tau,rho,ra);
    
    Qc = QQ1(m(ii).Cind,m(ii).Cind);
    
    alpha_vector = [alpha*ones(1,length(m(ii).Cind)/2) alpha*(1-tau)*ones(1,length(m(ii).Cind)/2)];
    
    f = -m(ii).I_vector(m(ii).Cind).*alpha_vector';
    
    sol=Qc\f;
    
    rstar_hh(ii) = sol(m(ii).Iind);
    
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
        
        part3 = zeros(1,Nmax);
        
        for jj=1:Nmax
            
            QQ1 = m(jj).input_params(beta_vec(jj),sigma,gamma,0,tau,rho,ra);
            
            part1_mat = QQ1(m(jj).Cind,m(jj).Cind) - speye(length(m(jj).Cind))*r;
            
            alpha_vector = [alpha*ones(1,length(m(jj).Cind)/2) alpha*(1-tau)*ones(1,length(m(jj).Cind)/2)];
            
            part2 = part1_mat\(-m(jj).I_vector(m(jj).Cind).*alpha_vector');
            
            part3(jj) = part2(m(jj).Iind);
            
        end
        
        result = part3*sb_dist' - 1;
        
    end


end