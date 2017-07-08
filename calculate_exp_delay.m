function [RR,rr] = calculate_exp_delay(N,alpha,beta,gamma,sigma,tau,rho,ra)
% Calculate R* (RR) and the early growth rate (rr) for the SEEIIR model,
% assuming an exponentially distributed delay.


m = SEEIIR_exp(N);


QQ1 = m.input_params(beta/(m.N-1),sigma,gamma,eps,tau,rho,ra);
Qc = QQ1(m.Cind,m.Cind);
id = speye(length(m.Cind));
alpha_vector = [alpha*ones(1,length(m.Cind)/2) alpha*(1-tau)*ones(1,length(m.Cind)/2)];
f = -m.I_vector(m.Cind).*alpha_vector';
sol=Qc\f;
RR = sol(m.Iind);


% If R* > 0 then calculate the early growth rate.
if RR > 0
    rr = fzero( @(r)eg_fun(r), [0, 6]);
else
    rr = 0;
end


    function [result] = eg_fun(r)
        % the function we want to minimize to find r
        
        part1_mat = Qc - id*r;
        part2 = part1_mat\f;
        
        result = part2(m.Iind)-1;
        
    end


end