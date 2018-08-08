function [RR,rr] = calculate_const_delay(N,alpha,beta,gamma,sigma,tau,rho,delay)
% Calculate R* (RR) and the early growth rate (rr) for the SEEIIR model,
% assuming a constant delay until the allocation of antivirals.


m = SEEIIR(N);

v = zeros(length(m.Cind),1);
u = m.totI(m.Cind)*alpha;
id = speye(length(m.Cind));

Q_pre = m.input_params(beta/(m.N-1),sigma,gamma,0);
Q_post = m.input_params(beta/(m.N-1)*(1-tau)*(1-rho),sigma,gamma,0);


% calculate R* and check it's greater than 1

part1 = phiv(delay,Q_pre(m.Cind,m.Cind),u,v);

init_second = mexpv(delay,Q_pre',m.init_cond);
init_second =  init_second(m.Cind);

f = -m.totI(m.Cind)*alpha*(1-tau); % post admin, alpha is also reduced
part2 = Q_post(m.Cind,m.Cind)\f;

RR = part1(m.Iind) + init_second'*part2;


%If R* > 0 then calculate the early growth rate.
if RR > 0
    rr = fzero( @(r)return_val_const(r), [0, 6]);
else
    rr = 0;
end



    function [result] = return_val_const(r)
        % the function we want to minimize to find r
        % for the constant delay model
        
        part1_mat = Q_pre(m.Cind,m.Cind)-id*r;
        
        part1 = phiv(delay,part1_mat,u,v);
        
        init_second = mexpv(delay,Q_pre',m.init_cond);
        init_second =  init_second(m.Cind);
        
        part2_mat = Q_post(m.Cind,m.Cind)-id*r;
        
        f = -m.totI(m.Cind)*alpha*(1-tau); % post admin, alpha is also reduced
        part2 = part2_mat\f;
        
        result = (part1(m.Iind) + (init_second'*part2)*exp(-r*delay))*(2*sigma/(r+2*sigma))^2 -1;
        
    end


end
