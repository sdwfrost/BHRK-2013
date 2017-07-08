function [RR,rr] = calculate_no_antivirals(N,alpha,beta,gamma,sigma)
% Calculate R* (RR) and the early growth rate (rr) for the SEEIIR model,
% assuming no interventions.


m = SEEIIR(N);


QQ1 = m.input_params(beta/(m.N-1),sigma,gamma,0);
Qc = QQ1(m.Cind,m.Cind);
id = speye(length(m.Cind));
f = -m.totI(m.Cind)*alpha;


% This calculates the household reproductive ratio R*
sol = Qc\f;
RR = sol(m.Eind);


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
        
        result = part2(m.Eind)-1;
        
    end


end