classdef SEEIIR_exp
    % class which combines all code for creating the SEEIIR matrix for the
    % model with an exponentially distributed delay.
    
    properties
        
        N;
        sps;
        Cind; Iind;
        I_vector;
        
    end
    
    
    properties (Access = private)
        
        react_latent;
        react_binf;
        react_infp
        react_rec
        react_inf_I;
        react_inf_S;
        react_ra;
        react_ra_diag
        
        rows; cols;
        totI; totE; S;
        
    end
    
    
    methods
        
        function m = SEEIIR_exp(N)
            % constructor. Build a template of size N.
            
            m.sps = round(1/120*(N+1)*(N+2)*(N+3)*(N+4)*(N+5));
            m.N = N;
            
            m.S = zeros(m.sps,1);
            E1 = zeros(m.sps,1);
            E2 = zeros(m.sps,1);
            I1 = zeros(m.sps,1);
            I2 = zeros(m.sps,1);
            R = zeros(m.sps,1);
            
            for ss=0:N
                for e1 = 0:(N-ss)
                    for e2 = 0:(N-ss-e1)
                        for i1 = 0:(N-ss-e1-e2)
                            for i2 = 0:(N-ss-e1-e2-i1)
                                
                                m.S(m.tfmSIR4(ss,e1,e2,i1,i2,N)) = ss;
                                E1(m.tfmSIR4(ss,e1,e2,i1,i2,N)) = e1;
                                E2(m.tfmSIR4(ss,e1,e2,i1,i2,N)) = e2;
                                I1(m.tfmSIR4(ss,e1,e2,i1,i2,N)) = i1;
                                I2(m.tfmSIR4(ss,e1,e2,i1,i2,N)) = i2;
                                R(m.tfmSIR4(ss,e1,e2,i1,i2,N)) = N-ss-e1-e2-i1-i2;
                                
                            end
                        end
                    end
                end
            end
            
            m.totI=(I1+I2);
            
            m.totE = (E1+E2);
            
            m.I_vector = [m.totI; m.totI];
            
            m.Cind = find(not(m.I_vector == 0 & [m.totE; m.totE] == 0));
            
            
            %% THIS VERSION STARTS FROM E1=1
            poss = find([E1, E1] == 1 & [m.S, m.S] == (N-1));
            m.Iind = m.Cind==poss(1);
            
            
            rows_E1 = find(E1>0);
            cols_E1 = find(E2>0);
            
            rows_E2 = find(E2>0);
            cols_E2 = find(I1>0);
            
            rows_I1 = find(I1>0);
            cols_I1 = find(I2>0);
            
            rows_I2 = find(I2>0);
            cols_I2 = find(R>0);
            
            rows_inf = find(m.S>0);
            cols_inf = find(E1>0);
            
            % Now need to find the transitions in rows_ra which are not allowed
            % and delete them.
            
            del1 = find(E1==1 & m.S==(N-1));
            
            del2 = find(E2==1 & m.S==(N-1));
            
            rows_ra = [1:del1-1 del1+1:del2-1 del2+1:m.sps]';
            
            cols_ra = rows_ra+m.sps;
            
            m.react_ra_diag = [ones(del1-1,1); 0; ones(del2-del1-1,1); 0; ones(m.sps-del2,1)];
            
            
            
            m.rows = [rows_E1 ; (rows_E1+m.sps); rows_E2 ;(rows_E2+m.sps);...
                rows_I1; (rows_I1+m.sps);...
                rows_I2; (rows_I2+m.sps); rows_inf; (rows_inf+m.sps);...
                rows_ra; (1:2*m.sps)'];
            
            m.cols = [cols_E1 ; (cols_E1+m.sps); cols_E2 ;(cols_E2+m.sps);...
                cols_I1; (cols_I1+m.sps);...
                cols_I2; (cols_I2+m.sps); cols_inf; (cols_inf+m.sps);...
                cols_ra; (1:2*m.sps)'];
            
            m.react_latent = [E1(E1>0);E1(E1>0)];
            m.react_binf = [E2(E2>0);E2(E2>0)];
            m.react_infp = [I1(I1>0);I1(I1>0)];
            m.react_rec = [I2(I2>0);I2(I2>0)];
            m.react_inf_I = m.totI(E1>0);
            m.react_inf_S = m.S(m.S>0);
            
            m.react_ra = ones(length(rows_ra),1);
            
        end
        
        
        function Q = input_params(m,beta,sigma,gamma,eps,tau,rho,ra)
            % input values into the matrix template.
            % the norm factor now needs to go into the parameters, it is no longer handled here.
            
            rates = [m.react_latent*2*sigma;...
                m.react_binf*2*sigma;...
                m.react_infp*2*gamma;...
                m.react_rec*2*gamma;...
                beta*(m.react_inf_I+eps).*m.react_inf_S;...
                beta*(1-tau)*(1-rho)*(m.react_inf_I+eps).*m.react_inf_S;...
                m.react_ra*ra;...
                -(2*sigma*m.totE+2*gamma*m.totI+beta*(m.totI+eps).*m.S+m.react_ra_diag*ra);...
                -(2*sigma*m.totE+2*gamma*m.totI+beta*(1-tau)*(1-rho)*(m.totI+eps).*m.S)];
            
            
            Q = sparse(m.rows,m.cols,rates);
            
        end
        
        function result = Rstar(m,alpha,beta,sigma,gamma,eps,tau,rho,ra)
            % return R*. Again not suitable for heavy lifting.
            
            QQ1 = m.input_params(beta/(m.N-1),sigma,gamma,eps,tau,rho,ra);
            
            Qc = QQ1(m.Cind,m.Cind);
            
            
            alpha_vector = [alpha*ones(1,length(m.Cind)/2) alpha*(1-tau)*ones(1,length(m.Cind)/2)];
            
            f = -m.I_vector(m.Cind).*alpha_vector';
            
            sol=Qc\f;
            
            result = sol(m.Iind);
            
        end
        
        
    end
    
    
    methods(Static)
        
        
        function y = tfmSIR4(m,n,p,q,r,N)
            % Convert five integers between 0 and N into a unique integer
            
            y  = (N-p-q-r)*n - 1/2*n*(n-3) + m + 1 + (1/6)*p*(11 - 6*p...
                + p^2 + 12*(N-q-r) - 3*p*(N-q-r) + 3*(N-q-r)^2)...
                - (1/24)*q*(-5 + q - 2*(N-r))*(q^2 - q*(5 + 2*(N-r))...
                + 2*(5 + (N-r)*(5 + (N-r))))...
                + (1/120)*r*(274 + r^4 - 5*r^3*(3 + N) + 5*N*(6 + N)*(15 + N*(6 + N))...
                - 5*r*(3 + N)*(15 + 2*N*(6 + N)) + 5*r^2*(17 + 2*N*(6 + N)));
            
            y = int64(y);
            
        end
        
        
    end
    
    
end