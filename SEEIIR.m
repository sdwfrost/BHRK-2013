classdef SEEIIR
    % class which combines all code for creating the SEEIIR matrix template and
    % inputing parameters in an efficient manner.
    
    % Usually in computations the template only needs to be built once (for a
    % given value of N), but differnt parameters will be input to it. So code
    % for creating the template is inefficient, but this is only called once.
    
    properties
        
        sps;
        N;
        totI; totE; S;
        init_cond;
        Cind; Iind; Eind
        
    end
    
    
    properties (Access = private)
        
        react_latent;
        react_binf;
        react_infp
        react_rec
        react_inf_I;
        react_inf_S;
        
        rows; cols;
        
    end
    
    
    methods
        
        function m = SEEIIR(N)
            % constructor. Build a template of size N.
            
            m.sps = round(1/120*(N+1)*(N+2)*(N+3)*(N+4)*(N+5));
            m.N = N;
            
            [m.S,E1,E2,I1,I2,R] = m.state_vectors(N);
            
            m.totI=(I1+I2);
            m.totE = (E1+E2);
            
            m.rows = [find(E1>0);find(E2>0);find(I1>0);find(I2>0);find(m.S>0); (1:m.sps)'];
            
            m.cols = [find(E2>0);find(I1>0);find(I2>0);find(R>0);find(E1>0); (1:m.sps)'];
            
            
            % this is their order in the rows and cols vectors.
            m.react_latent = E1(E1>0);
            m.react_binf = E2(E2>0);
            m.react_infp = I1(I1>0);
            m.react_rec = I2(I2>0);
            m.react_inf_I = m.totI(E1>0);
            m.react_inf_S = m.S(m.S>0);
            
            m.Cind = find(not((m.totI) == 0 & (m.totE) == 0));
            
            % I1=1 initial condition.
            poss = find(I1 == 1 & m.S == (N-1));
            m.Iind = m.Cind==poss(1);
            
            % E1=1 initial condition.
            poss = find(E1 == 1 & m.S == (N-1));
            m.Eind = m.Cind==poss(1);
            
            % initcond part for use with the constant delay algorthm.
            m.init_cond = (I1==1 & m.S==(N-1));
            m.init_cond = +m.init_cond;
            
        end
        
        
        function [g2] = input_params(m,beta,gamma,sigma,eps)
            % helper function to fill out the sparse matrix.
            
            ratesA = [m.react_latent*2*sigma; m.react_binf*2*sigma;...
                m.react_infp*2*gamma; m.react_rec*2*gamma;...
                beta*(m.react_inf_I+eps).*m.react_inf_S;...
                -(2*sigma*m.totE+2*gamma*m.totI+beta*(m.totI+eps).*m.S)];
            
            
            g2 = sparse(m.rows,m.cols,ratesA);
            
        end
        
        
        function result = Rstar(m,alpha,beta,gamma,sigma,eps)
            % calculate R*, not to be used for heavy lifting.
            
            QQ1 = m.input_params(beta/(m.N-1),sigma,gamma,eps);
            
            Qc = QQ1(m.Cind,m.Cind);
            
            f = -m.totI(m.Cind)*alpha;
            
            sol = Qc\f;
            
            result = sol(m.Iind);
            
        end
        
        
    end
    
    
    methods(Static)
        
        function [S E1 E2 I1 I2 R] = state_vectors(N)
            % Returns the state vectors for the SEEIIR model for use in building
            % matricies etc..
            
            DIM = round(1/120*(N+1)*(N+2)*(N+3)*(N+4)*(N+5));
            
            S = zeros(DIM,1);
            E1 = zeros(DIM,1);
            E2 = zeros(DIM,1);
            I1 = zeros(DIM,1);
            I2 = zeros(DIM,1);
            R = zeros(DIM,1);
            
            for ss=0:N
                for e1 = 0:(N-ss)
                    for e2 = 0:(N-ss-e1)
                        for i1 = 0:(N-ss-e1-e2)
                            for i2 = 0:(N-ss-e1-e2-i1)
                                
                                
                                S(tfmSIR4(ss,e1,e2,i1,i2,N)) = ss;
                                E1(tfmSIR4(ss,e1,e2,i1,i2,N)) = e1;
                                E2(tfmSIR4(ss,e1,e2,i1,i2,N)) = e2;
                                I1(tfmSIR4(ss,e1,e2,i1,i2,N)) = i1;
                                I2(tfmSIR4(ss,e1,e2,i1,i2,N)) = i2;
                                R(tfmSIR4(ss,e1,e2,i1,i2,N)) = N-ss-e1-e2-i1-i2;
                                
                            end
                        end
                    end
                end
            end
            
            
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
    
    
end