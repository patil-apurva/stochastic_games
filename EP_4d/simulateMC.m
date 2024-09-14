function S_tau_all = simulateMC(eps_t_all_1, eps_t_all_2, zt_prime_1, zt_prime_2, zt_prime_3, zt_prime_4, f_zt_prime_1, f_zt_prime_2, f_zt_prime_3, f_zt_prime_4, t, h, T, b, s, r, eta, d)
        
    eps_t_prime_1 = eps_t_all_1; %standard normal noise at t
    eps_t_prime_2 = eps_t_all_2; %standard normal noise at t
    
    S_tau = 0; %the cost-to-go of the state dependent cost of a sample path
    safe_flag_tau = 1;
    
    for t_prime = t:h:T-h % this loop is to compute S(tau_i)
                
        S_tau = S_tau + h*b*(zt_prime_1*zt_prime_1 + zt_prime_2*zt_prime_2); %add the state dependent running cost
        
        zt_prime_1 = zt_prime_1 + f_zt_prime_1*h; %move tau ahead
        zt_prime_2 = zt_prime_2 + f_zt_prime_2*h; %move tau ahead
        zt_prime_3 = zt_prime_3 + f_zt_prime_3*h + s*eps_t_prime_1*sqrt(h); %move tau ahead
        zt_prime_4 = zt_prime_4 + f_zt_prime_4*h + s*eps_t_prime_2*sqrt(h); %move tau ahead

        if ( (zt_prime_1*zt_prime_1 + zt_prime_2*zt_prime_2) <= (r*r) )%if yes means t_prime=t_exit
            S_tau = S_tau + eta; %add the boundary cost to S_tau
            safe_flag_tau = 0;
            break; %end this tau 
        end

        eps_t_prime_1 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        eps_t_prime_2 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        
        f_zt_prime_1 = zt_prime_3; %f_zt_prime_1 at new t_prime. Will be used in the next iteration 
        f_zt_prime_2 = zt_prime_4; %f_zt_prime_2 at new t_prime. Will be used in the next iteration 
        f_zt_prime_3 = 0; %f_zt_prime_3 at new t_prime. Will be used in the next iteration 
        f_zt_prime_4 = 0; %f_zt_prime_4 at new t_prime. Will be used in the next iteration 
        
    end

    if(safe_flag_tau==1) %if tau has not collided 
        S_tau = S_tau + d*(zt_prime_1*zt_prime_1 + zt_prime_2*zt_prime_2); %add the terminal cost to S_tau
    end

    S_tau_all = S_tau;