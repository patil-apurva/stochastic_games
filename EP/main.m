clearvars
close all
run('parameters.m')

% GAMMA2 = [1.1, 1.2, 1.3, 1.5, 2];
% GAMMA2 = [1.4, 1.7, 2.2, 2.5];
GAMMA2 = 2;

for gamma2 = GAMMA2
    gamma2
    lambda = s2*gamma2/(gamma2-1); %PDE linearization constant
    
%     ETA = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1];
    ETA = 0.2;
    PFAIL = zeros(1, length(ETA));
    iter = 0;

    for eta = ETA
%         eta
        iter = iter+1;
        seed = 1256;
        rng( seed ); % Reset the CPU random number generator.
        gpurng( seed ); % Reset the GPU random number generator.

        fail_cnt = 0; %number of trajectories failed
        flag1 = 0;
        flag2 = 0;

        % tic

        for traj_itr = 1:traj_num

            traj_itr
            Z = []; %to store all positions of this trajectory
            Z = [Z, z0]; %stack the initial position

            zt = z0; %start the state from the given initial position
            f_zt = [0; 0]; %initial f_zt
            safe_flag_traj = 1;

            for t = t0:h:T-h % this loop is to find u(t) v(t) and z(t) at each time step t => z(t+h) = z(t) + f(z(t)).h + G_u.u(t).h + G_v.v(t).h + Sigma*dw
                eps_t_all_1 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_1(t) at the start of each sample path starting at time t and state zt
                eps_t_all_2 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_2(t) at the start of each sample path starting at time t and state zt

                S_tau_all = arrayfun(@simulateMC, eps_t_all_1, eps_t_all_2, zt(1), zt(2), f_zt(1), f_zt(2), t, h, T, b, s, r, eta, d); %an array that stores S(tau) of each sample path starting at time t and state zt

                eps_t_all_arr = gather([eps_t_all_1; eps_t_all_2]); %concatenate eps_t_all_arr_1 and eps_t_all_arr_2 in an array

                denom_i = exp(-S_tau_all/lambda); %(size: (1 X runs))
                numer = eps_t_all_arr*(denom_i.'); %(size: (2 X 1))
                denom = sum(denom_i); %scalar

                ut = (gamma2*s*numer)/((gamma2-1)*sqrt(h)*denom); %the agent control input
                vt = (s*numer)/((gamma2-1)*sqrt(h)*denom); %the adversary control input

        %         if(any(isnan(ut(:))) || any(isnan(vt(:))))
        %             fprintf("error!")
        %             return
        %         end

                %move the trajectory forward
                eps = randn(2,1);
                zt = zt + f_zt*h + ut*h - vt*h + s*eps*sqrt(h); %update the position with the control input ut=> z(t+h) = z(t) + f.h + G_u.u(t).h + G_v.v(t).h + sigma*dw
                Z = [Z, zt]; %stack the new position

                if( (zt.'*zt) <= (r*r) ) %if yes means trajectory has crossed the safe set
                    fail_cnt = fail_cnt+1; 
                    safe_flag_traj = 0;
                    break;  %end this traj    
                end 

                f_zt = [0; 0]; %update f(z(t)) for the next t => t=t+h. Will be used in the next iteration 
            end

               if (safe_flag_traj==1)
        %             plot (Z(1, :), Z(2, :), 'g', 'LineWidth',1);
                    Z1 = Z;
                    flag1 = 1;
               else
        %             plot (Z(1, :), Z(2, :), 'b', 'LineWidth',1);
                    Z2 = Z;
                    flag2 = 1;
               end
    
                 if((flag1 && flag2) == 1)
                     plot (Z1(1, :), Z1(2, :), 'g', 'LineWidth',1);
                     plot (Z2(1, :), Z2(2, :), 'b', 'LineWidth',1);
                     plot (Z1(1, end), Z1(2, end), 'o', 'MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor','g') 
                     plot (Z2(1, end), Z2(2, end), 'o', 'MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor','b')
                     break
                 end
        end

        fail_prob = fail_cnt/traj_num;
        PFAIL(iter) = fail_prob;
    end
    PFAIL
end
% toc
 
plot(z0(1), z0(2), 'p', 'MarkerSize', 20, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y')
figname = ['new_eta=',num2str(eta),'.fig'];
saveas(gcf,figname)

save trajectories_rel.mat Z1 Z2 