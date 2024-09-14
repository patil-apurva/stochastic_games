clearvars
close all
run('parameters.m')

seed = 1234;
rng( seed ); % Reset the CPU random number generator.
gpurng( seed ); % Reset the GPU random number generator.

fail_cnt = 0; %number of trajectories failed

% J0_avg = 0; %average of all J0_traj
% time_steps = (T-t0)/h + 1;
% J0_avg_with_time = zeros(time_steps, 1); %average of all J0_traj with time 

for traj_itr = 1:traj_num
    traj_itr
    
%     J0_traj = 0; %J0 of this trajectory
%     J0_traj_with_time = zeros(time_steps, 1);%J0 of this trajectory with time
    
    X = []; %to store all positions of this trajectory
    X = [X, x0]; %stack the initial position

    xt = x0; %start the state from the given initial position
    % dynamics model is from the paper: https://arxiv.org/pdf/2009.14775.pdf
    f_xt = [k1*xt(1) + xt(3)*cos(xt(4)); k1*xt(2) + xt(3)*sin(xt(4)); k2*xt(3); k2*xt(4)]; %initial f_xt
    safe_flag_traj = 1;
    
%     iter = 1; %for the index of the time step
    for t = t0:h:T-h % this loop is to find u(t), v(t) and x(t) at each time step t => x(t+h) = x(t) + f(x(t)).h + G_u.u(t).h + G_v.v(t).h + Sigma*dw
      
        % Note that process noises are only applied to the last two states,
        % which are directly actuated
        eps_t_all_1 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_1(t) at the start of each sample path starting at time t and state xt
        eps_t_all_2 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_2(t) at the start of each sample path starting at time t and state xt

        S_tau_all = arrayfun(@simulateMC, eps_t_all_1, eps_t_all_2, xt(1), xt(2), xt(3), xt(4), f_xt(1), f_xt(2), f_xt(3), f_xt(4), t, h, T, b, s, xR1, xS1, yR1, yS1, xR2, xS2, yR2, yS2, xP, xQ, yP, yQ, eta, d, k1, k2); %an array that stores S(tau) of each sample path starting at time t and state xt

        eps_t_all_arr = gather([eps_t_all_1; eps_t_all_2]); %concatenate eps_t_all_arr_1 and eps_t_all_arr_2 in an array

        denom_i = exp(-S_tau_all/lambda); %(size: (1 X runs))
        numer = eps_t_all_arr*(denom_i.'); %(size: (2 X 1))
        denom = sum(denom_i); %scalar

        ut = (gamma2/(gamma2-1)) * (s/sqrt(h)) * (numer/denom); %the agent control input
        vt = -(1/(gamma2-1)) * (s/sqrt(h)) * (numer/denom); %the adversary control input

    %             if((any(isnan(ut(:)))) ||(any(isnan(vt(:)))))
    %                 fprintf("error!")
    %                 return
    %             end

%         J0_traj = J0_traj + (b*(xt.')*xt + 0.5*(ut.')*ut)*h; %add the running cost 
%         J0_traj_with_time(iter) = J0_traj_with_time(iter) + gather((b*(xt.')*xt + 0.5*(ut.')*ut)*h); %add the cost at time step iter
        
        %move the trajectory forward
        eps = randn(2,1);
        %update the position with the control inputs ut and vt=> x(t+h) = x(t) + f.h + G_u.u(t).h + G_v.v(t).h + Sigma*dw
        xt = xt + f_xt*h + G_u*(ut*h + vt*h + s*eps*sqrt(h));

        X = [X, xt]; %stack the new position

        if(((xt(1)>=xR1) && (xt(1)<=xS1) && (xt(2)>=yR1) && (xt(2)<=yS1)) || ((xt(1)>=xR2) && (xt(1)<=xS2) && (xt(2)>=yR2) && (xt(2)<=yS2)) || ((xt(1)<=xP) || (xt(1)>=xQ) || (xt(2)<=yP) || (xt(2)>=yQ))) %if yes means trajectory has crossed the safe set
%             J0_traj = J0_traj + eta; %add the boundary cost
%             J0_traj_with_time(iter) = J0_traj_with_time(iter) + eta; %add the cost at time step iter
            
            fail_cnt = fail_cnt+1; 
            safe_flag_traj = 0;
            break;  %end this traj    
        end 

        f_xt = [k1*xt(1) + xt(3)*cos(xt(4)); k1*xt(2) + xt(3)*sin(xt(4)); k2*xt(3); k2*xt(4)]; %update f(x(t)) for the next t => t=t+h. Will be used in the next iteration 
        
%         iter = iter + 1; %update the time step index
    end

     if (safe_flag_traj==1)
        plot (X(1, :), X(2, :), 'g', 'LineWidth',1);
    else
        plot (X(1, :), X(2, :), 'b', 'LineWidth',1);
     end
     
%      if(safe_flag_traj==1) %if trajectory has not collided 
%         J0_traj = J0_traj + d*(xt.')*xt; %add the terminal cost
%         J0_traj_with_time(iter) = J0_traj_with_time(iter) +  gather(d*(xt.')*xt); %add the cost at time step iter
%      end
%     
%      J0_avg = (J0_avg*(traj_itr-1) + J0_traj)/traj_itr;
%      J0_avg_with_time = (J0_avg_with_time*(traj_itr-1) + J0_traj_with_time)/traj_itr;
end

% sum(J0_avg_with_time)
% J0_avg
fail_prob = fail_cnt/traj_num

plot(x0(1),x0(2),'p', 'MarkerSize',40,'MarkerEdgeColor','k','MarkerFaceColor','y')
plot(0,0,'p', 'MarkerSize',40,'MarkerEdgeColor','k','MarkerFaceColor','m')
figname = ['u_safe_eta=',num2str(eta),'_gamma2=',num2str(gamma2),'_s2=',num2str(s2),'.fig'];
saveas(gcf,figname)

% J0_avg_with_time_safe = J0_avg_with_time;
% save timewise_cost_safe.mat J0_avg_with_time_safe