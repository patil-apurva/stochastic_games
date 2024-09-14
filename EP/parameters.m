r = 0.1; % radius of the adversary

s2 = 0.1 %Sigma=s.I, s2=s^2
eta = 0.2 %lagrange multiplier
% gamma2 = 2 %Rv = gamma^2

z0 = [0.3; 0.3]; %initial position
runs = 100000; %MC runs
traj_num = 300; %number of trajectories to plot

h = 0.01; % time step
t0 = 0.0; %initial time
T = 2; % final time

b = 0; %V=b.||z(t)||^2
d = 0; %psi=d.||z(T)||^2
%==========================================================================

% lambda = s2*gamma2/(gamma2-1); %PDE linearization constant
s = sqrt(s2); %Sigma=s.I

xP = -0.5;
yP = -0.5;
xQ = 1; 
yQ = 1;

figure (2)
hold on;
set(gca, 'FontName', 'Arial', 'FontSize', 18)
xlabel('$z_x$', 'Interpreter','latex', 'FontSize', 30); ylabel('$z_y$', 'Interpreter','latex','FontSize', 30); 
xticks(xP:0.2:xQ)
yticks(yP:0.2:yQ)
set(gca,'LineWidth',1)
ax = gca;
ax.LineWidth = 1;
ax.Color = 'w';
axis equal;
xlim([xP, xQ]);
ylim([yP, yQ]);
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', r);
plot (p, 'FaceColor', 'r')