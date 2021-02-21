clear all; close all; clc

%% fluid parameters
Pr      = 1;        % Prandtl number

%% numerical parameters
tol     = 1e-8;     % error tolerance for Newton-Raphson loop
etamax  = 50;       % domain height
itermax = 100;      % maximum iterations during Newton-Raphson loop

%% solve for velocity distribution

% initial guess of boundary values for shooting method
d2f_old = 0.5;
d2f_new = 0.51;

[eta,soln1] = ode45(@ode_Bl_f,[0,etamax],[0,0,d2f_old]);
[eta,soln2] = ode45(@ode_Bl_f,[0,etamax],[0,0,d2f_new]);

% error of initial guess
err_old      = soln1(end,2)-1;
err_new      = soln2(end,2)-1;

% Newton-Raphson loop
k = 1; % initialize loop counter
while(abs(err_new)>tol && k<itermax)
    delta   = - err_new/((err_new-err_old)/(d2f_new-d2f_old));  % boundary correction
    d2f_old = d2f_new;                                          % update boundary values
    d2f_new = d2f_new + delta;                                  % update boundary values
    [eta,soln]= ode45(@ode_Bl_f,[0 etamax],[0 0 d2f_new]);      % solve ODE
    err_old = err_new;
    err_new = soln(end,2)-1;
    k       = k+1;
end
d2f0 = soln(1,3);

%% solve coupled momentum and energy equation to acquire temperature solution

% initial guess of boundary values for shooting method
d1g_old     = -0.5;
d1g_new     = -0.51;

% solve for initial solution before Newton-Raphson iteration
[eta,soln1] = ode45(@(eta,f) ode_Bl_fg( Pr,eta,f ),[0,etamax],[0,0,d2f0,1,d1g_old]);
[eta,soln2] = ode45(@(eta,f) ode_Bl_fg( Pr,eta,f ),[0,etamax],[0,0,d2f0,1,d1g_new]);

% error of initial guess
err_old     = soln1(end,4);
err_new     = soln2(end,4);

% Newton-Raphson loop
k = 1;    % initialize loop counter
while(abs(err_new)>tol && k<itermax)
    delta   = - err_new/((err_new-err_old)/(d1g_new-d1g_old));          % boundary correction
    d1g_old = d1g_new;                                                  % update boundary values
    d1g_new = d1g_new + delta;                                          % update boundary values
    IC      = [0,0,d2f0,1,d1g_new];
    [eta,soln] = ode45(@(eta,f) ode_Bl_fg( Pr,eta,f ),[0,etamax],IC);   % solve ODE system
    err_old = err_new;
    err_new = soln(end,4);
    k       = k+1;
end

% extract velocity and temperature distribution
U   = soln(:,2);
T   = soln(:,4);

%% load White's data
[ ~, theta ]    = white_theta();    % read temperature solution by White ()
[etaW,~,UW,~]   = white_f();        % read temperature solution by White ()

%% plot solution and compare to White's data
figure(1), clf; hold on
a(1)=plot(U,eta,'k-','LineWidth',2);
a(2)=plot(T,eta,'b-','LineWidth',2);
plot(theta(3).etatheta(:,2),theta(3).etatheta(:,1),'r--','LineWidth',2)
plot(UW,etaW,'r--','LineWidth',2)
xlim([0,1.05]); ylim([0,10])
set(gca,'FontSize',20)
ylabel('$\eta$','Interpreter','Latex','FontSize',30)
leg=legend('f','$\theta$');
set(leg,'FontSize',30,'Interpreter','Latex')
