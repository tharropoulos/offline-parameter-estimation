clear
clc

%% Parameter set %%
poles = [0.1, 0.4; 0.5, 0.1; 0.2, 0.6; 0.8, 0.4; 0.9, 1; 2, 3; 2, 4; 5, 4; 5, 1;  6, 9]; 
% poles = [1:10; 1:10]';

tspan = 0:0.1:10;

syms msym bsym ksym;

m = 10;
b = 0.3;
k = 1.5;

y0 = [0; 0];

%% Function Generation %%
u = @(t) 10*sin(3*t) +5;
dydt = @(t,y) [y(2); (1/m)*(u(t)) - (b/m)*y(2) - (k/m)*y(1)];

us = u(tspan(:));


[t,y] = ode45(dydt, tspan, y0);

%% Y plot %%
figure;
plot(t, y(:,1), '-b');
hold on;
title('Displacement of mass $m$','interpreter','latex')
ax = gca;
ax.TitleFontSizeMultiplier = 4;
hold on;



figure;


%% Calculations %%
for i=1:length(poles)
    tic
    system = generate_system(poles(i,:));
    
    zeta(:,1) = lsim(system(1),y(:,1), tspan);
    zeta(:,2) = lsim(system(2), y(:,1), tspan);
    zeta(:,3) = lsim(system(3), us, tspan);
    
    theta = [bsym / msym - (poles(i,1)+poles(i,2)); ksym / msym - poles(i,1)*poles(i,2); 1/msym];
    
    theta_0 = optimal_theta(y,zeta);
    

    eqn = theta == optimal_theta(y,zeta);
    S = solve(eqn, [msym, bsym, ksym]);

    zetas(i).z = zeta;
    zetas(i).p = poles(i,:);
    zetas(i).l = system;
    zetas(i).theta = theta_0;
    zetas(i).mhat = double(S.msym);
    zetas(i).bhat = double(S.bsym);
    zetas(i).khat = double(S.ksym);
    
    mhat = double(S.msym);
    bhat = double(S.bsym);
    khat = double(S.ksym);
    
    dydt2 = @(t,y) [y(2); (1/mhat)*(u(t)) - (bhat/mhat)*y(2) - (khat/mhat)*y(1)];

    [t,yhat] = ode45(dydt2, tspan, y0);
    
    zetas(i).yhat = yhat(:,1);
    zetas(i).time = toc;
 
    %% Plotting
    subplot(5,2,i);
    plot(t,yhat(:,1) - y(:,1), '-r');
    title( sprintf('Error of approximation using the Least Squares Method for filter $(s+ %.1f)(s+%.1f)$', poles(i,1),poles(i,2)),'Interpreter', 'latex');
%     title( sprintf('Error of approximation using the Least Squares Method for filter $(s+ %d)^2$', poles(i,1)), 'Interpreter', 'latex');
    ax = gca;
    ax.TitleFontSizeMultiplier = 1.5;
    ylabel('$e=y - \hat{y}$',  'interpreter', 'latex', 'FontSize', 15);
    xlabel('$t$',  'interpreter', 'latex', 'FontSize', 15);
    legend('$e=y - \hat{y}$',  'interpreter', 'latex');

    
%     temp_var = strcat('zeta_', num2str(i));
%     eval(sprintf('%s = %g', temp_var, zeta));
end
