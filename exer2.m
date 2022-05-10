clear
clc


%% Parameter Set %%
poles = [1, 1; 2, 2; 5, 5; 10, 10; 50, 50; 100, 100; 150, 150; 200, 200; 250, 250; 300, 300];
% poles = [1, 5; 4, 7; 15, 30; 40, 20; 60 35; 20, 100; 100, 120; 150, 160; 250, 150; 350, 400]; 

t = (0:1e-5:5)';

tspan = (0:1e-5:5)';
%% Function Generation %%

u1 = @(t) 3*sin(2*t);
u1dot = @(t) 6*cos(2*t);
u1ddot = @(t) -12*sin(2*t);
u2 = @(t) 2;


u1s = u1((tspan(:)));
u2s=ones(length(tspan), 1)*2;

% theta = [1/(rhat*chat) - 2*rho;1/(lhat*chat) - rho*rho; 1/(rhat*chat); 0; 1/(rhat*chat); 1/(lhat*chat)];


[Vr, Vc] = v(tspan);


%% Wrong Values comment this out %%
% Vc(15000) = Vc(15000) + 100 * Vc(15000);
% Vc(26000) = Vc(26000) + 200 * Vc(26000);
% Vc(35000) = Vc(35000) + 150 * Vc(35000);
% Vr(15000) = Vr(15000) + 100 * Vr(15000);
% Vr(26000) = Vr(26000) + 200 * Vr(26000);
% Vr(35000) = Vr(35000) + 150 * Vr(35000);
%% Plotting Vr and Vc
figure;

subplot(2, 1,1);
plot(t, Vr);
title('Voltage drop across the resistor', 'interpreter', 'latex');
ax = gca;
ax.TitleFontSizeMultiplier = 4;
ylabel('$V_R$', 'interpreter', 'latex', 'FontSize', 20);
xlabel('$t$', 'interpreter', 'latex', 'FontSize', 20);
legend('$V_R$', 'interpreter', 'latex');

hold on;

subplot(2,1,2);
plot(t,Vc);
title('Voltage drop across the capacitor', 'interpreter', 'latex');
ax = gca;
ax.TitleFontSizeMultiplier = 4;
ylabel('$V_C$', 'interpreter', 'latex', 'FontSize', 20);
xlabel('$t$', 'interpreter', 'latex', 'FontSize', 20);
legend('$V_C$', 'interpreter', 'latex');
legend('Vc');

% Vr = Vr.';
% Vc = Vc.';




for i=1:length(poles)
    tic
    system = generate_system(poles(i,:));
    
    zeta(:,1) = lsim(system(1),Vc, tspan);
    zeta(:,2) = lsim(system(2), Vc, tspan);
    zeta(:,3) = lsim(-system(1),u1s, tspan);
    zeta(:,4) = lsim(-system(2),u1s, tspan);
    zeta(:,5) = lsim(-system(1),u2s, tspan);
    zeta(:,6) = lsim(-system(2),u2s, tspan);
    
    
    zetas(i).z = zeta;
    zetas(i).p = poles(i,:);
    zetas(i).l = system;
    
    theta_0 = optimal_theta(Vc, zeta);
    optimal = theta_0 + [poles(i,1)+poles(i,2); poles(i,1)*poles(i,2); zeros(4,1)];
    rc = 1/(optimal(1,:));
    lc = 1/(optimal(2,:));
    
    zetas(i).theta = theta_0;
    zetas(i).rc = rc;
    zetas(i).lc = lc;
    
    
    dvcdt = @(t, vc) [vc(2); (u2(t)/lc) + (u1dot(t))/rc - (1/lc)*vc(1) - (1/rc)*vc(2)];
    dvrdt = @(t, v) [v(2); -v(2)/rc - v(1)/lc + u1(t)/lc + u1ddot(t)];
    [t, Vrhat] = ode45(dvrdt, tspan,[v(0);0]);
    [t, Vchat] = ode45(dvcdt, tspan, [0;v(0)]);
    zetas(i).vrhat = Vrhat(:,1);
    zetas(i).vchat = Vchat(:,1);
    zetas(i).time = toc;
    
    %% Plotting
    subplot(length(poles)/2,2,i)
    plot(t, Vc(:) - zetas(i).vchat(:,1), '-r');
%     title( sprintf('Error of approximation using the Least Squares Method for filter $(s+ %.1f)(s+%.1f)$', poles(i,1),poles(i,2)),'Interpreter', 'latex');
    title( sprintf('Error of approximation using the Least Squares Method for filter $(s+ %d)^2$', poles(i,1)), 'Interpreter', 'latex');
    ax = gca;
    ax.TitleFontSizeMultiplier = 1.5;
    ylabel('$e=V_C- \hat{V}_C$',  'interpreter', 'latex', 'FontSize', 15);
    xlabel('$t$',  'interpreter', 'latex', 'FontSize', 15);
    legend('$e=V_C - \hat{V_C}$',  'interpreter', 'latex');

    
   
end
figure;
for i =1:length(zetas)
    subplot(length(poles)/2,2,i)
    plot(t, Vr(:) - zetas(i).vrhat(:,1), '-b');
     title( sprintf('Error of approximation using the Least Squares Method for filter $(s+ %d)^2$', poles(i,1)), 'Interpreter', 'latex');
%     title( sprintf('Error of approximation using the Least Squares Method for filter $(s+ %.1f)(s+%.1f)$', poles(i,1),poles(i,2)),'Interpreter', 'latex');
    ylabel('$e=V_R- \hat{V}_R$',  'interpreter', 'latex', 'FontSize', 15);
    xlabel('$t$',  'interpreter', 'latex', 'FontSize', 15);
    legend('$e=V_R - \hat{V_R}$',  'interpreter', 'latex');
end