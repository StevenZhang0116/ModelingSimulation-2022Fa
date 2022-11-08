% Generating examples of HH model.
Inj=50;  % Injection current
tspan=30;  % Simulating time
V0= -70; n0=0.3; m0=0.5; h0=0.06;  % Initial values for V, n (K gating variable),
                                  % m, h (Na gating variables)
plot_flag=1;         % 1 for generating plot. 0 for not
close all
[V,~,~,~,~] = hhrun(Inj, tspan, V0, m0, h0,n0,plot_flag);


% To compare several simulations, run the following code:

% V0=-100; n0=0.6;
% hhrun(Inj, tspan, V0, m0, h0,n0,plot_flag);
% 
% V0=-30; n0=0.55;
% hhrun(Inj, tspan, V0, m0, h0,n0,plot_flag);
% 
% V0=-60; n0=0.5;
% hhrun(Inj, tspan, V0, m0, h0,n0,plot_flag);
% 
% V0=50; n0=0.4;
% hhrun(Inj, tspan, V0, m0, h0,n0,plot_flag);
% 
% V0=40; n0=0.75;
% hhrun(Inj, tspan, V0, m0, h0,n0,plot_flag);

