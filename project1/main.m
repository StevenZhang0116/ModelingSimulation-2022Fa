clear
close all

format long

% load function registrator
ca_res = form_curr();

% v_in = -80 * 10^-3;      % voltage inside the cell [-70mV]

rest_vol = -80; % first level of voltage
acti_vol = 20; % second level of voltage

% define function of voltage inside the cell
steady_calsite = []; 
steady_open = [];
single_current_ca = [];

kk = 80:0.1:80;
for tttt = 1:length(kk)
    volt = kk(tttt);
    volt

    ptbound = [0,0.01,0.02,0.03];
    vbound = 10^(-3) * [rest_vol,volt,rest_vol];
    
    % v_in_func = @(t) 10^(-3) * ...
    %     (rest_vol .* (t >= 0) .* (t <= 0.01) ...
    %     + acti_vol .* (t > 0.01) .* (t <= 0.02) ...
    %     + -80 .* (t > 0.02) .* (t <= endtime));
    v_in_func = @(t) 0; 
    for ii = 1:length(vbound)
        v_in_func = @(t) (v_in_func(t) + vbound(ii) .* (t >= ptbound(ii)) .* (t < ptbound(ii+1)) );
    end
    
    endtime = ptbound(end); % end of simulation time
    
    
    gridbreak = 1e-5;
    timegrid = 0:gridbreak:endtime;
    v_in_set = v_in_func(timegrid);
    
    figure(1)
    plot(timegrid(1:end-1)*1000,v_in_set(1:end-1)*1000,'LineStyle','-','LineWidth',2);
    xlabel('Time (ms)','FontSize',12)
    ylabel('Voltage inside the cell (mV)','FontSize',12)
    xlim([0 endtime*1000])
    title('Voltage inside the cell vs. Time')
    
    v_out = 0 * 10^-3;       % voltage outside the cell [0mV]
    
    global nn ca2in ca2out alpha_v gamma_v_r alpha_s_0 beta_s a_0 r_0 l_0
    
    time_coll = []; % time vector
    state_coll = []; % state param vector
    fca_coll = []; % fca vector
    casite_coll = []; % ca_site vector
    
    v_in_set(end) = [];
    v_check = v_in_set;

    %%
    Inj=200;  % Injection current
    tspan=30;  % Simulating time
    V0= -70; n0=0.3; m0=0.5; h0=0.06;  % Initial values for V, n (K gating variable),
                                      % m, h (Na gating variables)
    plot_flag=1;         % 1 for generating plot. 0 for not
    close all
    [V,~,~,~,t] = hhrun(Inj, tspan, V0, m0, h0,n0,plot_flag);
    v_in_set = V/1000;
    %%
    
    for ind = 1:length(v_in_set)
        v_in = v_in_set(ind);
    
        ca2in = 0;               % [Ca++]_in, make approximation, P25, [~10-7 M]
        ca2out = 4 * 1e-3;       % [Ca++]_out, concentration of outside cell [~10^-3 M]
        a_0_r = 0.5 * 1e-9;      % cross sectional radius;
        a_0 = pi * (a_0_r).^2;    % cross sectional area of channel
        r_0 = 10 * 1e-6;          % distance from the open channel
        l_0 = 10 * 1e-9;         % length of specific channel [m]
        
        [nn,fca,casite] = ca_res.f_cal_ca(ca2in,ca2out,v_in,v_out,r_0,a_0,l_0);
        fca_coll(end+1) = fca;
        casite_coll(end+1) = casite;
        
        % all units in notes are given as ms-1, but to standardize it, we need to make it as
        % s-1
        gamma_v_r = 10 * 1e3;           % constant, P38, [s-1]
        alpha_v = 10 * 1e3;             % constant, P38, [s-1]
        alpha_s_0 = 20 * 1e3;           % constant, P26, [s-1]
        beta_s = 10 * 1e3;              % constant, P26, [s-1]
        
        t_0 = 0;
        t_end = gridbreak;
        
        if ind == 1
            f0 = zeros(1,12);
            f0(1) = 1;
        end
        
        opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-6);
        [t,statedyn] = ode45(ca_res.cal_complete_state,[t_0 t_end],f0,opts);
        
        t = t + (ind-1)*gridbreak;
    
        time_coll = [time_coll;t];
        state_coll = [state_coll;statedyn];
    
        f0 = statedyn(end,:); % update for the initial condition for the next round
    
    end
    
    figure(2)
    hold on
    plot(time_coll*1000,state_coll(:,1),'LineWidth',2);
    % plot(time_coll,state_coll(:,2),'LineWidth',2);
    % plot(time_coll,state_coll(:,3),'LineWidth',2);
    % plot(time_coll,state_coll(:,4),'LineWidth',2);
    % plot(time_coll,state_coll(:,5),'LineWidth',2);
    plot(time_coll*1000,state_coll(:,6),'LineWidth',2);
    % plot(time_coll,state_coll(:,7),'LineWidth',2);
    % plot(time_coll,state_coll(:,8),'LineWidth',2);
    % plot(time_coll,state_coll(:,9),'LineWidth',2);
    % plot(time_coll,state_coll(:,10),'LineWidth',2);
    % plot(time_coll,state_coll(:,11),'LineWidth',2);
    plot(time_coll*1000,state_coll(:,12),'LineWidth',2);
    legend({'X_0','X_5','Y_5'},'FontSize',12)
    
    xlabel('Time (ms)','FontSize',12)
    ylabel('State Val (dimensionless)','FontSize',12)
    xlim([0 endtime*1000])
    title('State variable vs. Time','FontSize',12)
    
    % opening rate
    openr = state_coll(:,6);
    steady_open(end+1) = openr(end);

    % single channel ca2+ current
    single_current_ca(end+1) = fca_coll(end);
    
    % 
    tr = timegrid(1:end-1);
    rtr = round(tr,5);
    casite_val = [];
    for ii = 1:length(ptbound)-1
        meant = 1/2*(ptbound(ii)+ptbound(ii+1));
        hh = find(rtr == meant);
        casite_val(end+1) = casite_coll(hh-1);
    end
    
    figure(3)
    plot(tr,casite_coll)
    
    casite_func = @(t) (0);
    for ii = 1:length(casite_val)
        casite_func = @(t) (casite_func(t) + casite_val(ii) .* (t >= ptbound(ii)) .* (t < ptbound(ii+1)));
    end
    
    calsite = casite_func(time_coll);
    calsite_total = openr .* calsite;
    
    figure(4)
    plot(time_coll,calsite_total)

    steady_calsite(end+1) = calsite_total(end-1);

    pause
    close all

end


figure()
plot(kk,steady_calsite,'LineWidth',2)
xlabel('Voltage (mV)','FontSize',12)
ylabel('[Ca^{2+}]_{site} (M)','FontSize',12)

figure()
plot(kk,steady_open,'LineWidth',2)
xlabel('Voltage (mV)','FontSize',12)
ylabel('Opening Rate','FontSize',12)

figure()
plot(kk,single_current_ca,'LineWidth',2)
xlabel('Voltage (mV)','FontSize',12)
ylabel('Single Channel Current','FontSize',12)


