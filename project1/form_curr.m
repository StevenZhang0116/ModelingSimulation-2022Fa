function ca_res = form_curr()
    ca_res.f_cal_ca = @f_cal_res;
    ca_res.cal_complete_state = @cal_complete_state;

end


function [nnval,f_ca,ca_site] = f_cal_res(ca_in,ca_out,v_in,v_out,r_0,a_0,l_0)
    % ca_in:   Ca2+ concentration inside cell                       [mol/L]
    % ca_out:  Ca2+ concentration outside cell                      [mol/L]
    % v_in:    voltage inside the cell                              [V]
    % v_out:   voltage outside the cell                             [V]
    % 
    % f_ca:    Flux of ca2+ for individual channel            
    % nnval:   Normalization constant = qv/kT
    % ca_site: [Ca2+]_site, the concentration at a distance r_0 for the open channel
    %
    % f_ca if the channel is open, otherwise 0
    % ca_site if the channel is open, otherwise 0

    % constants registeration
    q = 1.602176*10^(-19); % elementary charge                     [C]
    a = a_0;                % cross sectional area of channel      [m2] 
    l = l_0;                % channel length                       [m]
    D = 13.12*10^(-10);     % diffusion coefficient of Ca2+        [m2 s-1]
    k = 1.380649*10^(-23);  % Boltamann's constant                 [m2 kg s-2 K-1]
    T = 273 + 36;           % Absolute temperature                 [K]
    r = r_0;                % Distance from the open channel       [m]

    v = v_in - v_out;       % membrane potential (inside-outside)  [V]

    p1 = k*T/q;             % compare ratio: ~25mV, T mostly fixed                                                                                                                                                                                         
     
    % normalization value
    nnval = v/p1;

    % P27
    sc_per = a*D/l;         % single-channel Ca2+ permeability

    % P25
    f_ca = sc_per*ca_out*(2*nnval)/(exp(2*nnval)-1);

    % P36
    ca_site = ca_out*(a/(2*pi*r*l))*(2*nnval)/(exp(2*nnval)-1);

end

function df = cal_complete_state(t,f)
    % state diagram dynamics calculation
    % f = [X0,X1,X2,X3,X4,X5,Y0,Y1,Y2,Y3,Y4,Y5];
    % ind=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12];
    % df: time derivatives of respective state

    global nn ca2out alpha_v gamma_v_r alpha_s_0 beta_s a_0 r_0 l_0

    alpha_s = alpha_s_0*exp(nn); % P26
    gamma_v = gamma_v_r*ca2out*a_0/(2*pi*r_0*l_0)*(2*nn)/(exp(2*nn)-1) * 1e7; % P38

    df = zeros(length(f),1); % time deriative vector
    % detailed state transition
    % X_0
    df(1) = -5*alpha_s*f(1)                                + 1*beta_s*f(2) + alpha_v*f(7)                  ;
    % X_1
    df(2) = -4*alpha_s*f(2) - 1*beta_s*f(2)                + 2*beta_s*f(3) + alpha_v*f(8)  + 5*alpha_s*f(1);
    % X_2
    df(3) = -3*alpha_s*f(3) - 2*beta_s*f(3)                + 3*beta_s*f(4) + alpha_v*f(9)  + 4*alpha_s*f(2);
    % X_3
    df(4) = -2*alpha_s*f(4) - 3*beta_s*f(4)                + 4*beta_s*f(5) + alpha_v*f(10) + 3*alpha_s*f(3);
    % X_4
    df(5) = -1*alpha_s*f(5) - 4*beta_s*f(5)                + 5*beta_s*f(6) + alpha_v*f(11) + 2*alpha_s*f(4);
    % X_5
    df(6) =                 - 5*beta_s*f(6) - gamma_v*f(6)                 + alpha_v*f(12) + 1*alpha_s*f(5);
    % Y_0
    df(7) = -5*alpha_s*f(7)                                + 1*beta_s*f(8) - alpha_v*f(7)                  ;
    % Y_1
    df(8) = -4*alpha_s*f(8) - 1*beta_s*f(8)                + 2*beta_s*f(9) - alpha_v*f(8)  + 5*alpha_s*f(7);
    % Y_2
    df(9) = -3*alpha_s*f(9) - 2*beta_s*f(9)                + 3*beta_s*f(10)- alpha_v*f(9)  + 4*alpha_s*f(8);
    % Y_3
    df(10)= -2*alpha_s*f(10)- 3*beta_s*f(10)               + 4*beta_s*f(11)- alpha_v*f(10) + 3*alpha_s*f(9);
    % Y_4
    df(11)= -1*alpha_s*f(11)- 4*beta_s*f(11)               + 5*beta_s*f(12)- alpha_v*f(11) + 2*alpha_s*f(10);
    % Y_5
    df(12)=                 - 5*beta_s*f(12) + gamma_v*f(6)                - alpha_v*f(12) + 1*alpha_s*f(11);

end
