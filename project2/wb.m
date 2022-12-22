% Wang-Buzsaki Model

function [vW,inW,ikW,ilW] = wb(Iinj,dt)
    Smemb = 1000; % [um^2] surface area of the membrane 
    Cmemb = 1.0;  % [uF/cm^2] membrane capacitance density 
    Cm = Cmemb * Smemb * 1e-8; % [uF] membrane capacitance 
    
    % model-specific parameters
    pW = struct();
    pW.GL = 0.1; % [mS/cm^2] leak conductance density 
    pW.GN = 35.0; % [mS/cm^2] Na conductance density
    pW.GK = 15.0; % [mS/cm^2] K conductance density 
    pW.EN = +55.0; % [mV] Na reversal potential 
    pW.EK = -90.0; % [mV] K reversal potential 
    pW.EL = -65.0; % [mV] leak reversal potential 
    pW.gN = pW.GN * Smemb * 1e-8; % Na conductance [mS]
    pW.gK = pW.GK * Smemb * 1e-8; % K conductance [mS]
    pW.gL = pW.GL * Smemb * 1e-8; % leak conductance [mS]
   
    % vectors for storing variables
    vW = zeros(1, length(Iinj)); % [mV] membrane potential 
    mW = zeros(1, length(Iinj)); % Na activation variable
    hW = zeros(1, length(Iinj)); % Na inactivation variable
    nW = zeros(1, length(Iinj)); % K activation variable 
    inW = zeros(1, length(Iinj)); % Na current
    ikW = zeros(1, length(Iinj)); % K current
    ilW = zeros(1, length(Iinj)); % leak current 
    
    % initial values
    vW(:,1) = pW.EL;  % initial membrane potential (WB model) 
    mW(:,1) = WBalphaM(pW.EL) / (WBalphaM(pW.EL) + WBbetaM(pW.EL)) ; % initial m
    hW(:,1) = WBalphaH(pW.EL) / (WBalphaH(pW.EL) + WBbetaH(pW.EL)) ; % initial h
    nW(:,1) = WBalphaN(pW.EL) / (WBalphaN(pW.EL) + WBbetaN(pW.EL)) ; % initial n
    
    for j=1:length(Iinj)-1
        % ionic currents (WB model) g[mS] * V[mV] = I[uA]
        inW(j) = pW.gN .* mW(:,j).^3 .* hW(:,j) .* ( pW.EN - vW(j) ); 
        ikW(j) = pW.gK .* nW(:,j).^4            .* ( pW.EK - vW(j) ); 
        ilW(j) = pW.gL                          .* ( pW.EL - vW(j) ); 
    
        % derivatives I[uA] / C[uF] * dt[ms] = dv[mV]    
        dvW_dt = ( inW(j) + ikW(j) + ilW(j) + Iinj(j)*1e-6 ) ./ Cm;  
        dmW_dt = (1-mW(j)).* WBalphaM(vW(j)) - mW(j).*WBbetaM(vW(j));
        dhW_dt = (1-hW(j)).* WBalphaH(vW(j)) - hW(j).*WBbetaH(vW(j));
        dnW_dt = (1-nW(j)).* WBalphaN(vW(j)) - nW(j).*WBbetaN(vW(j));
    
        % calculate next step 
        vW(j+1) = dvW_dt * dt + vW(j); 
        mW(:,j+1) = mW(:,j) + dmW_dt * dt; 
        hW(:,j+1) = hW(:,j) + dhW_dt * dt; 
        nW(:,j+1) = nW(:,j) + dnW_dt * dt; 
    end 

function x = WBalphaM(v)
  x = 5 * 0.1 * (v+35) ./ ( 1 - exp(-(v+35)/10) ); 
function x = WBalphaH(v)
  x = 0.35 * exp(-(v+58)/20);
function x = WBalphaN(v)
  x = 0.05 * (v+34) ./ ( 1 - exp(-(v+34)/10) );
function x = WBbetaM(v)
  x = 5 * 4.0 * exp(-(v+60)/18);
function x = WBbetaH(v)
  x = 5.0 ./ ( 1 + exp(-(v+28)/10) ); 
function x = WBbetaN(v)
  x = 0.625 * exp(-(v+44)/80);
