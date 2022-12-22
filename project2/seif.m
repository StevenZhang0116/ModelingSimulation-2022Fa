function [vS,ilS,idS] = seif(Iinj,dt)

Smemb = 1000; % [um^2] surface area of the membrane 
Cmemb = 1.0;  % [uF/cm^2] membrane capacitance density 
Cm = Cmemb * Smemb * 1e-8; % [uF] membrane capacitance 

pS = struct();
pS.GL = 0.1; % [mS/cm^2] leak conductance density 
pS.EL = -65.3; % [mV] leak reversal potential 
pS.Vth = -60.2; % [mV] spike threshold 
pS.Kth = 3.5; % [mV] slope factor 
pS.Vsp = +15; % [mV] spike-detecting threshold 
pS.Tref = 2.8; % [ms] refractory period 
pS.gL = pS.GL * Smemb * 1e-8; % leak conductance [mS]

vS = zeros(1,length(Iinj)); % [mV] membrane potential 
ilS = zeros(1, length(Iinj)); % leak current 
idS = zeros(1, length(Iinj)); % leak current 
NrefS = round(pS.Tref/dt); % refractory period in steps 
rCountS = 0; % counter for refractory period 

vS(:,1) = pS.EL;  % initial membrane potential (sEIF model) 

for j=1:length(Iinj)-1
    % ionic currents
    ilS(j) = pS.gL * ( pS.EL-vS(j) ); 
    idS(j) = pS.gL * pS.Kth * exp( (vS(j)-pS.Vth)/pS.Kth ); 

    dvS_dt = ( ilS(j) + idS(j) + Iinj(j)*1e-6 ) ./ Cm;  
    vS(j+1) = dvS_dt * dt + vS(j); 

    % check for threshold crossing
    if(rCountS>0)
        rCountS = rCountS-1;  
        vS(j+1) = pS.EL;
    elseif(vS(j+1)>=pS.Vsp) 
        rCountS = NrefS; 
        vS(j+1) = pS.Vsp;
    end 
end 



