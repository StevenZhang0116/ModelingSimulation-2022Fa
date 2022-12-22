% test of myelinated axon models using WB model
% for illustration purpose
close all

Ndef = 141; % number of nodes

axonrange = 0.1:0.1:3;
valiter = [];

for i = 1:length(axonrange)
    Ddef = 2.0; % [um] axonal diameter 
    Ldef = 2.0; % [um] nodal length 
    Idef = 200; % [um] internodal length 
    
    % time vector
    Tleng = 25; % simulation time length [ms]
    TstimS = 10; % stimulus starting time [ms]
    TstimE = 11; % stimulus ending time [ms]
    t = 0:dt:Tleng; % time vector 
    
    iamp = 100; % [pA]
    Iinj = zeros(Ndef, length(t));
    Iinj(20,(t>=TstimS)) = iamp; % [pA]
    Iinj(20,(t>=TstimE)) = 0; 
    
    vWB = axonwb(Iinj,dt,Ndef,Ddef,Ldef,Idef);
    
    [~,i1] = max(vWB(40,:));
    [~,i2] = max(vWB(90,:));
    vel_mWB = ((Ldef+Idef)*50/1000)/(t(i2)-t(i1));
    
    % plotting
    plott = 1;
    if plott == 1
        figure()
        cla; hold on; 
        plot(t-TstimS,vWB(20,:),'-','LineWidth',2,'color',[0.5,0.5,0.5]); 
        plot(t-TstimS,vWB(30,:),'r-','LineWidth',2); 
        plot(t-TstimS,vWB(40,:),'y-','LineWidth',2); 
        plot(t-TstimS,vWB(50,:),'g-','LineWidth',2); 
        plot(t-TstimS,vWB(60,:),'c-','LineWidth',2); 
        plot(t-TstimS,vWB(70,:),'b-','LineWidth',2); 
        plot(t-TstimS,vWB(80,:),'m-','LineWidth',2); 
        plot(t-TstimS,vWB(90,:),'k-','LineWidth',2); 
        legend('20','30','40','50','60','70','80','90')
        xlim([-1,+6]);
        ylim([-100,40]);
        ylabel('potential (mV)');
        xlabel('Time (ms)')
        title(sprintf('WB model (myelinated axon): %.2f [m/s]',vel_mWB));
    end
    valiter(end+1) = vel_mWB;
end

figure()
plot(axonrange,valiter,'LineWidth',2)
xlabel('Axon diameter (µm)','FontSize',12)
ylabel('Conduction velocity (m/s)','FontSize',12)




