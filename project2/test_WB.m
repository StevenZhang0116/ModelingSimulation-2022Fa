% test of basic WB model
% e.g. try with different injected current
% mainly used to graph illustration

close all

dt = 0.001; % [us]

Tleng = 1000; % simulation time length [ms]
TstimS = 0; % stimulus starting time [ms]
TstimE = 1000; % stimulus ending time [ms]
t = 0:dt:Tleng; % time vector 

iibreak = 0.1;
iiiamp = 0:iibreak:200; % input current set
lll = length(iiiamp);
spike_cnt1 = zeros(lll,1);
spike_cnt2 = zeros(lll,1);

% figure()
% hold on
parfor (i = 1:length(iiiamp)) % parallel
    iamp = iiiamp(i); % [pA]
%     disp(iamp)
    Iinj = zeros(1, length(t));
    Iinj((t>=TstimS)) = iamp; 
    Iinj((t>=TstimE)) = 0; 
    [vW,inW,ikW,ilW] = wb(Iinj,dt);
    [vS,ils,ids] = seif(Iinj,dt);
%     plot(t-TstimS,vW,'LineWidth',2)
    cntspike1 = findpeaks(vW);
    spike_cnt1(i) = length(cntspike1);
    cntspike2 = findpeaks(vS);
    spike_cnt2(i) = length(cntspike2);
end
% legend('10','30','50','70','90','FontSize',12)
% xlabel('Time (ms)','FontSize',12)
% ylabel('Potential (mV)','FontSize',12);

figure()
hold on
plot(iiiamp/10,spike_cnt1,'LineWidth',2)
plot(iiiamp/10,spike_cnt2,'LineWidth',2)
legend('WB','SEIF','FontSize',12)
xlabel('I_{app} (\mu A/cm^2)','FontSize',12)
ylabel('f (Hz)','FontSize',12)
hold off

% F-I relationship curve fitting
x = iiiamp/10;
y = spike_cnt1; y = y';
f = @(b,x) b(1).*exp(b(2).*x)+b(3);  
B = fminsearch(@(b) norm(y - f(b,x)), [-200; -1; 100]);

plott = 0;
adj = 5;
if plott == 1
    % plotting
    figure(1); 
    cla; hold on; 
    plot(t-TstimS,vW,'LineWidth',2);
    plot(t-TstimS,vS,'LineWidth',2)
    xlim([-5,15]);
    ylim([-80,30]);
    legend('WB','SEIF','FontSize',12)
    ylabel('Potential (mV)','FontSize',12);
    xlabel('Time (ms)','FontSize',12);

    figure(2)
    cla; hold on;
    plot(t-TstimS,Iinj/10,'LineWidth',2)
    xlabel('Time (ms)','FontSize',12);
    ylabel('Current (nA)','FontSize',12)

    figure(3)
    cla; hold on;
    plot(t-TstimS,(inW+ikW)*1000,'LineWidth',2)
    xlim([3,8]);
    xlabel('Time (ms)','FontSize',12);
    ylabel('Current (nA)','FontSize',12)

end








