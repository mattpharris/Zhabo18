%Inputs
start = 0; %Seconds
stop = 600; %Seconds
stepSize = 10e-6; %Seconds
altstepSize = 1e-3;
intervals = (stop - start) / stepSize;
subunits = 12;

%Initial Conditions
ek = 20; %[CaMKII] @t=0, .1-30uM
CaM0 = 25; %[CaM] @t=0, 1uM 
ep0 = 1; %[PP1] @t=0, .01-1.2uM
I0 = 0; %[I0] @t=0, 0 or .1uM
Ca2r = .1; %[Ca2+] @t=0, uM
alpha0 = zeros(1,subunits + 8);
alpha0(1) = Ca2r;
alpha0(2) = CaM0;
alpha0(6:7) = [ep0,I0];
alpha0(8) = ek;
alpha = alpha0;

%Constants
kMPP1 = .4; %uM, .4-20uM
kHCaM4_CaMKII = .05; %uM, .05uM: Byrne et al 2009 (J Neurosci)
kHCa_Calc = .7; %uM, .3-1.4uM 

kCaMKII = .5; %s^-1, .5
kPP1 = 2.0; %s^-1, 2.0

kOnI1P = 1.0; %(uM)^-1(s)^-1
kOffI1P = .001; %s^-1
kOnN = 100; %(uM)^-1(s)^-1
kOffN = 750; %s^-1
kOnC = 4; %(uM)^-1(s)^-1
kOffC = 9.25; %s^-1

rateConst = [kMPP1, kHCaM4_CaMKII, kHCa_Calc, kCaMKII, kPP1, kOnI1P, kOffI1P, kOnN, kOffN, kOnC, kOffC];

%Tracking Variables
caSweep = [0:.05:1.5 1.45:-.05:0];
%caSweep = 10; %Only one calcium concentration? Comment...
resP = zeros(size(caSweep));
i = 0;
j = 0;

%Calcium Forcing Function Constants
f = 50; %Hz, Excitation freq range (5-100Hz)
A = 1; %uM, amplitude of Ca2 flux [1 = LTP, .4 = background]
tau = .2; %s, relaxation time of Ca2 decay
ca2Const = [f A tau];

%Calcium Function dy = d[Ca2+] @ time = t
ca2Func = @(t,y,ca2Const)ca2Const(2)/(exp(1/(ca2Const(1)*ca2Const(3))) - 1)*(1 - exp(-t/ca2Const(3)));
ca2Func = @(t,y,ca2Const)0; %Constant Calcium Concentration

for kPP1_ = [.5, 1, 2, 4, 8]
    j = j + 1;
    rateConst(5) = kPP1_;
    alpha = alpha0;
    clear('rxnOut');
    
    for constantCa = caSweep
        i = i + 1;
        if(exist('rxnOut','var'))
            alpha = rxnOut(end,:);
        end
        
        %Adjust Initial [Ca2+]
        alpha(1) = constantCa;
        
        %Execute Model
        [tVec,rxnOut] = Zhabotinsky2000_CaM4(start, stop, stepSize, altstepSize, subunits, ca2Func, ca2Const, alpha, rateConst);
        
        %Calculate total phosphorylated
        Pt = 0;
        for p = 1:subunits
            Pt = Pt + p*rxnOut(end,p+8);
        end
        
        resP(j,i) = Pt;
    end
end

FileName = strcat('OutputDutyCurves_',datetime('Format','MM-DD_HH-mm-ss'));
csvwrite(FileName,[caSweep' resP']);

%Plot The Duty Cycle for Ca2+ sweep
% figure(1)
% hold on
% plot(caSweep,resP)
% xlabel('[Ca2+] (uM)')
% ylabel('Total Phosphorylated Subunits (uM)')