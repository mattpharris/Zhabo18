function [time, rxnOutput] = Zhabotinsky2000_CaM4(a, b, h1, h2, u, ca2Func, ca2Const, alpha)
%base rk4
%differential equations that are written in the form

%dy1/dt = f1(t,y1,y2,...,ym)
%dy2/dt = f2(t,y1,y2,...,ym)
%.
%.
%.
%dym/dt = fm(t,y1,y2,...,ym)

%with t in the interval [a; b], intial conditions alpha, u = #subunits and
%N = #timesteps

%w coefficients
coef = zeros(1,u);
for coefCount = 1:1:u/2
    top = 0;
    bottom = 0;
    for j = 1:1:coefCount
        top = top + j*nchoosek(u/2+1,j);
        bottom = bottom + nchoosek(u/2+1,j);
    end
    coef(coefCount) = top/bottom;
    coef(u-coefCount) = top/bottom;
end

%Rotate input vector if incorrectly input
m = size(alpha,1);
if m == 1
   alpha = alpha';
end

%Calculate step and stuff
h = h1;        %the step size
t = a;
w = alpha;     %initial conditions
CaMstable = false;

%Looping
j = 1;
unstable = true;
while unstable 
    k1 = h*f(t(j), w(:,j),ca2Func,ca2Const,alpha,coef,CaMstable);
    k2 = h*f(t(j)+h/2, w(:,j)+0.5*k1,ca2Func,ca2Const,alpha,coef,CaMstable);
    k3 = h*f(t(j)+h/2, w(:,j)+0.5*k2,ca2Func,ca2Const,alpha,coef,CaMstable);
    k4 = h*f(t(j)+h, w(:,j)+k3,ca2Func,ca2Const,alpha,coef,CaMstable);
    w(:,j+1) = w(:,j) + (k1 + 2*k2 + 2*k3 + k4)/6;
    t(j+1) = t(j) + h;
    if(j > 1000 && sum((w(8:u+8,j+1) - w(8:u+8,j)).^2) < 1e-20)
        unstable = false;
        fprintf('System Stabilized\n')
    elseif(h == h1 && j > 1000 && sum((w(2:u-3,j+1) - w(2:u-3,j)).^2) < 1e-10)
        h = h2;
        CaMstable = true;
        fprintf('CaM Stabilized\n')
    elseif(b && t(j+1) > b)
        unstable = false;
        fprintf('System Failed to Stabilize\n')
    end
    j = j+1;
end

%Plot Rxn Curves
% figure;
% plot(t',w(1:end,:)')
% leg = {'Ca2+','CaM00','CaM2N','CaM2C','CaM4','PP1','I1P'};
% for L = 0:1:u
%     leg{end+1} = strcat('P',num2str(L));
% end
% legend(leg);
% ylim([0,30]);
%drawnow;

time = t';
rxnOutput = w';

function dy = f(t, y, ca2Func, ca2Const, alpha, w, CaMstable)
%% ODE FUNCTION
%v(1) = v1
%v(2) = v2
%v(3) = v3
%v(4) = vPKA
%v(5) = vCaN

%dy(1)/dt = d[Ca2+]/dt

%dy(6)/dt = dep/dt
%dy(7)/dt = dI/dt
%dy(8)/dt = dP0/dt
%...
%dy(8+u)/dt = dPu/dt

%Constants
kMPP1 = .4; %uM, .4-20uM
kHCaM4_CaMKII = .05; %uM, .05uM: Byrne et al 2009 (J Neurosci)
kHCa_Calc = .7; %uM, .3-1.4uM 

kCaMKII = .2; %s^-1, .5
kPP1 = 6.0; %s^-1, 2.0

kOnI1P = 1.0; %(uM)^-1(s)^-1
kOffI1P = .001; %s^-1
kOnN = 100; %(uM)^-1(s)^-1
kOffN = 750; %s^-1
kOnC = 4; %(uM)^-1(s)^-1
kOffC = 9.25; %s^-1

%Determine number of subunits (10,12,14)
u = size(y);
u = max(u)-8;

%Total Phosphorylated subunits
Pt = 0;
for Pcount = 1:1:u
  Pt = Pt + y(Pcount+8)*Pcount;
end        

%% ODEs
dy = zeros(u+8,1);

%%%%%%%%Ca2+%%%%%%%%
%d[Ca2+]
dy(1) = ca2Func(t,y(1),ca2Const);

%Calcium Held Constant or Stoich Responsive
constant = true;
if ~constant
    dy(1) = dy(1) - 2*kOnN*(y(2) + y(4))*y(1)^2 - 2*kOnC*(y(2) + y(3))*y(1)^2;
    dy(1) = dy(1) + 2*kOffN*(y(3) + y(5)) + 2*kOffC*(y(4) + y(5));
end

%%%%%%%%CaM%%%%%%%%
%dCaM0-dCaM4
if ~CaMstable
dy(2:5) = [
    kOffN*y(3) + kOffC*y(4) - kOnN*y(2)*y(1)^2 - kOnC*y(2)*y(1)^2,
    kOnN*y(2)*y(1)^2 + kOffC*y(5) - kOnC*y(3)*y(1)^2 - kOffN*y(3),
    kOnC*y(2)*y(1)^2 + kOffN*y(5) - kOnN*y(4)*y(1)^2 - kOffC*y(4),
    (kOnN*y(4) + kOnC*y(3))*y(1)^2 - (kOffC + kOffN)*y(5);
    ];
else
    dy(2:5) = [0,0,0,0];
end
%%%%%%%REUSED_COEFFICIENTS%%%%%%%%
F = (y(5)/kHCaM4_CaMKII) / (1 + y(5)/kHCaM4_CaMKII);

v = [
    u*kCaMKII*F^2*y(8); %Initiation of Autophosphorylation, probability two neighbors are bound to CaM
    kCaMKII*F; %Propogation of phosphorylation
    (kPP1*y(6))/(kMPP1 + Pt); %Dephosphorylation by M-M
    1; %vCaN
    1  %vPKA
    ];

%%%%%%%%PP1%%%%%%%%
dy(6:7) = [
     -kOnI1P*y(6)*y(7) + kOffI1P * (alpha(6) - y(6));
     -kOnI1P*y(6)*y(7) + kOffI1P * (alpha(6) - y(6)) + v(4)*(alpha(7) - y(7)) - (v(5)*(y(1)/kHCa_Calc)^3*y(7))/(1 + (y(1)/kHCa_Calc)^3);
     ];

%%%%%%%%CaMKII%%%%%%%%
%CaM exclusion (Comment out to remove rule), Essentially Determines a
%fraction "F" from the publication of CaMKII units occupied by CaM/Ca2+, He
%assumes that this proportion is essentially uneffected by the
%phosphorylation state of the subunit, and as such can be a simple scale
%factor. This equation thus reduces the calculated rate of
%dephosphorylation, v3, scaling by a factor of (1-F).
%v(3) = v(3)*(1-(((Ca2/kH1)^4)/(1 + (Ca2/kH1)^4)));

%P0-P1
dy(8:9) = [
     -v(1) + v(3)*y(9);
     v(1) - v(3)*y(9) - v(2)*y(9) + 2*v(3)*y(10)
     ];

%P2-P(U-1)
for p = 1:(u-2)
     dy(p + 9) = w(p)*v(2)*y(p+8) - (p+1)*v(3)*y(p+9) - w(p+1)*v(2)*y(p+9) + (p+2)*v(3)*y(p+10);
end

%PUz
dy(end) = v(2)*y(u+7) - u*v(3)*y(u+8);
