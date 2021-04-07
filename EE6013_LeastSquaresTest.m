%% EE6013_LeastSquaresTest.m

%% Clean Up

clear all;
close all;
clc 

%% Variables

N = 100;
variance = 0.05;
normUpper = 0.95;
normLower = 0.2;

usableCap = 9.3;        %kWh
totalCap = 9.8;         %kWh 100% SoE
dischargeV = 350;       %V
dischargeA = 14.3;      %A
chargeV = 420;          %V
chargeA = 11.9;         %A

%% Battery Aggregation

% Got frustrated trying to code an OO battery since you cannot make a
% matrix of battery objects. Therefore since the application is still
% simple the functions where just added to the end of the code.
battMatrix = batt(N,normUpper,normLower,variance);


%% Other Variables

Ceol = -800000;
ybol = 0.6;
iterationMax = 5;

%% Optimization
% Using Simultaneous OCP (Optimization Control ?)
% Will involve the use of the KKT


% Assume for right now that the cycles will be 1;
f_d = S_t*S_SoC*S_T + S_delta*S_SoC*S_T;

DoD = exp(-f_d);


deltaInit = zeros(N,1);

for i = 1:N
    
    deltaInit(i) = 0.75 ;
    
end

%%

for i = 1:iterationMax
    
    
    
    
    
end



%% Functions


function KKTout = KKT()

    

end


function lagrange = lagrangeCalc(F,lambda,c)

    lagrange = F;
    for i = 1:length(c)
        
       lagrange = lagrange + lambda(i)*c(i);
       
    end

end


function f_d = degredation(S_SoC,S_t,S_delta,cycle_index)

    f_d = S_t*S_SoC*S_T;
    for i = 1:cycle_index
        f_d = f_d + S_delta*S_SoC*S_T;
    end
    
end


function S_t = stressTimeFactor(t)
    
    kt = 1.49*10^-6;
    S_t = kt*t;
    
end


function S_SoC = stressSOCFactor(SoC)
    
    SoCref = 0.5;
    kSoC = 1.04;
    S_SoC = exp(kSoC*(SoC - SoCref));
    
end


function S_delta = stressDODFactor(delta)
    
    e = 0.3;
    d = 0.1;
    c = 0.09;
    b = 0.03;
    a = 0.01;
    S_delta = a*delta^4 + b*delta^3 + c*delta^2 + d*delta + e;
    
end


% 
% function tempStress = stressTempFactor()
% 
% 
%     tempStress = exp(kT*(T - Tref)*(Tref/T));
% 
% 
% end


function loadPattern = cappedLoad(precision,upLim,downLim)

    x = 0:precision:2*pi;
    for i = length(x)
        y = cos(x);
        
        if y > upLim
            
            y = upLim;
            
        elseif y < downLim
            
            y = downLim;
            
        end
        
    end
    
    loadPattern = y;
end


function batteryValues = batt(N,upp,low,variance)

    batteryValues = zeros(N,2);
    for i = 1:N
        
        batteryValues(i,1) = randLim(variance,low);
        batteryValues(i,2) = randLim(variance,upp);
        
    end

end


function limit = randLim(precision,norm)

    limit = norm*(precision*rand);

end



%% Edit
%
% [ECM] 30/01/2021
%   Added basic functions to allow for battery aggregation. 
%