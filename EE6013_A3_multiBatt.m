%% EE6013_LeastSquaresTest.m

%% Clean Up

clear all;
close all;
clc 

%% Variables

N = 1;
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
maxIteration = 2;
S_T = 1*10^-14;

kt = 1.49*10^-6;
SoCref = 0.5;
% Tref = 0.25;
% kT = 0.0693;
kSoC = 1.04;
e = 6.16*10^-2;
d = 5.37*10^-1;
c = 3.3209;
b = -6.8292;
a = 5.79;

init = [20;25;1;1];

%% Optimization
% Using Simultaneous OCP (Optimization Control ?)
% Will involve the use of the KKT

syms SoC delta DoD_dot t

var = {SoC;delta;t;DoD_dot};

S_t = kt;
S_SoC = exp(kSoC*(SoC - SoCref));
S_delta = a*delta^4 + b*delta^3 + c*delta^2 + d*delta + e;
% S_T = exp(kT*(T - Tref)*(Tref/T));
% Assume for right now that the cycles will be 1;
f_d = S_t*S_SoC*S_T + S_delta*S_SoC*S_T;

DoD = exp(-f_d);
F = DoD_dot + kt*S_SoC*exp(-f_d);

deltaInit = zeros(N,1);
SoCInit = zeros(N,1);
timeInit = zeros(N,maxIteration);
DoD_dot_Init = DoD_dot*ones(maxIteration,1);
% TInit = zeros(N,1);
si = length(var)-1 + maxIteration;
lambda = 2*10^(-3)*ones(si,1);
for i = 1:N
    
    deltaInit(i) = 0.75 ;
    SoCInit(i) = 0.5;
    timeInit(i,1) = 0.1;
%     timeInit(i,maxIteration) = 60*60*24;
    
end



deltaInit = delta - deltaInit;
SoCInit = SoC - SoCInit;
timeInit = t - timeInit;

% for i = 1:maxIteration
%     
%     DoD_dotInit(i) = DoD_dot - F
%     
% end

constraints = [SoCInit;deltaInit;t;DoD_dot_Init];

sumVal = (lambda.*constraints);

lagrange = F + sum(sumVal);

hess = hessian( lagrange ,[delta, SoC, t, DoD_dot]);

grad = gradient( DoD_dot ,[delta, SoC, t, DoD_dot]);

jac = jacobian( constraints ,[delta, SoC, t, DoD_dot]);

numZeros = size(jac);
zer = zeros(numZeros(1),numZeros(1));

largeMat = [hess,jac';jac,zer];

outputs = zeros(maxIteration,length(constraints));

for i = 1:maxIteration
    
    f_d = S_t*S_SoC*S_T + S_delta*S_SoC*S_T;
    
%     for k = 1:i
%         f_d = f_d + S_delta*S_SoC*S_T;
%     end
    
    F = DoD_dot - kt*S_SoC*exp(-f_d);
    
%     constraints((length(var) - 1) + maxIteration) = DoD_dot - F;
        
    sumVal = (lambda.*constraints);
    lagrange = F + sum(sumVal);
    hess = hessian( lagrange ,[delta, SoC, t, DoD_dot]);
    grad = gradient( F ,[delta, SoC, t, DoD_dot]);
    jac = jacobian( constraints ,[delta, SoC, t, DoD_dot]);
    
    if i == 1
        
        aPrev = init;
        
        h = double(runSim(hess,var,aPrev));
        j = double(runSim(jac,var,aPrev));
        g = double(runSim(grad,var,aPrev));
        con = double(runSim(constraints,var,aPrev));
        KKTout = KKT(j,h,g,con);
        
    else
        
        h = double(runSim(hess,var,aPrev));
        j = double(runSim(jac,var,aPrev));
        g = double(runSim(grad,var,aPrev));
        con = double(runSim(constraints,var,aPrev));
        KKTout = KKT(j,h,g,con);
        
    end
    
    for j = 1:length(init)
        aNext(j) = init(j) + KKTout(j);
    end
    for j = 1:length(lambda)
        lambdaNext(j) = KKTout(j + length(init));
    end
    
    SoH_val = exp(-f_d);
    SoH_dot = -kt*S_SoC*exp(-f_d);
    SoH_val = double(runSim(SoH_val,var,aPrev));
    SoH_dot = double(runSim(SoH_dot,var,aPrev));
    outputs(i,1) = aPrev(1);
    outputs(i,2) = aPrev(2);
    outputs(i,3) = SoH_val;
    outputs(i,4) = SoH_dot;
    aPrev = aNext';
    lambda = lambdaNext';
end

optimVaL = double(runSim(f_d,var,aPrev));
SoH = exp(-optimVaL);



%% Plots


figure('name','Mesh - Surface Plot')
[SoC,delta] = meshgrid(-5:0.05:5,-5:0.05:5);

S_SoC = exp(kSoC*(SoC - SoCref));
S_delta = a*delta.^4 + b*delta.^3 + c*delta.^2 + d*delta + e;
f_d = S_t*S_SoC*S_T + S_delta*S_SoC*S_T;
SoH = exp(-f_d);
F = - kt*S_SoC*exp(-f_d);

C = F;
surf(SoC,delta,SoH,C)
hold on
xlabel('SoC');
ylabel('DoD');
zlabel('SoH');
colorbar
plot3(outputs(:,1),outputs(:,2),outputs(:,3),'r',outputs(1,1),outputs(1,2),outputs(1,3),'LineWidth',2);
hold off


figure('name','Mesh - Surface Plot Change in SoH')
[SoC,delta] = meshgrid(-50:0.5:50,-50:0.5:50);

S_SoC = exp(kSoC*(SoC - SoCref));
S_delta = a*delta.^4 + b*delta.^3 + c*delta.^2 + d*delta + e;
f_d = S_t*S_SoC*S_T + S_delta*S_SoC*S_T;
DoD = exp(-f_d);
F = - kt*S_SoC*exp(-f_d);

C = DoD;
surf(SoC,delta,F,C)
hold on
xlabel('SoC');
ylabel('DoD');
zlabel('Delta SoH');
plot3(outputs(:,1),outputs(:,2),outputs(:,4),'r','LineWidth',2);
hold off


%% Functions
function [aFinal,lambda] = newtonMethodKarushKahnTucker(delta,c,aPrev,cOriginal,var,lamPrev,maxIteration)
%newtonMethod Performs newton's method on a given set of equations

    syms x y z k v u 
    
    hessian = ([2,0;0,2]);
    kahn = x^2 + y^2 - lamPrev*(x - y + 5 );
    gradient = ([2*x-lamPrev; 2*y - lamPrev]);
    jac = jacobian([ kahn ],[x, y]);
    outT = [hessian, jac.'; jac, 0];
    useMat= inv(outT);

    cPrev = double(runSim(cOriginal,var, aPrev));

    for i = 1:maxIteration
        if (i > 1)&&((cPrev == 0)||((cPrev > 0)&&(cPrev <= delta))||((cPrev < 0)&&(cPrev >= -delta)))
%             converge = i + 1;
            aFinal = aPrev;
            lambda = lamPrev;
        else
            useMat = double(runSim(useMat,var,aPrev));
            g = double(runSim(gradient,var,aPrev));
            c = double(runSim(c,var,aPrev));
            other = [g.'; c];
            
            created = useMat.*other;
            
            for i = 1:3
                out(i) = sum(created(i,:));
            end
            delta_a = -out(1);
            lamPrev = [out(2),out(3)];
            aNext = aPrev + delta_a';
            cNext = runSim(cOriginal, var, aNext);
            cPrev = double(cNext);
            aPrev = double(aNext);
            aFinal = double(aPrev);
            lambda = lamPrev;
            
            gradient = ([2*x-lamPrev; 2*y - lamPrev]);
            kahn = x^2+y^2-lamPrev*(x-y+5);
            jac = jacobian([ kahn ],[x, y]);
            
            outT = [hessian, jac.'; jac, 0];
            useMat= inv(outT);
        end
    end

%     printOut(converge,var,aFinal,cOriginal)

end


function KKTout = KKT(jac,hess,grad,c)

    syms delta SoC t DoD_dot

    numZeros = size(jac);
    if numZeros(1) ~= length(hess)
        useVar = numZeros(1);
    else
        useVar = numZeros(2);
    end
    zer = zeros(useVar,useVar);

    largeMat = [hess,jac';jac,zer];
    
    dim = size(grad);
    if dim(1) ~= 1
        g = [grad;c];
    else
        g = [grad';c];
    end
    
    KKTout = inv(largeMat)*g;
    
end


function cNext = runSim(cOriginal, var, aNext)
    
    cNext = subs(cOriginal,var,aNext);
%     cNext = cNext';

end


function printOut(converge,var,aFinal,cOriginal)

    fprintf('The Newton method converged in %d iterations.\n\n',converge);
    fprintf('The resulting values\n\n');
    for in = 1:length(aFinal)
        fprintf('%d:\t%6.2f\n',in,aFinal(in));
    end

    fprintf('\nTo confirm that the found values work.\nSubstitute the above into the equations.\n');
    cCheck = runSim(cOriginal, var, aFinal);
    fprintf('The result is:\n');
    for in = 1:length(cCheck)
        fprintf('C%d:\t%6.2f\n',in,cCheck(in));
    end

end


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