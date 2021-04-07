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
% battMatrix = batt(N,normUpper,normLower,variance);


%% Other Variables

Ceol = -800000;
ybol = 0.6;
maxIteration = 50;


kt = 1.49*10^-8;
SoCref = 0.7;
Tref = 0.25;
kT = 0.0693;
kSoC = 1.04;
e = 3*10^-4;
d = 0.4*10^-4;
c = 10*10^-6;
b = 9*10^-6;
a = 5*10^-6;

S_T = exp(kT*((1E-4)-Tref)*(Tref/(1E-4)));

init = [0.15;0.19;0.001;0.001];

%% Optimization
% Using Simultaneous OCP (Optimization Control ?)
% Will involve the use of the KKT

N = 1;
maxIteration = 10;

functionVal = zeros(maxIteration,N);

delta = zeros(maxIteration,N);
init = [0.75];
delta(1,:) = init(:);

SoC = zeros(maxIteration,N);
init = [0.6];
SoC(1,:) = init(:);

time = zeros(maxIteration,N);
init = 0.1:(1-0.1)/maxIteration:1;
time(:,1) = init(:);

deltaSoH = zeros(maxIteration,N);
init = [0.6];
deltaSoH(1,:) = init(:);

%% Running

for i = 1:maxIteration                                                      % Recursively update 10x Note that since my costs are quadratic, it should only take 1 iteration to be optimal
        
        [x1,x2,u] = splitx(x);                                              % Split up the variable x into states x1 and x2, and control u
        c = calc_c(x1,x2,u,dt);                                             % Calculate the constraints
        G = calc_G(x1,dt);                                                  % Calculate the Jacobian of the constraints with respect to the optimizable elements
        g = calc_g(x1,lam,u,dt);                                            % Calculate the gradient of the Lagrangian
        H = calc_H(x1,lam,dt);                                              % Calculate the Hessian of the Lagrangian

        KKT = [H G'; G zeros(2*N+2)];                                       % Form the KKT matrix
        a = -KKT\[g;c];                                                     % Solve for [-delta_x; lambda]

        dx = a(1:3*(N-1)+2);                                                % Seperate answer into negative delta x and lambda
        
        dlam = a(3*(N-1)+3:end);
        x = x + dx;                                                         % Update optimizable elements
        lam = lam + dlam;
        
%         plot(x1); hold on;                                                % Show the updated position
 end


%% Plots


figure('name','Mesh - Surface Plot')
[SoC,delta] = meshgrid(-50:0.5:50,-50:0.5:50);

S_SoC = exp(kSoC*(SoC - SoCref));
S_delta = a*delta.^4 + b*delta.^3 + c*delta.^2 + d*delta + e;
f_d = S_t*S_SoC*S_T + S_delta*S_SoC*S_T;
SoH = exp(-f_d);
F = SoH - kt*S_SoC*exp(-f_d);

C = SoC.*delta;
surf(SoC,delta,SoH,C)
hold on
xlabel('SoC');
ylabel('DoD');
zlabel('SoH');
% plot3(outputs(:,1),outputs(:,2),outputs(:,3),'r','LineWidth',2);
hold off


%% Functions
function c = calc_c(x1,x2,u,dt)
% Constraints

    N = length(x1);
    c = zeros(2*(N-1)+4,1);                                                 % 2 constraints per time-step, plus 4 boundary constraints
    
    for i = 1:N-1
        j = 2*(i-1) + 1;
        c(j)   = x1(i+1) - x1(i) - x2(i)*dt;                                % Enforce the positions to match
        c(j+1) = x2(i+1) - x2(i) - sin(x1(i))*dt - u(i)*dt;                 % Enforce the velocity to match
    end
     
    j = j+2;                                                                % Initial condition
    c(j) = x1(1) - pi;
    c(j+1) = x2(1);
    
    j = j+2;                                                                % Terminal condition
    c(j) = x1(end);
    c(j+1) = x2(end);
end
    

function G = calc_G(x1,dt)
% Jacobian of constraints with respect to optimizable variables - note that it's rectangular

    N = length(x1);
    rs = 2*(N-1) + 4;                                                       % Number of rows = 2*number of time-steps + 4 boundary conditions
    cs = 3*N-1;                                                             % Number of columns = number of x1 (m) + x2 (m) + u (m-1)
    G = zeros(rs,cs);

    for i = 1:N-1                                                           % Note how this process creates a banded matrix
        n = 3*(i-1) + 1;
        j = 2*(i-1) + 1;

        G(j,  n:n+3) = [-1             -dt  0  1];
        G(j+1,n:n+4) = [-cos(x1(i))*dt -1  -dt 0 1];
    end
    j = j+2;        % Initial conditions
    G(j,1) = 1;     
    G(j+1,2) = 1;
    
    j = j+2;        % Terminal conditions
    G(j,end-1) = 1; 
    G(j+1,end) = 1;
end


function g = calc_g(x1, lam, u, dt)
% Gradient of the Lagrangian
    N = length(x1);
    rc = 3*N - 1;
    g1 = zeros(rc,1);
    g2 = g1;
    g3 = g1;
    

    % Components caused by cost-function F:
    for i = 1:N-1
        n = 3*(i-1) + 1;
        g3(n) = 2*x1(i);
        g3(n+2) = 2*u(i);
    end

    % Components caused by dynamic constraints
    i=1;
    j = 2*(i-1) + 1;
    n = 3*(i-1) + 1;
    g1(n) = -lam(j) + lam(j+1)*dt*cos(x1(i));
    g1(n+1) = -lam(j)*dt - lam(j+1);
    g1(n+2) = -lam(j+1)*dt;
    
    for i = 2:N-1
        j = 2*(i-1) + 1;
        n = 3*(i-1) + 1;
        g1(n) = lam(j-2)-lam(j) + lam(j+1)*dt*cos(x1(i));
        g1(n+1) = lam(j-1)-lam(j)*dt - lam(j+1);
        g1(n+2) = -lam(j+1)*dt;
    end
    
    % Components caused by boundary constraints
    g2(1) = lam(2*N-1);
    g2(2) = lam(2*N);
    g2(end-1) = lam(2*N+1);
    g2(end) = lam(2*N+2);

    % Sum them all up
    g = g1 + g2 + g3;
end     
   

function H = calc_H(x1,lam,dt)
% Hessian of the Lagrangian
    N = length(x1);
    rc = 3*N - 1;
    
    H = zeros(rc);
    
    for i = 1:N-1 % Note how sparse this diagonal matrix is
        n = 3*(i-1) + 1;
        j = 2*(i-1) + 1;
        
        H(n,n) = 2 - lam(j+1)*dt*sin(x1(i));
        H(n+2,n+2) = 2;
    end    
end


function [x1,x2,u] = splitx(x)
% Split the optimizable variables up into their components
    N = length(x);
    N = (N-2)/3 + 1;
    x1 = [];
    x2 = [];
    u = [];
    for i = 1:N-1
        n = 3*(i-1) + 1;
        x1 = [x1 x(n)];
        x2 = [x2 x(n+1)];
        u =  [u  x(n+2)];
    end
    i = N;
    n = 3*(i-1) + 1;
    x1 = [x1 x(n)];
    x2 = [x2 x(n+1)];
end


%% Edit
%
% [ECM] 30/01/2021
%   Added basic functions to allow for battery aggregation. 
%