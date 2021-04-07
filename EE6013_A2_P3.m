%% EE6013_A2_P3.m

%% Clean 

clc
close all
clear all



%% Using Newton's Method

syms x y

eq1 = -x+exp(y);

eq2 = (-x^4-y.^4+4);

var = {x, y};

cOriginal = ([eq1,eq2]);

g = jacobian([eq1,eq2],[x, y]);

maxIteration = 20;
start = 0.05;
aPrev = double(start*ones(1,length(g)));

aFinal = newtonMethod(g,cOriginal,var,aPrev,maxIteration);

%% Variables

num = -3:0.01:3;
y1 = -3:0.01:1;
y2 = num;

eq1_sketch = exp(y1);

eq2_sketch = (4 - y2.^4).^(1/4);

%% Plot

figure('name','Mapped Equations');
plot(eq1_sketch,y1,eq2_sketch,y2);
xlabel('X');
ylabel('Y');
legend('x = exp(y)','x^4 + y^4 = 4');
grid on;



