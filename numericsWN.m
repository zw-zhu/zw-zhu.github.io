%% Numerical implemations - Section 2.4, 2.5, 2.6
% One program to do it all: Euler, Heun and Runge-Kutta
%
% Instructions: 
% (1) set your initial data below
% (2) set the step size and number of iterations, or alternatively set the ending time
% t_end to get N
% (3) set/modify the function definition "f" to fit the DE problem
% (4) run the program (click the green triangle) to call the three
% numerical algorithms to plot results. Label and save your graphs
% appropriately.


function [] = numericsWN()
clc

%%% initial data
t0 = 0;
y0 = 0;

%%% set number of steps and step size
h = .05;
%N = 2;
%t_end = N*h
%% or set tend to find N:
t_end = 5;
N = round(t_end/h);


%%% Call numerical algorithm, either Euler, Heun, or Runge-Kutta, or all of
%%% them!

[tE,yE] = Euler(@f,t0,y0,h,N);
[tH,yH] = Heun(@f,t0,y0,h,N);
[tRK,yRK] = RK(@f,t0,y0,h,N);

figure(1)
plot(tE,yE,'k', tH,yH, 'r',tRK,yRK,'b','linewidth',2)
xlabel('t','FontSize',16)
ylabel('y','FontSize',16)
set(gca,'FontSize',16)
title('Label your output appropriately')
legend('Euler method', 'Improved Euler', 'Runge-Kutta method')
hold on

yRK(end)

%%%%% sub-function definitions


function [fOut] = f(t,y)
fOut = 32-0.2*y-0.15*y^(1.5);

function [err] = Error(ya,yb)
err = abs((ya-yb)./yb);




%% Euler method
% Solve the IVP y'(x)=f(x,y) with y(y0)=x0. The function f is given as a
% symbolic function.
% Given: step size h and number of steps. 

function [xE,yE] = Euler(f,x0,y0,h,N)

xE = zeros(1,N+1); yE = zeros(1,N+1);   % Initialization
xE(1) = x0; yE(1) = y0;                 % Initial condition

for i=1:N                               % For loop over the steps
    xE(i+1) = xE(i)+h;                  % Next x
    yE(i+1) = yE(i)+h*f(xE(i),yE(i));   % Euler iteration - next y
end


%% Heun method
% Solve the IVP y'(x)=f(x,y) with y(y0)=x0. The function f is given as a
% symbolic function.
% Given: step size h and number of steps. 

function [xH,yH] = Heun(f,x0,y0,h,N)

xH = zeros(1,N+1); yH = zeros(1,N+1);   % Initialization 
xH(1) = x0; yH(1) = y0;                 % Initial condition

for i=1:N                               % For loop over the steps
    k1 = f(xH(i),yH(i));                % k1
    k2 = f(xH(i)+h,yH(i)+h*k1);         % k2
    k = (k1+k2)/2;                      % k
    xH(i+1) = xH(i)+h;                  % Next x
    yH(i+1) = yH(i)+h*k;                % Heun iteration - next y
end


%% RK method
% Solve the IVP y'(x)=f(x,y) with y(y0)=x0. The function f is given as a
% symbolic function.
% Given: step size h and number of steps. 

function [xRK,yRK] = RK(f,x0,y0,h,N)

xRK = zeros(1,N+1); yRK = zeros(1,N+1); % Initialization
xRK(1) = x0; yRK(1) = y0;               % Initial condition

for i=1:N                               % For loop over the steps
    k1 = f(xRK(i),yRK(i));              % k1
    k2 = f(xRK(i)+h/2,yRK(i)+h*k1/2);   % k2       
    k3 = f(xRK(i)+h/2,yRK(i)+h*k2/2);   % k3
    k4 = f(xRK(i)+h,yRK(i)+h*k3);       % k4
    k = (k1+2*k2+2*k3+k4)/6;            % k
    xRK(i+1) = xRK(i)+h;                % Next x
    yRK(i+1) = yRK(i)+h*k;              % Runge Kutta iteration - next y
end
