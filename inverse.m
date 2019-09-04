function inverse
% main program for the parameter estimation using the DE algorithm.
% subroutines: differential_evolution.m (DE algorithm)
%              eval_objective.m (objective function)
%              graphs (generate the graphs of the simulated TGA and DTG)
% The objective of this program is to estimate the cinematic parameters 
% (k and E) for the thermal decomposition of lignocellulosic material in [
% term os if pseudocomponents: hemicellulose (h), cellulose (c) and lignin (L).
% k is the pre-exponential constant and E stands for the activation Energy.

clc
clear all
clear global
close all
format long

initial_time = cputime; 
VTR = 1.e-10; 
D = 12;            % quantidade de variáveis de projeto (número parâmetros estimados)

xref = [  7.8669337e+007  1.1092011e+005  5.2106049e+015  2.0695492e+005  5.2609923e+003  7.6461744e+004  2.1971877e-001  3.8957752e-001  3.8437544e-001  1.0000000e+000  1.0000000e+000  3.0000000e+000];
% FIT for total TG = 0.267898
% FIT for total DTG = 1.053254

 XVmin = [1*xref(1) 1*xref(2) 1*xref(3) 1*xref(4) 1*xref(5) 1*xref(6) 1*xref(7) 1*xref(8) 1*xref(9) 1*xref(10) 1*xref(11) 1*xref(12)];
 XVmax = [1*xref(1) 1*xref(2) 1*xref(3) 1*xref(4) 1*xref(5) 1*xref(6) 1*xref(7) 1*xref(8) 1*xref(9) 1*xref(10) 1*xref(11) 1*xref(12)];  

%  XVmin = [0.95*xref(1) 0.95*xref(2) 0.95*xref(3) 0.95*xref(4) 0.95*xref(5) 0.95*xref(6) 1*xref(7) 1*xref(8) 1*xref(9) 1*xref(10) 1*xref(11) 1*xref(12)];
%  XVmax = [1.05*xref(1) 1.05*xref(2) 1.05*xref(3) 1.05*xref(4) 1.05*xref(5) 1.05*xref(6) 1*xref(7) 1*xref(8) 1*xref(9) 1*xref(10) 1*xref(11) 1*xref(12)];

%  XVmin = [0.8*xref(1) 0.9*xref(2) 0.8*xref(3) 0.9*xref(4) 0.9*xref(5) 0.9*xref(6) 1*xref(7) 1*xref(8) 1*xref(9) 1*xref(10) 1*xref(11) 1*xref(12)];
%  XVmax = [1.2*xref(1) 1.1*xref(2) 1.2*xref(3) 1.1*xref(4) 1.1*xref(5) 1.1*xref(6) 1*xref(7) 1*xref(8) 1*xref(9) 1*xref(10) 1*xref(11) 1*xref(12)];

NP = 20;           % population size
itermax = 200;     % max number of iterations
F = 0.8;           % perturbation ratio
CR = 0.8;          % crossover probability
strategy=7;        % mutation strategy 
refresh=10;

load data_carpel_5.txt;   %set of experimental data
                          %(obs: this file is built from the load.m
                          %program)
                   
y=data_carpel_5;


% DE algorithm
[X,FO,NF]=differential_evolution('eval_objective',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh);

limits=[XVmin; X; XVmax;];
save limits.txt limits -ASCII
open limits.txt

final_time = 1/60*(cputime-initial_time);
fprintf(1,'  time of execution (in minutes) = %f\n',final_time);

graphsDTG(X,y)
