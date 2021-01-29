function [t,y,rate,inhibition] = ADM1_R4_mass_output(t,x,s,input,parameter) 

% -------------------------------------------------------------------------
%
% Simplified and mass-based
% Anaerobic Digestion Model No. 1 (ADM1)
%
% ADM1-R4
%
% Implementation in Matlab 
% R2019b (The MathWorks, Inc.)
%
% Matlab ODE function 
%
% Version 1.0
%
% https://github.com/soerenweinrich/ADM1
%
% Copyright (c) 2021 Sören Weinrich
% E-Mail: soeren.weinrich@dbfz.de
%
% Citation example:
%
% Weinrich, S.; Nelles, M. (2021).
% Systematic simplification of the Anaerobic Digestion Model No. 1 (ADM1) -
% Model development and stoichiometric analysis.
% Submitted to Bioresource Technology.
%
% -------------------------------------------------------------------------

% System parameters  
 
V_liq = s(1); 
V_gas = s(2); 
p_atm = s(4); 
 
% Model parameters 
 
K_H_ch4 =  parameter(1); 
K_H_co2 =  parameter(2); 
R =  parameter(3); 
T =  parameter(4); 
k_La =  parameter(5); 
k_ch =  parameter(6); 
k_dec =  parameter(7); 
k_li =  parameter(8); 
k_p =  parameter(9); 
k_pr =  parameter(10); 
p_h2o =  parameter(11); 
 
% Define algebraic equations 
 
p_ch4 = x(9)*R*T/16; 
p_co2 = x(10)*R*T/44; 
p_gas = p_ch4 + p_co2 + p_h2o ; 
q_gas = k_p * ( p_gas - p_atm) *p_gas/p_atm; 
 
% Define output (states) 
 
for i = 1:10 
    y(i) = x(i);  
end 
 
% Define output (algebraic components 
 
y(11) = p_ch4; 
y(12) = p_co2; 
y(13) = p_gas; 
y(14) = q_gas; 
 
% Define inhibtion functions 
 
inhibition = 0; 
 
% Define rate equations 
 
rate(1) = k_ch * x(5); 
rate(2) = k_pr * x(6); 
rate(3) = k_li * x(7); 
rate(4) = k_dec * x(8); 
rate(5) = k_La*(x(1) -16*(K_H_ch4*p_ch4)); 
rate(6) = k_La*(x(2) - 44*(K_H_co2*p_co2)); 
