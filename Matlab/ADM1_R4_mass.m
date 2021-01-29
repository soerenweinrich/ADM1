function dx = ADM1_R4_mass(t,x,s,input,parameter) 

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

% Interpolate input parameters 
 
[x_in,T_op,q_in]= input_function(input,t); 
 
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
 
% Define inhibtion functions 
 
inhibition = 0; 
 
% Define rate equations 
 
rate(1) = k_ch * x(5); 
rate(2) = k_pr * x(6); 
rate(3) = k_li * x(7); 
rate(4) = k_dec * x(8); 
rate(5) = k_La*(x(1) -16*(K_H_ch4*p_ch4)); 
rate(6) = k_La*(x(2) - 44*(K_H_co2*p_co2)); 
 
% Define process equations 
 
process(1) =  0.24819 * rate(1) +  0.32208 * rate(2) +  0.63928 * rate(3) - rate(5); 
process(2) =  0.68087 * rate(1) +  0.79543 * rate(2) +  0.58172 * rate(3) - rate(6); 
process(3) =  -0.02065 * rate(1) +  0.16892 * rate(2) -0.034418 * rate(3); 
process(4) =  -0.045576 * rate(1) -0.45876 * rate(2) -0.41518 * rate(3); 
process(5) =  - rate(1) +  0.18 * rate(4); 
process(6) =  - rate(2) +  0.77 * rate(4); 
process(7) =  - rate(3) +  0.05 * rate(4); 
process(8) =  0.13716 * rate(1) +  0.17233 * rate(2) +  0.2286 * rate(3) - rate(4); 
process(9) =  (V_liq/V_gas) * rate(5); 
process(10) =  (V_liq/V_gas) * rate(6); 
 
% Define differential equations 
 
dx(1,1) = q_in*(x_in(1) - x(1))/V_liq + process(1); 
dx(2,1) = q_in*(x_in(2) - x(2))/V_liq + process(2); 
dx(3,1) = q_in*(x_in(3) - x(3))/V_liq + process(3); 
dx(4,1) = q_in*(x_in(4) - x(4))/V_liq + process(4); 
dx(5,1) = q_in*(x_in(5) -x(5))/V_liq + process(5); 
dx(6,1) = q_in*(x_in(6) - x(6))/V_liq + process(6); 
dx(7,1) = q_in*(x_in(7) - x(7))/V_liq + process(7); 
dx(8,1) = q_in*(x_in(8) - x(8))/V_liq + process(8); 
dx(9,1) = - x(9) * q_gas / V_gas + process(9); 
dx(10,1) = - x(10) * q_gas / V_gas + process(10); 
