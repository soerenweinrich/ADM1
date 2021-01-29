function dx = ADM1_R3_mass(t,x,s,input,parameter) 

% -------------------------------------------------------------------------
%
% Simplified and mass-based
% Anaerobic Digestion Model No. 1 (ADM1)
%
% ADM1-R3
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
K_I_IN =  parameter(3); 
K_I_nh3 =  parameter(4); 
K_a_IN =  parameter(5); 
K_a_ac =  parameter(6); 
K_a_co2 =  parameter(7); 
K_ac =  parameter(8); 
K_w =  parameter(9); 
R =  parameter(10); 
T =  parameter(11); 
k_AB_IN =  parameter(12); 
k_AB_ac =  parameter(13); 
k_AB_co2 =  parameter(14); 
k_La =  parameter(15); 
k_ch =  parameter(16); 
k_dec =  parameter(17); 
k_li =  parameter(18); 
k_m_ac =  parameter(19); 
k_p =  parameter(20); 
k_pr =  parameter(21); 
pK_l_ac =  parameter(22); 
pK_u_ac =  parameter(23); 
p_h2o =  parameter(24); 
 
% Define algebraic equations 
 
S_nh4_i = x(4) - x(15); 
S_co2 = x(3) - x(14); 
phi = x(11) + S_nh4_i/17 - x(14)/44 - x(13)/60  - x(12); 
S_H = -phi*0.5 + 0.5*sqrt(phi*phi+4*K_w); 
pH = -log10(S_H); 
p_ch4 = x(16)*R*T/16; 
p_co2 = x(17)*R*T/44; 
p_gas = p_ch4 + p_co2 + p_h2o ; 
q_gas = k_p * ( p_gas - p_atm) *p_gas/p_atm; 
 
% Define inhibtion functions 
 
inhibition(1) = x(4) / (x(4) + K_I_IN); 
inhibition(2) = 10^(-(3/(pK_u_ac - pK_l_ac))*(pK_l_ac+pK_u_ac)/2) / (S_H^(3/(pK_u_ac - pK_l_ac)) + 10^(-(3/(pK_u_ac - pK_l_ac))*(pK_l_ac+pK_u_ac)/2)); 
inhibition(3) = K_I_nh3 / (K_I_nh3 + x(15)); 
 
% Define rate equations 
 
rate(1) = k_ch * x(6); 
rate(2) = k_pr * x(7); 
rate(3) = k_li * x(8); 
rate(4) = k_m_ac * x(1) / (K_ac + x(1))*x(10) * inhibition(1) * inhibition(2) * inhibition(3); 
rate(5) = k_dec * x(9); 
rate(6) = k_dec * x(10); 
rate(7) = k_AB_ac*(x(13) * (K_a_ac + S_H) - K_a_ac * x(1)); 
rate(8) = k_AB_co2*(x(14) * (K_a_co2 + S_H) - K_a_co2 * x(3)); 
rate(9) = k_AB_IN*(x(15) * (K_a_IN + S_H) - K_a_IN * x(4)); 
rate(10) = k_La*(x(2) -16*(K_H_ch4*p_ch4)); 
rate(11) = k_La*(S_co2 - 44*(K_H_co2*p_co2)); 
 
% Define process equations 
 
process(1) =  0.6555 * rate(1) +  0.9947 * rate(2) +  1.7651 * rate(3) -26.5447 * rate(4); 
process(2) =  0.081837 * rate(1) +  0.069636 * rate(2) +  0.19133 * rate(3) +  6.7367 * rate(4) - rate(10); 
process(3) =  0.2245 * rate(1) +  0.10291 * rate(2) -0.64716 * rate(3) +  18.4808 * rate(4) - rate(11); 
process(4) =  -0.016932 * rate(1) +  0.17456 * rate(2) -0.024406 * rate(3) -0.15056 * rate(4); 
process(5) =  -0.057375 * rate(1) -0.47666 * rate(2) -0.44695 * rate(3) +  0.4778 * rate(4); 
process(6) =  - rate(1) +  0.18 * rate(5) +  0.18 * rate(6); 
process(7) =  - rate(2) +  0.77 * rate(5) +  0.77 * rate(6); 
process(8) =  - rate(3) +  0.05 * rate(5) +  0.05 * rate(6); 
process(9) =  0.11246 * rate(1) +  0.13486 * rate(2) +  0.1621 * rate(3) - rate(5); 
process(10) = rate(4) - rate(6); 
process(11) = 0; 
process(12) = 0; 
process(13) =  - rate(7); 
process(14) =  - rate(8); 
process(15) =  - rate(9); 
process(16) =  (V_liq/V_gas) * rate(10); 
process(17) =  (V_liq/V_gas) * rate(11); 
 
% Define differential equations 
 
dx(1,1) = q_in*(x_in(1) - x(1))/V_liq + process(1); 
dx(2,1) = q_in*(x_in(2) - x(2))/V_liq + process(2); 
dx(3,1) = q_in*(x_in(3) - x(3))/V_liq + process(3); 
dx(4,1) = q_in*(x_in(4) - x(4))/V_liq + process(4); 
dx(5,1) = q_in*(x_in(5) - x(5))/V_liq + process(5); 
dx(6,1) = q_in*(x_in(6) - x(6))/V_liq + process(6); 
dx(7,1) = q_in*(x_in(7) - x(7))/V_liq + process(7); 
dx(8,1) = q_in*(x_in(8) - x(8))/V_liq + process(8); 
dx(9,1) = q_in*(x_in(9) - x(9))/V_liq + process(9); 
dx(10,1) = q_in*(x_in(10) - x(10))/V_liq + process(10); 
dx(11,1) = q_in*(x_in(11) - x(11))/V_liq + process(11); 
dx(12,1) = q_in*(x_in(12) - x(12))/V_liq + process(12); 
dx(13,1) = process(13); 
dx(14,1) = process(14); 
dx(15,1) = process(15); 
dx(16,1) = - x(16) * q_gas / V_gas + process(16); 
dx(17,1) = - x(17) * q_gas / V_gas + process(17); 
