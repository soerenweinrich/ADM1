function dx = ADM1_R2_mass(t,x,s,input,parameter) 

% -------------------------------------------------------------------------
%
% Simplified and mass-based
% Anaerobic Digestion Model No. 1 (ADM1)
%
% ADM1-R2
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
K_a_bu =  parameter(7); 
K_a_co2 =  parameter(8); 
K_a_pro =  parameter(9); 
K_a_va =  parameter(10); 
K_ac =  parameter(11); 
K_bu =  parameter(12); 
K_pro =  parameter(13); 
K_va =  parameter(14); 
K_w =  parameter(15); 
R =  parameter(16); 
T =  parameter(17); 
k_AB_IN =  parameter(18); 
k_AB_ac =  parameter(19); 
k_AB_bu =  parameter(20); 
k_AB_co2 =  parameter(21); 
k_AB_pro =  parameter(22); 
k_AB_va =  parameter(23); 
k_La =  parameter(24); 
k_ch =  parameter(25); 
k_dec =  parameter(26); 
k_li =  parameter(27); 
k_m_ac =  parameter(28); 
k_m_bu =  parameter(29); 
k_m_pro =  parameter(30); 
k_m_va =  parameter(31); 
k_p =  parameter(32); 
k_pr =  parameter(33); 
pK_l_aa =  parameter(34); 
pK_l_ac =  parameter(35); 
pK_u_aa =  parameter(36); 
pK_u_ac =  parameter(37); 
p_h2o =  parameter(38); 
 
% Define algebraic equations 
 
S_nh4_i = x(7) - x(24); 
S_co2 = x(6) - x(23); 
phi = x(17) + S_nh4_i/17 - x(23)/44 - x(22)/60 - x(21)/74 - x(20)/88 - x(19)/102 - x(18); 
S_H = -phi*0.5 + 0.5*sqrt(phi*phi+4*K_w); 
pH = -log10(S_H); 
p_ch4 = x(25)*R*T/16; 
p_co2 = x(26)*R*T/44; 
p_gas = p_ch4 + p_co2 + p_h2o ; 
q_gas = k_p * ( p_gas - p_atm) *p_gas/p_atm; 
 
% Define inhibtion functions 
 
inhibition(1) = x(7) / (x(7) + K_I_IN); 
inhibition(2) = 10^(-(3/(pK_u_aa - pK_l_aa))*(pK_l_aa+pK_u_aa)/2) / (S_H^(3/(pK_u_aa - pK_l_aa)) + 10^(-(3/(pK_u_aa - pK_l_aa))*(pK_l_aa+pK_u_aa)/2)); 
inhibition(3) = 10^(-(3/(pK_u_ac - pK_l_ac))*(pK_l_ac+pK_u_ac)/2) / (S_H^(3/(pK_u_ac - pK_l_ac)) + 10^(-(3/(pK_u_ac - pK_l_ac))*(pK_l_ac+pK_u_ac)/2)); 
inhibition(4) = K_I_nh3 / (K_I_nh3 + x(24)); 
 
% Define rate equations 
 
rate(1) = k_ch * x(9); 
rate(2) = k_pr * x(10); 
rate(3) = k_li * x(11); 
rate(4) = k_m_va * x(1) / (K_va + x(1))*x(13)*x(1)/(x(2) + x(1) + 1e-8) * inhibition(1) * inhibition(2); 
rate(5) = k_m_bu * x(2) / (K_bu + x(2))*x(14)*x(2)/(x(1) + x(2) + 1e-8) * inhibition(1) * inhibition(2); 
rate(6) = k_m_pro * x(3) / (K_pro + x(3))*x(15) * inhibition(1) * inhibition(2); 
rate(7) = k_m_ac * x(4) / (K_ac + x(4))*x(16) * inhibition(1) * inhibition(3) * inhibition(4); 
rate(8) = k_dec * x(12); 
rate(9) = k_dec * x(13); 
rate(10) = k_dec * x(14); 
rate(11) = k_dec * x(15); 
rate(12) = k_dec * x(16); 
rate(13) = k_AB_va*(x(19)  * (K_a_va + S_H) - K_a_va*x(1)); 
rate(14) = k_AB_bu*(x(20) * (K_a_bu + S_H) - K_a_bu * x(2)); 
rate(15) = k_AB_pro*(x(21) * (K_a_pro + S_H) - K_a_pro * x(3)); 
rate(16) = k_AB_ac*(x(22) * (K_a_ac + S_H) - K_a_ac * x(4)); 
rate(17) = k_AB_co2*(x(23) * (K_a_co2 + S_H) - K_a_co2 * x(6)); 
rate(18) = k_AB_IN*(x(24) * (K_a_IN + S_H) - K_a_IN * x(7)); 
rate(19) = k_La*(x(5) -16*(K_H_ch4*p_ch4)); 
rate(20) = k_La*(S_co2 - 44*(K_H_co2*p_co2)); 
 
% Define process equations 
 
process(1) =  0.15883 * rate(2) -10.1452 * rate(4); 
process(2) =  0.076292 * rate(1) +  0.20135 * rate(2) +  0.0092572 * rate(3) -10.9274 * rate(5); 
process(3) =  0.19032 * rate(1) +  0.046509 * rate(2) +  0.023094 * rate(3) +  6.9368 * rate(4) -14.4449 * rate(6); 
process(4) =  0.41 * rate(1) +  0.52784 * rate(2) +  1.7353 * rate(3) +  5.6494 * rate(4) +  14.0023 * rate(5) +  11.2133 * rate(6) -26.5447 * rate(7); 
process(5) =  0.047712 * rate(1) +  0.019882 * rate(2) +  0.18719 * rate(3) +  0.68644 * rate(4) +  0.87904 * rate(5) +  2.1242 * rate(6) +  6.7367 * rate(7) - rate(19); 
process(6) =  0.22553 * rate(1) +  0.18347 * rate(2) -0.64703 * rate(3) -2.6138 * rate(4) -3.0468 * rate(5) +  1.5366 * rate(6) +  18.4808 * rate(7) - rate(20); 
process(7) =  -0.013897 * rate(1) +  0.18131 * rate(2) -0.024038 * rate(3) -0.15056 * rate(4) -0.15056 * rate(5) -0.15056 * rate(6) -0.15056 * rate(7); 
process(8) =  -0.028264 * rate(1) -0.40923 * rate(2) -0.44342 * rate(3) -1.363 * rate(4) -1.7566 * rate(5) -1.2786 * rate(6) +  0.4778 * rate(7); 
process(9) =  - rate(1) +  0.18 * rate(8) +  0.18 * rate(9) +  0.18 * rate(10) +  0.18 * rate(11) +  0.18 * rate(12); 
process(10) =  - rate(2) +  0.77 * rate(8) +  0.77 * rate(9) +  0.77 * rate(10) +  0.77 * rate(11) +  0.77 * rate(12); 
process(11) =  - rate(3) +  0.05 * rate(8) +  0.05 * rate(9) +  0.05 * rate(10) +  0.05 * rate(11) +  0.05 * rate(12); 
process(12) =  0.092305 * rate(1) +  0.090036 * rate(2) +  0.15966 * rate(3) - rate(8); 
process(13) = rate(4) - rate(9); 
process(14) = rate(5) - rate(10); 
process(15) = rate(6) - rate(11); 
process(16) = rate(7) - rate(12); 
process(17) = 0; 
process(18) = 0; 
process(19) =  - rate(13); 
process(20) =  - rate(14); 
process(21) =  - rate(15); 
process(22) =  - rate(16); 
process(23) =  - rate(17); 
process(24) =  - rate(18); 
process(25) =  (V_liq/V_gas) * rate(19); 
process(26) =  (V_liq/V_gas) * rate(20); 
 
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
dx(13,1) = q_in*(x_in(13) - x(13))/V_liq + process(13); 
dx(14,1) = q_in*(x_in(14) - x(14))/V_liq + process(14); 
dx(15,1) = q_in*(x_in(15) - x(15))/V_liq + process(15); 
dx(16,1) = q_in*(x_in(16) - x(16))/V_liq + process(16); 
dx(17,1) = q_in*(x_in(17) - x(17))/V_liq + process(17); 
dx(18,1) = q_in*(x_in(18) - x(18))/V_liq + process(18); 
dx(19,1) = process(19); 
dx(20,1) = process(20); 
dx(21,1) = process(21); 
dx(22,1) = process(22); 
dx(23,1) = process(23); 
dx(24,1) = process(24); 
dx(25,1) = -x(25) * q_gas / V_gas + process(25); 
dx(26,1) = -x(26) * q_gas / V_gas + process(26); 
