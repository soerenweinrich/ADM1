function dx = ADM1_R1_mass(t,x,s,input,parameter) 

% -------------------------------------------------------------------------
%
% Simplified (mass-based)
% Anaerobic Digestion Model No. 1 (ADM1)
%
% ADM1-R1
%
% Implementation in Matlab 
% R2019b (The MathWorks, Inc.)
%
% Matlab ODE function 
%
% Version 1.1
%
% https://github.com/soerenweinrich/ADM1
%
% Copyright (c) 2021 Sören Weinrich
% E-Mail: soeren.weinrich@dbfz.de
%
% Additional information (citation example):
%
% Weinrich, S.; Nelles, M. (2021).
% Systematic simplification of the Anaerobic Digestion Model No. 1 (ADM1) -
% Model development and stoichiometric analysis. Bioresource Technology. 
% In press. https://doi.org/10.1016/j.biortech.2021.125124.
%
% -------------------------------------------------------------------------

% Interpolate input parameters 
 
[x_in,q_in]= ADM1_interp_input(input,t); 
 
% System parameters  
 
V_liq = s(1); 
V_gas = s(2); 
p_atm = s(3); 
 
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
K_aa =  parameter(11); 
K_ac =  parameter(12); 
K_bu =  parameter(13); 
K_fa =  parameter(14); 
K_pro =  parameter(15); 
K_su =  parameter(16); 
K_va =  parameter(17); 
K_w =  parameter(18); 
R =  parameter(19); 
T =  parameter(20); 
k_AB_IN =  parameter(21); 
k_AB_ac =  parameter(22); 
k_AB_bu =  parameter(23); 
k_AB_co2 =  parameter(24); 
k_AB_pro =  parameter(25); 
k_AB_va =  parameter(26); 
k_La =  parameter(27); 
k_ch =  parameter(28); 
k_dec =  parameter(29); 
k_li =  parameter(30); 
k_m_aa =  parameter(31); 
k_m_ac =  parameter(32); 
k_m_bu =  parameter(33); 
k_m_fa =  parameter(34); 
k_m_pro =  parameter(35); 
k_m_su =  parameter(36); 
k_m_va =  parameter(37); 
k_p =  parameter(38); 
k_pr =  parameter(39); 
pK_l_aa =  parameter(40); 
pK_l_ac =  parameter(41); 
pK_u_aa =  parameter(42); 
pK_u_ac =  parameter(43); 
p_h2o =  parameter(44); 
 
% Define algebraic equations 
 
S_nh4_i = x(10) - x(29); 
S_co2 = x(9) - x(28); 
phi = x(22) + S_nh4_i/17 - x(28)/44 - x(27)/60 - x(26)/74 - x(25)/88 - x(24)/102 - x(23); 
S_H = -phi*0.5 + 0.5*sqrt(phi*phi+4*K_w); 
pH = -log10(S_H); 
p_ch4 = x(30)*R*T/16; 
p_co2 = x(31)*R*T/44; 
p_gas = p_ch4 + p_co2 + p_h2o ; 
q_gas = k_p * ( p_gas - p_atm) *p_gas/p_atm; 
 
% Define inhibtion functions 
 
inhibition(1) = x(10) / (x(10) + K_I_IN); 
inhibition(2) = 10^(-(3/(pK_u_aa - pK_l_aa))*(pK_l_aa+pK_u_aa)/2) / (S_H^(3/(pK_u_aa - pK_l_aa)) + 10^(-(3/(pK_u_aa - pK_l_aa))*(pK_l_aa+pK_u_aa)/2)); 
inhibition(3) = 10^(-(3/(pK_u_ac - pK_l_ac))*(pK_l_ac+pK_u_ac)/2) / (S_H^(3/(pK_u_ac - pK_l_ac)) + 10^(-(3/(pK_u_ac - pK_l_ac))*(pK_l_ac+pK_u_ac)/2)); 
inhibition(4) = K_I_nh3 / (K_I_nh3 + x(29)); 
 
% Define rate equations 
 
rate(1) = k_ch * x(12); 
rate(2) = k_pr * x(13); 
rate(3) = k_li * x(14); 
rate(4) = k_m_su * x(1) / (K_su + x(1))*x(15) * inhibition(1) * inhibition(2); 
rate(5) = k_m_aa * x(2) / (K_aa + x(2))*x(16) * inhibition(1) * inhibition(2); 
rate(6) = k_m_fa * x(3) / (K_fa + x(3))*x(17) * inhibition(1) * inhibition(2); 
rate(7) = k_m_va * x(4) / (K_va + x(4))*x(18)*x(4)/(x(5) + x(4) + 1e-8) * inhibition(1) * inhibition(2); 
rate(8) = k_m_bu * x(5) / (K_bu + x(5))*x(19)*x(5)/(x(4) + x(5) + 1e-8) * inhibition(1) * inhibition(2); 
rate(9) = k_m_pro * x(6) / (K_pro + x(6))*x(20) * inhibition(1) * inhibition(2); 
rate(10) = k_m_ac * x(7) / (K_ac + x(7))*x(21) * inhibition(1) * inhibition(3) * inhibition(4); 
rate(11) = k_dec * x(15); 
rate(12) = k_dec * x(16); 
rate(13) = k_dec * x(17); 
rate(14) = k_dec * x(18); 
rate(15) = k_dec * x(19); 
rate(16) = k_dec * x(20); 
rate(17) = k_dec * x(21); 
rate(18) = k_AB_va*(x(24)  * (K_a_va + S_H) - K_a_va*x(4)); 
rate(19) = k_AB_bu*(x(25) * (K_a_bu + S_H) - K_a_bu * x(5)); 
rate(20) = k_AB_pro*(x(26) * (K_a_pro + S_H) - K_a_pro * x(6)); 
rate(21) = k_AB_ac*(x(27) * (K_a_ac + S_H) - K_a_ac * x(7)); 
rate(22) = k_AB_co2*(x(28) * (K_a_co2 + S_H) - K_a_co2 * x(9)); 
rate(23) = k_AB_IN*(x(29) * (K_a_IN + S_H) - K_a_IN * x(10)); 
rate(24) = k_La*(x(8) -16*(K_H_ch4*p_ch4)); 
rate(25) = k_La*(S_co2 - 44*(K_H_co2*p_co2)); 
 
% Define process equations 
 
process(1) =  1.1111 * rate(1) +  0.13482 * rate(3) -12.0373 * rate(4); 
process(2) = rate(2) -11.1067 * rate(5); 
process(3) =  0.95115 * rate(3) -6.4068 * rate(6); 
process(4) =  1.764 * rate(5) -10.1452 * rate(7); 
process(5) =  0.82651 * rate(4) +  2.2363 * rate(5) -10.9274 * rate(8); 
process(6) =  2.0619 * rate(4) +  0.51657 * rate(5) +  6.9368 * rate(7) -14.4449 * rate(9); 
process(7) =  4.4418 * rate(4) +  5.8626 * rate(5) +  11.3536 * rate(6) +  5.6494 * rate(7) +  14.0023 * rate(8) +  11.2133 * rate(9) -26.5447 * rate(10); 
process(8) =  0.51689 * rate(4) +  0.22083 * rate(5) +  1.2219 * rate(6) +  0.68644 * rate(7) +  0.87904 * rate(8) +  2.1242 * rate(9) +  6.7367 * rate(10) - rate(24); 
process(9) =  -0.02933 * rate(3) +  2.4433 * rate(4) +  2.0378 * rate(5) -4.3451 * rate(6) -2.6138 * rate(7) -3.0468 * rate(8) +  1.5366 * rate(9) +  18.4808 * rate(10) - rate(25); 
process(10) =  -0.15056 * rate(4) +  2.0137 * rate(5) -0.15056 * rate(6) -0.15056 * rate(7) -0.15056 * rate(8) -0.15056 * rate(9) -0.15056 * rate(10); 
process(11) =  -0.11111 * rate(1) -0.056636 * rate(3) +  0.89752 * rate(4) -4.5451 * rate(5) -2.673 * rate(6) -1.363 * rate(7) -1.7566 * rate(8) -1.2786 * rate(9) +  0.4778 * rate(10); 
process(12) =  - rate(1) +  0.18 * rate(11) +  0.18 * rate(12) +  0.18 * rate(13) +  0.18 * rate(14) +  0.18 * rate(15) +  0.18 * rate(16) +  0.18 * rate(17); 
process(13) =  - rate(2) +  0.77 * rate(11) +  0.77 * rate(12) +  0.77 * rate(13) +  0.77 * rate(14) +  0.77 * rate(15) +  0.77 * rate(16) +  0.77 * rate(17); 
process(14) =  - rate(3) +  0.05 * rate(11) +  0.05 * rate(12) +  0.05 * rate(13) +  0.05 * rate(14) +  0.05 * rate(15) +  0.05 * rate(16) +  0.05 * rate(17); 
process(15) = rate(4) - rate(11); 
process(16) = rate(5) - rate(12); 
process(17) = rate(6) - rate(13); 
process(18) = rate(7) - rate(14); 
process(19) = rate(8) - rate(15); 
process(20) = rate(9) - rate(16); 
process(21) = rate(10) - rate(17); 
process(22) = 0; 
process(23) = 0; 
process(24) =  - rate(18); 
process(25) =  - rate(19); 
process(26) =  - rate(20); 
process(27) =  - rate(21); 
process(28) =  - rate(22); 
process(29) =  - rate(23); 
process(30) =  (V_liq/V_gas) * rate(24); 
process(31) =  (V_liq/V_gas) * rate(25); 
 
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
dx(19,1) = q_in*(x_in(19) - x(19))/V_liq + process(19); 
dx(20,1) = q_in*(x_in(20) - x(20))/V_liq + process(20); 
dx(21,1) = q_in*(x_in(21) - x(21))/V_liq + process(21); 
dx(22,1) = q_in*(x_in(22) - x(22))/V_liq + process(22); 
dx(23,1) = q_in*(x_in(23) - x(23))/V_liq + process(23); 
dx(24,1) = process(24); 
dx(25,1) = process(25); 
dx(26,1) = process(26); 
dx(27,1) = process(27); 
dx(28,1) = process(28); 
dx(29,1) = process(29); 
dx(30,1) = -x(30) * q_gas / V_gas + process(30); 
dx(31,1) = -x(31) * q_gas / V_gas + process(31); 
