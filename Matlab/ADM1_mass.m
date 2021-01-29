function dx = ADM1_mass(t,x,s,input,parameter) 

% -------------------------------------------------------------------------
%
% Complete and mass-based
% Anaerobic Digestion Model No. 1 (ADM1)
%
% ADM1
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
K_H_h2 =  parameter(3); 
K_I_IN =  parameter(4); 
K_I_c4 =  parameter(5); 
K_I_fa =  parameter(6); 
K_I_nh3 =  parameter(7); 
K_I_pro =  parameter(8); 
K_a_IN =  parameter(9); 
K_a_ac =  parameter(10); 
K_a_bu =  parameter(11); 
K_a_co2 =  parameter(12); 
K_a_pro =  parameter(13); 
K_a_va =  parameter(14); 
K_aa =  parameter(15); 
K_ac =  parameter(16); 
K_bu =  parameter(17); 
K_fa =  parameter(18); 
K_h2 =  parameter(19); 
K_pro =  parameter(20); 
K_su =  parameter(21); 
K_va =  parameter(22); 
K_w =  parameter(23); 
R =  parameter(24); 
T =  parameter(25); 
k_AB_IN =  parameter(26); 
k_AB_ac =  parameter(27); 
k_AB_bu =  parameter(28); 
k_AB_co2 =  parameter(29); 
k_AB_pro =  parameter(30); 
k_AB_va =  parameter(31); 
k_La =  parameter(32); 
k_ch =  parameter(33); 
k_dec =  parameter(34); 
k_li =  parameter(35); 
k_m_aa =  parameter(36); 
k_m_ac =  parameter(37); 
k_m_bu =  parameter(38); 
k_m_fa =  parameter(39); 
k_m_h2 =  parameter(40); 
k_m_pro =  parameter(41); 
k_m_su =  parameter(42); 
k_m_va =  parameter(43); 
k_p =  parameter(44); 
k_pr =  parameter(45); 
pK_l_aa =  parameter(46); 
pK_l_ac =  parameter(47); 
pK_l_h2 =  parameter(48); 
pK_u_aa =  parameter(49); 
pK_u_ac =  parameter(50); 
pK_u_h2 =  parameter(51); 
p_h2o =  parameter(52); 
 
% Define algebraic equations 
 
S_nh4_i = x(11) - x(31); 
S_co2 = x(10) - x(30); 
phi = x(24) + S_nh4_i/17 - x(30)/44 - x(29)/60 - x(28)/74 - x(27)/88 - x(26)/102 - x(25); 
S_H = -phi*0.5 + 0.5*sqrt(phi*phi+4*K_w); 
pH = -log10(S_H); 
p_h2 = x(32)*R*T/2; 
p_ch4 = x(33)*R*T/16; 
p_co2 = x(34)*R*T/44; 
p_gas = p_h2 + p_ch4 + p_co2 + p_h2o ; 
q_gas = k_p * ( p_gas - p_atm) *p_gas/p_atm; 
 
% Define inhibtion functions 
 
inhibition(1) = x(11) / (x(11) + K_I_IN); 
inhibition(2) = K_I_fa / (K_I_fa + x(8)); 
inhibition(3) = K_I_c4 / (K_I_c4 + x(8)); 
inhibition(4) = K_I_pro / (K_I_pro + x(8)); 
inhibition(5) = 10^(-(3/(pK_u_aa - pK_l_aa))*(pK_l_aa+pK_u_aa)/2) / (S_H^(3/(pK_u_aa - pK_l_aa)) + 10^(-(3/(pK_u_aa - pK_l_aa))*(pK_l_aa+pK_u_aa)/2)); 
inhibition(6) = 10^(-(3/(pK_u_ac - pK_l_ac))*(pK_l_ac+pK_u_ac)/2) / (S_H^(3/(pK_u_ac - pK_l_ac)) + 10^(-(3/(pK_u_ac - pK_l_ac))*(pK_l_ac+pK_u_ac)/2)); 
inhibition(7) = 10^(-(3/(pK_u_h2 - pK_l_h2))*(pK_l_h2+pK_u_h2)/2) / (S_H^(3/(pK_u_h2 - pK_l_h2)) + 10^(-(3/(pK_u_h2 - pK_l_h2))*(pK_l_h2+pK_u_h2)/2)); 
inhibition(8) = K_I_nh3 / (K_I_nh3 + x(31)); 
 
% Define rate equations 
 
rate(1) = k_ch * x(13); 
rate(2) = k_pr * x(14); 
rate(3) = k_li * x(15); 
rate(4) = k_m_su * x(1) / (K_su + x(1))*x(16) * inhibition(1) * inhibition(5); 
rate(5) = k_m_aa * x(2) / (K_aa + x(2))*x(17) * inhibition(1) * inhibition(5); 
rate(6) = k_m_fa * x(3) / (K_fa + x(3))*x(18) * inhibition(1) * inhibition(2) * inhibition(5); 
rate(7) = k_m_va * x(4) / (K_va + x(4))*x(19)*x(4)/(x(5) + x(4) + 1e-8) * inhibition(1) * inhibition(3) * inhibition(5); 
rate(8) = k_m_bu * x(5) / (K_bu + x(5))*x(20)*x(5)/(x(4) + x(5) + 1e-8) * inhibition(1) * inhibition(3) * inhibition(5); 
rate(9) = k_m_pro * x(6) / (K_pro + x(6))*x(21) * inhibition(1) * inhibition(4) * inhibition(5); 
rate(10) = k_m_ac * x(7) / (K_ac + x(7))*x(22) * inhibition(1) * inhibition(6) * inhibition(8); 
rate(11) = k_m_h2 * x(8) / (K_h2 + x(8))*x(23) * inhibition(1) * inhibition(7); 
rate(12) = k_dec * x(16); 
rate(13) = k_dec * x(17); 
rate(14) = k_dec * x(18); 
rate(15) = k_dec * x(19); 
rate(16) = k_dec * x(20); 
rate(17) = k_dec * x(21); 
rate(18) = k_dec * x(22); 
rate(19) = k_dec * x(23); 
rate(20) = k_AB_va*(x(26)  * (K_a_va + S_H) - K_a_va*x(4)); 
rate(21) = k_AB_bu*(x(27) * (K_a_bu + S_H) - K_a_bu * x(5)); 
rate(22) = k_AB_pro*(x(28) * (K_a_pro + S_H) - K_a_pro * x(6)); 
rate(23) = k_AB_ac*(x(29) * (K_a_ac + S_H) - K_a_ac * x(7)); 
rate(24) = k_AB_co2*(x(30) * (K_a_co2 + S_H) - K_a_co2 * x(10)); 
rate(25) = k_AB_IN*(x(31) * (K_a_IN + S_H) - K_a_IN * x(11)); 
rate(26) = k_La*(x(8) -2*(K_H_h2*p_h2)); 
rate(27) = k_La*(x(9) -16*(K_H_ch4*p_ch4)); 
rate(28) = k_La*(S_co2 - 44*(K_H_co2*p_co2)); 
 
% Define process equations 
 
process(1) =  1.1111 * rate(1) +  0.13482 * rate(3) -13.2724 * rate(4); 
process(2) = rate(2) -11.5665 * rate(5); 
process(3) =  0.95115 * rate(3) -8.2136 * rate(6); 
process(4) =  1.8371 * rate(5) -11.5757 * rate(7); 
process(5) =  0.91131 * rate(4) +  2.3289 * rate(5) -12.9817 * rate(8); 
process(6) =  2.2734 * rate(4) +  0.53795 * rate(5) +  7.9149 * rate(7) -23.3892 * rate(9); 
process(7) =  4.8975 * rate(4) +  6.1053 * rate(5) +  14.5554 * rate(6) +  6.4459 * rate(7) +  16.6347 * rate(8) +  18.1566 * rate(9) -26.5447 * rate(10); 
process(8) =  0.30475 * rate(4) +  0.12297 * rate(5) +  0.83761 * rate(6) +  0.41881 * rate(7) +  0.55841 * rate(8) +  1.8392 * rate(9) -2.9703 * rate(11) - rate(26); 
process(9) =  6.7367 * rate(10) +  5.5548 * rate(11) - rate(27); 
process(10) =  -0.02933 * rate(3) +  4.4571 * rate(4) +  2.8335 * rate(5) -0.72457 * rate(6) -0.55945 * rate(7) -0.38907 * rate(8) +  13.1283 * rate(9) +  18.4808 * rate(10) -17.1839 * rate(11) - rate(28); 
process(11) =  -0.15056 * rate(4) +  2.1033 * rate(5) -0.15056 * rate(6) -0.15056 * rate(7) -0.15056 * rate(8) -0.15056 * rate(9) -0.15056 * rate(10) -0.15056 * rate(11); 
process(12) =  -0.1111 * rate(1) -0.05664 * rate(3) -0.4211 * rate(4) -5.3025 * rate(5) -7.3043 * rate(6) -3.4939 * rate(7) -4.6718 * rate(8) -10.5843 * rate(9) +  0.47776 * rate(10) +  13.75 * rate(11); 
process(13) =  - rate(1) +  0.18 * rate(12) +  0.18 * rate(13) +  0.18 * rate(14) +  0.18 * rate(15) +  0.18 * rate(16) +  0.18 * rate(17) +  0.18 * rate(18) +  0.18 * rate(19); 
process(14) =  - rate(2) +  0.77 * rate(12) +  0.77 * rate(13) +  0.77 * rate(14) +  0.77 * rate(15) +  0.77 * rate(16) +  0.77 * rate(17) +  0.77 * rate(18) +  0.77 * rate(19); 
process(15) =  - rate(3) +  0.05 * rate(12) +  0.05 * rate(13) +  0.05 * rate(14) +  0.05 * rate(15) +  0.05 * rate(16) +  0.05 * rate(17) +  0.05 * rate(18) +  0.05 * rate(19); 
process(16) = rate(4) - rate(12); 
process(17) = rate(5) - rate(13); 
process(18) = rate(6) - rate(14); 
process(19) = rate(7) - rate(15); 
process(20) = rate(8) - rate(16); 
process(21) = rate(9) - rate(17); 
process(22) = rate(10) - rate(18); 
process(23) = rate(11) - rate(19); 
process(24) = 0; 
process(25) = 0; 
process(26) =  - rate(20); 
process(27) =  - rate(21); 
process(28) =  - rate(22); 
process(29) =  - rate(23); 
process(30) =  - rate(24); 
process(31) =  - rate(25); 
process(32) =  (V_liq/V_gas) * rate(26); 
process(33) =  (V_liq/V_gas) * rate(27); 
process(34) =  (V_liq/V_gas) * rate(28); 
 
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
dx(24,1) = q_in*(x_in(24) - x(24))/V_liq + process(24); 
dx(25,1) = q_in*(x_in(25) - x(25))/V_liq + process(25); 
dx(26,1) = process(26); 
dx(27,1) = process(27); 
dx(28,1) = process(28); 
dx(29,1) = process(29); 
dx(30,1) = process(30); 
dx(31,1) = process(31); 
dx(32,1) = -x(32) * q_gas / V_gas + process(32); 
dx(33,1) = -x(33) * q_gas / V_gas + process(33); 
dx(34,1) = -x(34) * q_gas / V_gas + process(34); 
