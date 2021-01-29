function [t,y,rate,inhibition] = ADM1_R1_mass_output(t,x,s,input,parameter) 

% -------------------------------------------------------------------------
%
% Simplified and mass-based
% Anaerobic Digestion Model No. 1 (ADM1)
%
% ADM1-R1
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
 
% Define output (states) 
 
for i = 1:31
    y(i) = x(i);  
end 
 
% Define output (algebraic components)
 
y(32) = S_nh4_i; 
y(33) = S_co2; 
y(34) = phi; 
y(35) = S_H; 
y(36) = pH; 
y(37) = p_ch4; 
y(38) = p_co2; 
y(39) = p_gas; 
y(40) = q_gas; 
 
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
