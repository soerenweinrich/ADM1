// ------------------------------------------------------------------------
//
// Original (mass-based)
// Anaerobic Digestion Model No. 1 (ADM1)
//
// ADM1
//
// Implementation in Matlab/Simulink 
// R2019b (The MathWorks, Inc.)
//
// Level 2 S-Function
//
// Version 1.1
//
// https://github.com/soerenweinrich/ADM1
//
// Copyright (c) 2021 Sören Weinrich
// E-Mail: soeren.weinrich@dbfz.de
//
// Implementation based on the detailed C-Mex template
// provided by The MathWorks, Inc.
//
// Additional information (citation example):
//
// Weinrich, S.; Nelles, M. (2021).
// Systematic simplification of the Anaerobic Digestion Model No. 1 (ADM1) -
// Model development and stoichiometric analysis. Bioresource Technology.
// In press. https://doi.org/10.1016/j.biortech.2021.125124.
//
// ------------------------------------------------------------------------

#define S_FUNCTION_NAME ADM1_mass
#define S_FUNCTION_LEVEL 2 

#include "simstruc.h"
#include <math.h>

#define INITIAL(S)          ssGetSFcnParam(S,0)
#define PARAMETER(S)        ssGetSFcnParam(S,1)
#define SYSTEM(S)           ssGetSFcnParam(S,2)
#define NINPUTS(S)          ssGetSFcnParam(S,3)
#define PARAMETER_INPUT(S)  ssGetSFcnParam(S,4)
#define NOUTPUTS 1     

static void mdlInitializeSizes(SimStruct *S)
{
   int  i;  
   int_T input = mxGetPr(NINPUTS(S))[0];  

   if (!ssSetNumInputPorts(S, input)) return;  

   for (i = 0; i < input; i++) {  
      ssSetInputPortWidth(S, i,35);  
      ssSetInputPortDirectFeedThrough(S, i, 1);  
      ssSetInputPortRequiredContiguous(S, i, 1);  
      ssSetInputPortSampleTime(S, i, INHERITED_SAMPLE_TIME);  
      ssSetInputPortOffsetTime(S, i, 0.0);  
      ssSetInputPortOverWritable(S, i, 0);  
   }

   if (!ssSetNumOutputPorts(S, NOUTPUTS)) return;  

   for (i = 0; i < NOUTPUTS; i++) {  
      ssSetOutputPortWidth(S, i,45);  
      ssSetOutputPortSampleTime(S, i, INHERITED_SAMPLE_TIME);  
      ssSetOutputPortOffsetTime(S, i, 0.0);  
   }

   ssSetNumContStates(    S, 34); 
   ssSetNumDiscStates(    S, 0); 
   ssSetNumSampleTimes(   S, 1); 
   ssSetNumSFcnParams(    S, 5); 
   ssSetNumRWork(         S, 0); 
   ssSetNumIWork(         S, 0); 
   ssSetNumPWork(         S, 0); 
   ssSetNumModes(         S, 0); 
   ssSetNumNonsampledZCs( S, 0);  
   
   ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE); 

   ssSetOptions(S, 0);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
   ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
   ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_INITIALIZE_CONDITIONS 
#if defined(MDL_INITIALIZE_CONDITIONS) 

static void mdlInitializeConditions(SimStruct *S)
{
   real_T *x0 = ssGetContStates(S);
   int i;

   for (i = 0; i < 34; i++) {
      x0[i] = mxGetPr(INITIAL(S))[i];
   }
}

#endif 
#undef MDL_START  

static void mdlOutputs(SimStruct *S, int_T tid)
{
   real_T *y = ssGetOutputPortRealSignal(S,0);
   real_T *x = ssGetContStates(S);

   // System parameters

   double V_liq = mxGetPr(SYSTEM(S))[0];
   double V_gas = mxGetPr(SYSTEM(S))[1];
   double p_atm = mxGetPr(SYSTEM(S))[2];
   double q_in = 0;

   // Model parameters

   double K_H_ch4 =  mxGetPr(PARAMETER(S))[0];
   double K_H_co2 =  mxGetPr(PARAMETER(S))[1];
   double K_H_h2 =  mxGetPr(PARAMETER(S))[2];
   double K_I_IN =  mxGetPr(PARAMETER(S))[3];
   double K_I_c4 =  mxGetPr(PARAMETER(S))[4];
   double K_I_fa =  mxGetPr(PARAMETER(S))[5];
   double K_I_nh3 =  mxGetPr(PARAMETER(S))[6];
   double K_I_pro =  mxGetPr(PARAMETER(S))[7];
   double K_a_IN =  mxGetPr(PARAMETER(S))[8];
   double K_a_ac =  mxGetPr(PARAMETER(S))[9];
   double K_a_bu =  mxGetPr(PARAMETER(S))[10];
   double K_a_co2 =  mxGetPr(PARAMETER(S))[11];
   double K_a_pro =  mxGetPr(PARAMETER(S))[12];
   double K_a_va =  mxGetPr(PARAMETER(S))[13];
   double K_aa =  mxGetPr(PARAMETER(S))[14];
   double K_ac =  mxGetPr(PARAMETER(S))[15];
   double K_bu =  mxGetPr(PARAMETER(S))[16];
   double K_fa =  mxGetPr(PARAMETER(S))[17];
   double K_h2 =  mxGetPr(PARAMETER(S))[18];
   double K_pro =  mxGetPr(PARAMETER(S))[19];
   double K_su =  mxGetPr(PARAMETER(S))[20];
   double K_va =  mxGetPr(PARAMETER(S))[21];
   double K_w =  mxGetPr(PARAMETER(S))[22];
   double R =  mxGetPr(PARAMETER(S))[23];
   double T =  mxGetPr(PARAMETER(S))[24];
   double k_AB_IN =  mxGetPr(PARAMETER(S))[25];
   double k_AB_ac =  mxGetPr(PARAMETER(S))[26];
   double k_AB_bu =  mxGetPr(PARAMETER(S))[27];
   double k_AB_co2 =  mxGetPr(PARAMETER(S))[28];
   double k_AB_pro =  mxGetPr(PARAMETER(S))[29];
   double k_AB_va =  mxGetPr(PARAMETER(S))[30];
   double k_La =  mxGetPr(PARAMETER(S))[31];
   double k_ch =  mxGetPr(PARAMETER(S))[32];
   double k_dec =  mxGetPr(PARAMETER(S))[33];
   double k_li =  mxGetPr(PARAMETER(S))[34];
   double k_m_aa =  mxGetPr(PARAMETER(S))[35];
   double k_m_ac =  mxGetPr(PARAMETER(S))[36];
   double k_m_bu =  mxGetPr(PARAMETER(S))[37];
   double k_m_fa =  mxGetPr(PARAMETER(S))[38];
   double k_m_h2 =  mxGetPr(PARAMETER(S))[39];
   double k_m_pro =  mxGetPr(PARAMETER(S))[40];
   double k_m_su =  mxGetPr(PARAMETER(S))[41];
   double k_m_va =  mxGetPr(PARAMETER(S))[42];
   double k_p =  mxGetPr(PARAMETER(S))[43];
   double k_pr =  mxGetPr(PARAMETER(S))[44];
   double pK_l_aa =  mxGetPr(PARAMETER(S))[45];
   double pK_l_ac =  mxGetPr(PARAMETER(S))[46];
   double pK_l_h2 =  mxGetPr(PARAMETER(S))[47];
   double pK_u_aa =  mxGetPr(PARAMETER(S))[48];
   double pK_u_ac =  mxGetPr(PARAMETER(S))[49];
   double pK_u_h2 =  mxGetPr(PARAMETER(S))[50];
   double p_h2o =  mxGetPr(PARAMETER(S))[51];

   // Algebraic components

   double S_nh4_i;
   double S_co2;
   double phi;
   double S_H;
   double pH;
   double p_h2;
   double p_ch4;
   double p_co2;
   double p_gas;
   double q_gas;

   int i, n; 

   // Define input 

   n = ssGetNumInputPorts(S);  

   for (i = 0; i < n; i++) {
      const real_T *u = ssGetInputPortRealSignal(S,i);
      q_in = q_in + u[0]; 
   }

   // Define output (states) 

   y[0] = q_in; 

   for (i = 0; i < 34; i++){
      y[i+1] = x[i];
   }

   // Define algebraic equations

   S_nh4_i = x[10] - x[30];
   S_co2 = x[9] - x[29];
   phi = x[23] + S_nh4_i/17 - x[29]/44 - x[28]/60 - x[27]/74 - x[26]/88 - x[25]/102 - x[24];
   S_H = -phi*0.5 + 0.5*sqrt(phi*phi+4*K_w);
   pH = -log10(S_H);
   p_h2 = x[31]*R*T/2;
   p_ch4 = x[32]*R*T/16;
   p_co2 = x[33]*R*T/44;
   p_gas = p_h2 + p_ch4 + p_co2 + p_h2o ;
   q_gas = k_p * ( p_gas - p_atm) *p_gas/p_atm;

   // Define output (algebraic components)

   y[35] = S_nh4_i;
   y[36] = S_co2;
   y[37] = phi;
   y[38] = S_H;
   y[39] = pH;
   y[40] = p_h2;
   y[41] = p_ch4;
   y[42] = p_co2;
   y[43] = p_gas;
   y[44] = q_gas;
}

#define MDL_DERIVATIVES 
#undef MDL_UPDATE 
#if defined(MDL_DERIVATIVES) 

static void mdlDerivatives(SimStruct *S)
{
   real_T *dx = ssGetdX(S);
   real_T *x = ssGetContStates(S);

   // System parameters

   double V_liq = mxGetPr(SYSTEM(S))[0];
   double V_gas = mxGetPr(SYSTEM(S))[1];
   double p_atm = mxGetPr(SYSTEM(S))[2];
   double q_in;

   // Model parameters

   double K_H_ch4 =  mxGetPr(PARAMETER(S))[0];
   double K_H_co2 =  mxGetPr(PARAMETER(S))[1];
   double K_H_h2 =  mxGetPr(PARAMETER(S))[2];
   double K_I_IN =  mxGetPr(PARAMETER(S))[3];
   double K_I_c4 =  mxGetPr(PARAMETER(S))[4];
   double K_I_fa =  mxGetPr(PARAMETER(S))[5];
   double K_I_nh3 =  mxGetPr(PARAMETER(S))[6];
   double K_I_pro =  mxGetPr(PARAMETER(S))[7];
   double K_a_IN =  mxGetPr(PARAMETER(S))[8];
   double K_a_ac =  mxGetPr(PARAMETER(S))[9];
   double K_a_bu =  mxGetPr(PARAMETER(S))[10];
   double K_a_co2 =  mxGetPr(PARAMETER(S))[11];
   double K_a_pro =  mxGetPr(PARAMETER(S))[12];
   double K_a_va =  mxGetPr(PARAMETER(S))[13];
   double K_aa =  mxGetPr(PARAMETER(S))[14];
   double K_ac =  mxGetPr(PARAMETER(S))[15];
   double K_bu =  mxGetPr(PARAMETER(S))[16];
   double K_fa =  mxGetPr(PARAMETER(S))[17];
   double K_h2 =  mxGetPr(PARAMETER(S))[18];
   double K_pro =  mxGetPr(PARAMETER(S))[19];
   double K_su =  mxGetPr(PARAMETER(S))[20];
   double K_va =  mxGetPr(PARAMETER(S))[21];
   double K_w =  mxGetPr(PARAMETER(S))[22];
   double R =  mxGetPr(PARAMETER(S))[23];
   double T =  mxGetPr(PARAMETER(S))[24];
   double k_AB_IN =  mxGetPr(PARAMETER(S))[25];
   double k_AB_ac =  mxGetPr(PARAMETER(S))[26];
   double k_AB_bu =  mxGetPr(PARAMETER(S))[27];
   double k_AB_co2 =  mxGetPr(PARAMETER(S))[28];
   double k_AB_pro =  mxGetPr(PARAMETER(S))[29];
   double k_AB_va =  mxGetPr(PARAMETER(S))[30];
   double k_La =  mxGetPr(PARAMETER(S))[31];
   double k_ch =  mxGetPr(PARAMETER(S))[32];
   double k_dec =  mxGetPr(PARAMETER(S))[33];
   double k_li =  mxGetPr(PARAMETER(S))[34];
   double k_m_aa =  mxGetPr(PARAMETER(S))[35];
   double k_m_ac =  mxGetPr(PARAMETER(S))[36];
   double k_m_bu =  mxGetPr(PARAMETER(S))[37];
   double k_m_fa =  mxGetPr(PARAMETER(S))[38];
   double k_m_h2 =  mxGetPr(PARAMETER(S))[39];
   double k_m_pro =  mxGetPr(PARAMETER(S))[40];
   double k_m_su =  mxGetPr(PARAMETER(S))[41];
   double k_m_va =  mxGetPr(PARAMETER(S))[42];
   double k_p =  mxGetPr(PARAMETER(S))[43];
   double k_pr =  mxGetPr(PARAMETER(S))[44];
   double pK_l_aa =  mxGetPr(PARAMETER(S))[45];
   double pK_l_ac =  mxGetPr(PARAMETER(S))[46];
   double pK_l_h2 =  mxGetPr(PARAMETER(S))[47];
   double pK_u_aa =  mxGetPr(PARAMETER(S))[48];
   double pK_u_ac =  mxGetPr(PARAMETER(S))[49];
   double pK_u_h2 =  mxGetPr(PARAMETER(S))[50];
   double p_h2o =  mxGetPr(PARAMETER(S))[51];

   // Algebraic components

   double S_nh4_i;
   double S_co2;
   double phi;
   double S_H;
   double pH;
   double p_h2;
   double p_ch4;
   double p_co2;
   double p_gas;
   double q_gas;
   double inhibition[8];
   double rate[28];
   double process[34];
   double in_out[34];

   int i, j, n; 

   // Define algebraic equations

   S_nh4_i = x[10] - x[30];
   S_co2 = x[9] - x[29];
   phi = x[23] + S_nh4_i/17 - x[29]/44 - x[28]/60 - x[27]/74 - x[26]/88 - x[25]/102 - x[24];
   S_H = -phi*0.5 + 0.5*sqrt(phi*phi+4*K_w);
   pH = -log10(S_H);
   p_h2 = x[31]*R*T/2;
   p_ch4 = x[32]*R*T/16;
   p_co2 = x[33]*R*T/44;
   p_gas = p_h2 + p_ch4 + p_co2 + p_h2o ;
   q_gas = k_p * ( p_gas - p_atm) *p_gas/p_atm;

   // Define inhibtion functions

   inhibition[0] = x[10] / (x[10] + K_I_IN);
   inhibition[1] = K_I_fa / (K_I_fa + x[7]);
   inhibition[2] = K_I_c4 / (K_I_c4 + x[7]);
   inhibition[3] = K_I_pro / (K_I_pro + x[7]);
   inhibition[4] = pow(10,(-(3/(pK_u_aa - pK_l_aa))*(pK_l_aa+pK_u_aa)/2)) / (pow(S_H,(3/(pK_u_aa - pK_l_aa))) + pow(10,(-(3/(pK_u_aa - pK_l_aa))*(pK_l_aa+pK_u_aa)/2)));
   inhibition[5] = pow(10,(-(3/(pK_u_ac - pK_l_ac))*(pK_l_ac+pK_u_ac)/2)) / (pow(S_H,(3/(pK_u_ac - pK_l_ac))) + pow(10,(-(3/(pK_u_ac - pK_l_ac))*(pK_l_ac+pK_u_ac)/2)));
   inhibition[6] = pow(10,(-(3/(pK_u_h2 - pK_l_h2))*(pK_l_h2+pK_u_h2)/2)) / (pow(S_H,(3/(pK_u_h2 - pK_l_h2))) + pow(10,(-(3/(pK_u_h2 - pK_l_h2))*(pK_l_h2+pK_u_h2)/2)));
   inhibition[7] = K_I_nh3 / (K_I_nh3 + x[30]);

   // Define rate equations

   rate[0] = k_ch * x[12];
   rate[1] = k_pr * x[13];
   rate[2] = k_li * x[14];
   rate[3] = k_m_su * x[0] / (K_su + x[0])*x[15] * inhibition[0] * inhibition[4];
   rate[4] = k_m_aa * x[1] / (K_aa + x[1])*x[16] * inhibition[0] * inhibition[4];
   rate[5] = k_m_fa * x[2] / (K_fa + x[2])*x[17] * inhibition[0] * inhibition[1] * inhibition[4];
   rate[6] = k_m_va * x[3] / (K_va + x[3])*x[18]*x[3]/(x[4] + x[3] + 1e-8) * inhibition[0] * inhibition[2] * inhibition[4];
   rate[7] = k_m_bu * x[4] / (K_bu + x[4])*x[19]*x[4]/(x[3] + x[4] + 1e-8) * inhibition[0] * inhibition[2] * inhibition[4];
   rate[8] = k_m_pro * x[5] / (K_pro + x[5])*x[20] * inhibition[0] * inhibition[3] * inhibition[4];
   rate[9] = k_m_ac * x[6] / (K_ac + x[6])*x[21] * inhibition[0] * inhibition[5] * inhibition[7];
   rate[10] = k_m_h2 * x[7] / (K_h2 + x[7])*x[22] * inhibition[0] * inhibition[6];
   rate[11] = k_dec * x[15];
   rate[12] = k_dec * x[16];
   rate[13] = k_dec * x[17];
   rate[14] = k_dec * x[18];
   rate[15] = k_dec * x[19];
   rate[16] = k_dec * x[20];
   rate[17] = k_dec * x[21];
   rate[18] = k_dec * x[22];
   rate[19] = k_AB_va*(x[25]  * (K_a_va + S_H) - K_a_va*x[3]);
   rate[20] = k_AB_bu*(x[26] * (K_a_bu + S_H) - K_a_bu * x[4]);
   rate[21] = k_AB_pro*(x[27] * (K_a_pro + S_H) - K_a_pro * x[5]);
   rate[22] = k_AB_ac*(x[28] * (K_a_ac + S_H) - K_a_ac * x[6]);
   rate[23] = k_AB_co2*(x[29] * (K_a_co2 + S_H) - K_a_co2 * x[9]);
   rate[24] = k_AB_IN*(x[30] * (K_a_IN + S_H) - K_a_IN * x[10]);
   rate[25] = k_La*(x[7] -2*(K_H_h2*p_h2));
   rate[26] = k_La*(x[8] -16*(K_H_ch4*p_ch4));
   rate[27] = k_La*(S_co2 - 44*(K_H_co2*p_co2));

   // Define process equations

   process[0] =  1.1111 * rate[0] +  0.13482 * rate[2] -13.2724 * rate[3];
   process[1] = rate[1] -11.5665 * rate[4];
   process[2] =  0.95115 * rate[2] -8.2136 * rate[5];
   process[3] =  1.8371 * rate[4] -11.5757 * rate[6];
   process[4] =  0.91131 * rate[3] +  2.3289 * rate[4] -12.9817 * rate[7];
   process[5] =  2.2734 * rate[3] +  0.53795 * rate[4] +  7.9149 * rate[6] -23.3892 * rate[8];
   process[6] =  4.8975 * rate[3] +  6.1053 * rate[4] +  14.5554 * rate[5] +  6.4459 * rate[6] +  16.6347 * rate[7] +  18.1566 * rate[8] -26.5447 * rate[9];
   process[7] =  0.30475 * rate[3] +  0.12297 * rate[4] +  0.83761 * rate[5] +  0.41881 * rate[6] +  0.55841 * rate[7] +  1.8392 * rate[8] -2.9703 * rate[10] - rate[25];
   process[8] =  6.7367 * rate[9] +  5.5548 * rate[10] - rate[26];
   process[9] =  -0.02933 * rate[2] +  4.4571 * rate[3] +  2.8335 * rate[4] -0.72457 * rate[5] -0.55945 * rate[6] -0.38907 * rate[7] +  13.1283 * rate[8] +  18.4808 * rate[9] -17.1839 * rate[10] - rate[27];
   process[10] =  -0.15056 * rate[3] +  2.1033 * rate[4] -0.15056 * rate[5] -0.15056 * rate[6] -0.15056 * rate[7] -0.15056 * rate[8] -0.15056 * rate[9] -0.15056 * rate[10];
   process[11] =  -0.1111 * rate[0] -0.05664 * rate[2] -0.4211 * rate[3] -5.3025 * rate[4] -7.3043 * rate[5] -3.4939 * rate[6] -4.6718 * rate[7] -10.5843 * rate[8] +  0.47776 * rate[9] +  13.75 * rate[10];
   process[12] =  - rate[0] +  0.18 * rate[11] +  0.18 * rate[12] +  0.18 * rate[13] +  0.18 * rate[14] +  0.18 * rate[15] +  0.18 * rate[16] +  0.18 * rate[17] +  0.18 * rate[18];
   process[13] =  - rate[1] +  0.77 * rate[11] +  0.77 * rate[12] +  0.77 * rate[13] +  0.77 * rate[14] +  0.77 * rate[15] +  0.77 * rate[16] +  0.77 * rate[17] +  0.77 * rate[18];
   process[14] =  - rate[2] +  0.05 * rate[11] +  0.05 * rate[12] +  0.05 * rate[13] +  0.05 * rate[14] +  0.05 * rate[15] +  0.05 * rate[16] +  0.05 * rate[17] +  0.05 * rate[18];
   process[15] = rate[3] - rate[11];
   process[16] = rate[4] - rate[12];
   process[17] = rate[5] - rate[13];
   process[18] = rate[6] - rate[14];
   process[19] = rate[7] - rate[15];
   process[20] = rate[8] - rate[16];
   process[21] = rate[9] - rate[17];
   process[22] = rate[10] - rate[18];
   process[23] = 0;
   process[24] = 0;
   process[25] =  - rate[19];
   process[26] =  - rate[20];
   process[27] =  - rate[21];
   process[28] =  - rate[22];
   process[29] =  - rate[23];
   process[30] =  - rate[24];
   process[31] =  (V_liq/V_gas) * rate[25];
   process[32] =  (V_liq/V_gas) * rate[26];
   process[33] =  (V_liq/V_gas) * rate[27];

   // Define input equations

   n = ssGetNumInputPorts(S); 

   for (i = 0; i < 34; i++) {
      in_out[i] = 0;  
   } 

   for (j = 0; j < n; j++) {
      const real_T *u = ssGetInputPortRealSignal(S,j); 
      q_in = u[0]; 

      in_out[0] = in_out[0] + q_in*(u[1] - x[0])/V_liq;
      in_out[1] = in_out[1] + q_in*(u[2] - x[1])/V_liq;
      in_out[2] = in_out[2] + q_in*(u[3] - x[2])/V_liq;
      in_out[3] = in_out[3] + q_in*(u[4] - x[3])/V_liq;
      in_out[4] = in_out[4] + q_in*(u[5] - x[4])/V_liq;
      in_out[5] = in_out[5] + q_in*(u[6] - x[5])/V_liq;
      in_out[6] = in_out[6] + q_in*(u[7] - x[6])/V_liq;
      in_out[7] = in_out[7] + q_in*(u[8] - x[7])/V_liq;
      in_out[8] = in_out[8] + q_in*(u[9] - x[8])/V_liq;
      in_out[9] = in_out[9] + q_in*(u[10] - x[9])/V_liq;
      in_out[10] = in_out[10] + q_in*(u[11] - x[10])/V_liq;
      in_out[11] = in_out[11] + q_in*(u[12] - x[11])/V_liq;
      in_out[12] = in_out[12] + q_in*(u[13] - x[12])/V_liq;
      in_out[13] = in_out[13] + q_in*(u[14] - x[13])/V_liq;
      in_out[14] = in_out[14] + q_in*(u[15] - x[14])/V_liq;
      in_out[15] = in_out[15] + q_in*(u[16] - x[15])/V_liq;
      in_out[16] = in_out[16] + q_in*(u[17] - x[16])/V_liq;
      in_out[17] = in_out[17] + q_in*(u[18] - x[17])/V_liq;
      in_out[18] = in_out[18] + q_in*(u[19] - x[18])/V_liq;
      in_out[19] = in_out[19] + q_in*(u[20] - x[19])/V_liq;
      in_out[20] = in_out[20] + q_in*(u[21] - x[20])/V_liq;
      in_out[21] = in_out[21] + q_in*(u[22] - x[21])/V_liq;
      in_out[22] = in_out[22] + q_in*(u[23] - x[22])/V_liq;
      in_out[23] = in_out[23] + q_in*(u[24] - x[23])/V_liq;
      in_out[24] = in_out[24] + q_in*(u[25] - x[24])/V_liq;
      in_out[25] = in_out[25] ;
      in_out[26] = in_out[26] ;
      in_out[27] = in_out[27] ;
      in_out[28] = in_out[28] ;
      in_out[29] = in_out[29] ;
      in_out[30] = in_out[30] ;
      in_out[31] = in_out[31] -x[31] * q_gas / V_gas;
      in_out[32] = in_out[32] -x[32] * q_gas / V_gas;
      in_out[33] = in_out[33] -x[33] * q_gas / V_gas;
   }

   // Define differential equations

   dx[0] = in_out[0] + process[0];
   dx[1] = in_out[1] + process[1];
   dx[2] = in_out[2] + process[2];
   dx[3] = in_out[3] + process[3];
   dx[4] = in_out[4] + process[4];
   dx[5] = in_out[5] + process[5];
   dx[6] = in_out[6] + process[6];
   dx[7] = in_out[7] + process[7];
   dx[8] = in_out[8] + process[8];
   dx[9] = in_out[9] + process[9];
   dx[10] = in_out[10] + process[10];
   dx[11] = in_out[11] + process[11];
   dx[12] = in_out[12] + process[12];
   dx[13] = in_out[13] + process[13];
   dx[14] = in_out[14] + process[14];
   dx[15] = in_out[15] + process[15];
   dx[16] = in_out[16] + process[16];
   dx[17] = in_out[17] + process[17];
   dx[18] = in_out[18] + process[18];
   dx[19] = in_out[19] + process[19];
   dx[20] = in_out[20] + process[20];
   dx[21] = in_out[21] + process[21];
   dx[22] = in_out[22] + process[22];
   dx[23] = in_out[23] + process[23];
   dx[24] = in_out[24] + process[24];
   dx[25] = in_out[25] + process[25];
   dx[26] = in_out[26] + process[26];
   dx[27] = in_out[27] + process[27];
   dx[28] = in_out[28] + process[28];
   dx[29] = in_out[29] + process[29];
   dx[30] = in_out[30] + process[30];
   dx[31] =  -x[31] * q_gas / V_gas + process[31];
   dx[32] =  -x[32] * q_gas / V_gas + process[32];
   dx[33] =  -x[33] * q_gas / V_gas + process[33];
}

#endif 

static void mdlTerminate(SimStruct *S)
{
}

#ifdef	MATLAB_MEX_FILE    
#include "simulink.c"      
#else
#include "cg_sfun.h"      
#endif