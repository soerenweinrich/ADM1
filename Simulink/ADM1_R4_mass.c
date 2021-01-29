// Simplified and mass-based
// Anaerobic Digestion Model No. 1 (ADM1)
//
// ADM1-R4
//
// Implementation in Matlab/Simulink 
// R2019b (The MathWorks, Inc.)
//
// Level 2 S-Function
//
// Version 1.0
//
// https://github.com/soerenweinrich/ADM1
//
// Copyright (c) 2021 Sören Weinrich
// E-Mail: soeren.weinrich@dbfz.de
//
// Implementation based on the detailed C-Mex template
// provided by The MathWorks, Inc.
//
// Citation example:
//
// Weinrich, S.; Nelles, M. (2021).
// Systematic simplification of the Anaerobic Digestion Model No. 1 (ADM1) -
// Model development and stoichiometric analysis.
// Submitted to Bioresource Technology.

#define S_FUNCTION_NAME ADM1_R4_mass
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
      ssSetInputPortWidth(S, i,11);  
      ssSetInputPortDirectFeedThrough(S, i, 1);  
      ssSetInputPortRequiredContiguous(S, i, 1);  
      ssSetInputPortSampleTime(S, i, INHERITED_SAMPLE_TIME);  
      ssSetInputPortOffsetTime(S, i, 0.0);  
      ssSetInputPortOverWritable(S, i, 0);  
   }

   if (!ssSetNumOutputPorts(S, NOUTPUTS)) return;  

   for (i = 0; i < NOUTPUTS; i++) {  
      ssSetOutputPortWidth(S, i,15);  
      ssSetOutputPortSampleTime(S, i, INHERITED_SAMPLE_TIME);  
      ssSetOutputPortOffsetTime(S, i, 0.0);  
   }

   ssSetNumContStates(    S, 10); 
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

   for (i = 0; i < 10; i++) {
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
   double R =  mxGetPr(PARAMETER(S))[2];
   double T =  mxGetPr(PARAMETER(S))[3];
   double k_La =  mxGetPr(PARAMETER(S))[4];
   double k_ch =  mxGetPr(PARAMETER(S))[5];
   double k_dec =  mxGetPr(PARAMETER(S))[6];
   double k_li =  mxGetPr(PARAMETER(S))[7];
   double k_p =  mxGetPr(PARAMETER(S))[8];
   double k_pr =  mxGetPr(PARAMETER(S))[9];
   double p_h2o =  mxGetPr(PARAMETER(S))[10];

   // Algebraic components

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

   for (i = 0; i < 10; i++){
      y[i+1] = x[i];
   }

   // Define algebraic equations

   p_ch4 = x[8]*R*T/16;
   p_co2 = x[9]*R*T/44;
   p_gas = p_ch4 + p_co2 + p_h2o ;
   q_gas = k_p * ( p_gas - p_atm) *p_gas/p_atm;

   // Define output (algebraic components)

   y[11] = p_ch4;
   y[12] = p_co2;
   y[13] = p_gas;
   y[14] = q_gas;
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
   double R =  mxGetPr(PARAMETER(S))[2];
   double T =  mxGetPr(PARAMETER(S))[3];
   double k_La =  mxGetPr(PARAMETER(S))[4];
   double k_ch =  mxGetPr(PARAMETER(S))[5];
   double k_dec =  mxGetPr(PARAMETER(S))[6];
   double k_li =  mxGetPr(PARAMETER(S))[7];
   double k_p =  mxGetPr(PARAMETER(S))[8];
   double k_pr =  mxGetPr(PARAMETER(S))[9];
   double p_h2o =  mxGetPr(PARAMETER(S))[10];

   // Algebraic components

   double p_ch4;
   double p_co2;
   double p_gas;
   double q_gas;
   double rate[6];
   double process[10];
   double in_out[10];

   int i, j, n; 

   // Define algebraic equations

   p_ch4 = x[8]*R*T/16;
   p_co2 = x[9]*R*T/44;
   p_gas = p_ch4 + p_co2 + p_h2o ;
   q_gas = k_p * ( p_gas - p_atm) *p_gas/p_atm;

   // Define rate equations

   rate[0] = k_ch * x[4];
   rate[1] = k_pr * x[5];
   rate[2] = k_li * x[6];
   rate[3] = k_dec * x[7];
   rate[4] = k_La*(x[0] -16*(K_H_ch4*p_ch4));
   rate[5] = k_La*(x[1] - 44*(K_H_co2*p_co2));

   // Define process equations

   process[0] =  0.24819 * rate[0] +  0.32208 * rate[1] +  0.63928 * rate[2] - rate[4];
   process[1] =  0.68087 * rate[0] +  0.79543 * rate[1] +  0.58172 * rate[2] - rate[5];
   process[2] =  -0.02065 * rate[0] +  0.16892 * rate[1] -0.034418 * rate[2];
   process[3] =  -0.045576 * rate[0] -0.45876 * rate[1] -0.41518 * rate[2];
   process[4] =  - rate[0] +  0.18 * rate[3];
   process[5] =  - rate[1] +  0.77 * rate[3];
   process[6] =  - rate[2] +  0.05 * rate[3];
   process[7] =  0.13716 * rate[0] +  0.17233 * rate[1] +  0.2286 * rate[2] - rate[3];
   process[8] =  (V_liq/V_gas) * rate[4];
   process[9] =  (V_liq/V_gas) * rate[5];

   // Define input equations

   n = ssGetNumInputPorts(S); 

   for (i = 0; i < 10; i++) {
      in_out[i] = 0;  
   } 

   for (j = 0; j < n; j++) {
      const real_T *u = ssGetInputPortRealSignal(S,j); 
      q_in = u[0]; 

      in_out[0] = in_out[0] + q_in*(u[1] - x[0])/V_liq;
      in_out[1] = in_out[1] + q_in*(u[2] - x[1])/V_liq;
      in_out[2] = in_out[2] + q_in*(u[3] - x[2])/V_liq;
      in_out[3] = in_out[3] + q_in*(u[4] - x[3])/V_liq;
      in_out[4] = in_out[4] + q_in*(u[5] -x[4])/V_liq;
      in_out[5] = in_out[5] + q_in*(u[6] - x[5])/V_liq;
      in_out[6] = in_out[6] + q_in*(u[7] - x[6])/V_liq;
      in_out[7] = in_out[7] + q_in*(u[8] - x[7])/V_liq;
      in_out[8] = in_out[8] - x[8] * q_gas / V_gas;
      in_out[9] = in_out[9] - x[9] * q_gas / V_gas;
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
   dx[8] =  - x[8] * q_gas / V_gas + process[8];
   dx[9] =  - x[9] * q_gas / V_gas + process[9];
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