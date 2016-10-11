/*
*
* CDDL HEADER START
*
* The contents of this file are subject to the terms of the Common Development
* and Distribution License Version 1.0 (the "License").
*
* You can obtain a copy of the license at
* http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
* specific language governing permissions and limitations under the License.
*
* When distributing Covered Code, include this CDDL HEADER in each file and
* include the License file in a prominent location with the name LICENSE.CDDL.
* If applicable, add the following below this CDDL HEADER, with the fields
* enclosed by brackets "[]" replaced with your own identifying information:
*
* Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
*
* CDDL HEADER END
*

*
* Copyright (c) 2012, Regents of the University of Minnesota.  All rights reserved.
*
* Contributors:
*    Amit Singh
*/

/*******************************************************************************
*
*  Four_Body_Mistriotis_Flytzanis_Farantos 
*
*  Mistriotis-Flytzanis-Farantos Four-Body potential KIM Model Driver
*
*  Language: C
*
*  Release: 
*
*******************************************************************************/
/*******************************************************************************

* For four-body potential any potential energy function for N-particle system can be written in the
* following form of two-body, three-body and four-body terms:
     F(1, 2, ...., N)   = sum_{i; 1 <= i \neq j <= N} 0.5*phi_two(r; r = r_ij) 
                          + sum_{i; 1 <= i \neq j < k <= N} phi_three(r_ij, r_ik, r_jk);
                          + sum_{i; 1 <= i \neq j < k < l <= N} phi_four(r_ij, r_ik, r_il, r_jk, r_jl, r_kl);
                          
* For modified Stillinger-Weber proposed by "Mistriotis et al, Physical Review B 39, 1212 (1989)", 
* the two-body term phi_two can be written as:
 
     phi_two(r) = epsilon * f_2(r_cap; r_cap = r/sigma); where f_2 is 
 
     f_2(r_cap) = A * ( B*r_cap^(-p) - r_cap^-q ) * exp(1/(r_cap - a)); when  r_cap <  a, where a = cutoff/sigma,  
                = 0   when r_cap >= a

* And three-body term phi_three can be written as:
     
     phi_three(r_ij, r_ik, r_jk) = epsilon * f_3(r1_cap, r2_cap, r3_cap); where
     
     r1_cap = r_ij/sigma, r2_cap = r_ik/sigma, r3_cap = r_jk/sigma,      and

     f_3(r1_cap, r2_cap, r3_cap) = lambda * (exp(gamma((1/(r1_cap - a)) + (1/(r2_cap - a))))) 
                                 * {1 - exp[-Q((costheta_jik - costheta_0)^2)]};  when r1_cap < a && r2_cap < a
                                 = 0;                                             otherwise

     costheta_jik = (r_ij^2 + r_ik^2 - r_jk^2) / (2*r_ij*r_ik).
  
* and four-body term phi_four can be written as:
  
     phi_four(r_ij, r_ik, r_il, r_jk, r_jl, r_kl) = epsilon * f_4(r1_cap, r2_cap, r3_cap, r4_cap, r5_cap, r6_cap); where
     
     r1_cap = r_ij/sigma, r2_cap = r_ik/sigma, r3_cap = r_il/sigma, r4_cap = r_jk/sigma, r5_cap = r_jl/sigma, r6_cap = r_kl/sigma    and

     f_4(r1_cap, r2_cap, r3_cap, r4_cap, r5_cap, r6_cap) = g(r1_cap, r2_cap, r3_cap, costheta_jik, costheta_jil, costheta_kil); where            

     g(r1_cap, r2_cap, r3_cap, costheta_jik, costheta_jil, costheta_kil) = lambda_2 * (exp(gamma((1/(r1_cap - a)) + (1/(r2_cap - a)) + (1/(r3_cap - a)))))
                                                                         * {1 - exp[-Q((costheta_jik - costheta_0)^2 + (costheta_jil - costheta_0)^2 + (costheta_kil - costheta_0)^2)]}; when r1_cap < a && r2_cap < a && r3_cap < a
                                                                         = 0;                                                                                                       otherwise

     costheta_jik = (r_ij^2 + r_ik^2 - r_jk^2) / (2*r_ij*r_ik);
     costheta_jil = (r_ij^2 + r_il^2 - r_jl^2) / (2*r_ij*r_il);
     costheta_kil = (r_ik^2 + r_il^2 - r_kl^2) / (2*r_ik*r_il).
     
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "KIM_API_C.h"
#include "KIM_API_status.h"

/******************************************************************************
* Below are the definitions and values of all Model parameters
*******************************************************************************/
#define DIM 3       /* dimensionality of space */
#define SPEC1 1  /* internal species code */
#define SPEC2 2  /* internal species code */

/* Define prototypes for Model Driver init */
/* must be all lowercase to be compatible with the KIM API (to support Fortran Tests) */
/**/
int model_driver_init(void* km, char* paramfile_names, int* nmstrlen, int* numparamfiles);

/* Define prototypes for Model (Driver) reinit, compute, and destroy */
/* defined as static to avoid namespace clashes with other Models    */
/**/
static int reinit(void* km);
static int destroy(void* km);
static int compute(void* km);
/**/
static void calc_phi_two(double* A, double* B, double* p, double* q, double* a, double* sigma, double* epsilon, 
                         double r, double* phi);
static void calc_phi_dphi_two(double* A, double* B, double* p, double* q, double* a, double* sigma, double* epsilon,
                              double r, double* phi, double* dphi);
static void calc_phi_d2phi_two(double* A, double* B, double* p, double* q, double* a, double* sigma, double* epsilon,
                               double r, double* phi, double* dphi, double* d2phi);

static void calc_phi_three(double* a, double* lambda, double* gamma, double* sigma, double* epsilon, double* Q, double* costhetat,
                           double rij, double rik, double rjk, double* phi);
static void calc_phi_dphi_three(double* a, double* lambda, double* gamma, double* sigma, double* epsilon, double* Q, double* costhetat, 
                                double rij, double rik, double rjk, double* phi, double* dphi);
static void calc_phi_d2phi_three(double* a, double* lambda, double* gamma, double* sigma, double* epsilon, double* Q, double* costhetat, 
                                 double rij, double rik, double rjk, double* phi, double* dphi, double* d2phi);

static void calc_phi_four(double* a, double* lambda_2, double* gamma, double* sigma, double* epsilon, double* Q, double* costhetat,
                          double rij, double rik, double ril, double rjk, double rjl, double rkl, double* phi);
static void calc_phi_dphi_four(double* a, double* lambda_2, double* gamma, double* sigma, double* epsilon, double* Q, double* costhetat, 
                               double rij, double rik, double ril, double rjk, double rjl, double rkl, double* phi, double* dphi);
static void calc_phi_d2phi_four(double* a, double* lambda_2, double* gamma, double* sigma, double* epsilon, double* Q, double* costhetat, 
                                double rij, double rik, double ril, double rjk, double rjl, double rkl, double* phi, double* dphi, double* d2phi);

/* Define model_buffer structure */
struct model_buffer {
   int NBC;
   int IterOrLoca;
   int energy_ind;
   int forces_ind;
   int particleEnergy_ind;
   int process_dEdr_ind;
   int process_d2Edr2_ind;
   int model_index_shift;
   int numberOfParticles_ind;
   int numberOfSpecies_ind;
   int particleSpecies_ind;
   int coordinates_ind;
   int boxSideLengths_ind;
   int get_neigh_ind;
   int cutoff_ind;

   
   double* Pcutoff;
   double* cutsq;
   double* A;
   double* B;
   double* p;
   double* q;
   double* a;
   double* lambda;
   double* lambda_2;
   double* gamma;
   double* sigma;
   double* epsilon;
   double* Q;
   double* costhetat;
};

/* Calculate two-body term phi_two(r) */
static void calc_phi_two(double* A, double* B, double* p, double* q, double* a, double* sigma, double* epsilon,
                         double r, double* phi)
{
   /* Local variables */
   double r_cap;
   
   r_cap = r/(*sigma);

   if (r_cap >=  *a)
   {
      *phi = 0.0;
   }
   else
   {
      *phi = (*epsilon) * (*A) * ( (*B) * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) ) * exp(1/(r_cap - *a)); 
   }

   return;
}

/* Calculate two-body term phi_two(r) and its derivative dphi_two(r) */
static void calc_phi_dphi_two(double* A, double* B, double* p, double* q, double* a, double* sigma, double* epsilon,
                              double r, double* phi, double* dphi)
{
   /* Local variables */
   double r_cap;
   
   r_cap = r/(*sigma);

   if (r_cap >= *a)
   {
      *phi = 0.0;
      *dphi = 0.0;
   }
   else
   {
      *phi = (*epsilon) * (*A) * ( (*B) * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) ) * exp(1/(r_cap - *a)); 

      *dphi =  ( (*q) * pow(r_cap,-((*q)+1)) - (*p * (*B)) * pow(r_cap,-((*p)+1)) ) 
                  - ( (*B) * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) ) * pow((r_cap - *a),-2);
      *dphi *= (*epsilon / *sigma) * (*A) * exp(1/(r_cap - *a));
   }

   return;
}

/* Calculate two-body term phi_two(r) and its 1st & 2nd derivatives dphi_two(r), d2phi_two(r) */
static void calc_phi_d2phi_two(double* A, double* B, double* p, double* q, double* a, double* sigma, double* epsilon,
                               double r, double* phi, double* dphi, double* d2phi)
{
   /* Local variables */
   double r_cap;
   
   r_cap = r/(*sigma);

   if (r_cap >= *a)
   {
      *phi = 0.0;
      *dphi = 0.0;
      *d2phi = 0.0;
   }
   else
   {
      *phi = (*epsilon) * (*A) * ( (*B) * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) ) * exp(1/(r_cap - *a)); 

      *dphi =  ( (*q) * pow(r_cap,-((*q)+1)) - (*p * *B) * pow(r_cap,-((*p)+1)) ) 
                  - ( (*B) * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) ) * pow((r_cap - *a),-2);
      *dphi *= (*epsilon / *sigma) * (*A) * exp(1/(r_cap - *a));
      
      *d2phi = ( *B * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) ) * ( pow((r_cap - *a),-4) + 2*pow((r_cap - *a),-3) )  
                  + 2 * ( *p * *B *pow(r_cap,-(*p + 1)) - *q * pow(r_cap,-(*q + 1)) ) * pow((r_cap - *a),-2)
                  + ( *p * (*p + 1) *  *B * pow(r_cap,-(*p + 2)) - *q * (*q + 1) * pow(r_cap,-(*q + 2)) );
      *d2phi *= (*epsilon / (*sigma * *sigma)) * (*A) * exp(1/(r_cap - *a));
   }

   return;

}

/* Calculate  three-body term phi_three(r_ij, r_ik, r_jk) */
static void calc_phi_three(double* a, double* lambda, double* gamma, double* sigma, double* epsilon, double* Q, double* costhetat,
                           double rij, double rik, double rjk, double* phi)
{
   /* local variables */
   double c1;
    
   double rij_cap;
   double rik_cap;

   double costhetajik;

   double diff_costhetajik;

   double exp_ij_ik;

   c1 = *lambda * *epsilon;

   rij_cap = rij/(*sigma);
   rik_cap = rik/(*sigma);

   costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);

   diff_costhetajik = costhetajik - *costhetat; 

   exp_ij_ik = exp((*gamma) * (1/(rij_cap - *a) + 1/(rik_cap - *a)));
   
   if ((rij_cap < *a) && (rik_cap < *a)) 
   {
     *phi  = c1 * exp_ij_ik * (1 - exp(- *Q *pow(diff_costhetajik,2))); 
   }
     
   return;
}

/* Calculate three-body term phi_three(r_ij, r_ik, r_jk) and its 1st derivative 
  dphi_three(r_ij, r_ik, r_jk)  
 
  dphi has three components as derivatives of phi w.r.t. r_ij, r_ik, r_jk 
*/
static void calc_phi_dphi_three(double* a, double* lambda, double* gamma, double* sigma, double* epsilon, double* Q, double* costhetat,
                                double rij, double rik, double rjk, double* phi, double* dphi)
{
   /* local variables */
   double c1;
   double c2;
    
   double rij_cap;
   double rik_cap;

   double costhetajik;

   double diff_costhetajik;

   double costhetajik_ij;
   double costhetajik_ik;
   double costhetajik_jk;
   
   double exp_ij_ik;

   double d_ij;
   double d_ik;
   
   double g1, g2, g3, g4;
   
   c1 = *lambda * *epsilon;
   c2 = c1 / *sigma;

   rij_cap = rij/(*sigma);
   rik_cap = rik/(*sigma);

   costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);

   diff_costhetajik = costhetajik - *costhetat;

   /* Derivatives of cosines w.r.t rij, rik, rjk */
   costhetajik_ij = (*sigma)*(pow(rij,2) - pow(rik,2) + pow(rjk,2))/(2*rij*rij*rik);
   costhetajik_ik = (*sigma)*(pow(rik,2) - pow(rij,2) + pow(rjk,2))/(2*rij*rik*rik);
   costhetajik_jk = -(*sigma)*rjk/(rij*rik);

   /* Variables for simplifying terms */
   exp_ij_ik = exp((*gamma) * (1/(rij_cap - *a) + 1/(rik_cap - *a)));

   d_ij = -(*gamma)*pow((rij_cap - *a),-2);
   d_ik = -(*gamma)*pow((rik_cap - *a),-2);

   
   if ((rij_cap < *a) && (rik_cap < *a))
   {
     g1  =  pow(diff_costhetajik,2);
     g2  =  exp(- *Q * g1);
     g3  =  1 - g2;
     g4  =  2 * *Q * g2;        

     *phi  = c1 * exp_ij_ik * g3;
     
     dphi[0] = c2 * exp_ij_ik *  (g3 * d_ij + g4 * (diff_costhetajik * costhetajik_ij)); /* w.r.t. ij */
     dphi[1] = c2 * exp_ij_ik *  (g3 * d_ik + g4 * (diff_costhetajik * costhetajik_ik)); /* w.r.t. ik */
     dphi[2] = c2 * exp_ij_ik *  (g4 * (diff_costhetajik * costhetajik_jk));             /* w.r.t. jk */
   }
   
   return;
}

/* Calculate three-body term phi_three(r_ij, r_ik, r_jk) and its 1st & 2nd derivatives 
  dphi_three(r_ij, r_ik, r_jk), d2phi_three(r_ij, r_ik, r_jk) 
 
  dphi has three components as derivatives of phi w.r.t. r_ij, r_ik, r_jk 

  d2phi as symmetric Hessian matrix of phi has six components: [0]=(ij,ij), [3]=(ij,ik), [4]=(ij,jk)
                                                                            [1]=(ik,ik), [5]=(ik,jk)
                                                                                         [2]=(jk,jk)
*/
static void calc_phi_d2phi_three(double* a, double* lambda, double* gamma, double* sigma, double* epsilon, double* Q, double* costhetat,
                                 double rij, double rik, double rjk, double* phi, double* dphi, double* d2phi)
{
    /* local variables */
   double c1;
   double c2;
   double c3;
    
   double rij_cap;
   double rik_cap;

   double costhetajik;

   double diff_costhetajik;
   
   double costhetajik_ij;
   double costhetajik_ik;
   double costhetajik_jk;
 
   double costhetajik_ij_ij;
   double costhetajik_ik_ik;
   double costhetajik_jk_jk;
   double costhetajik_ij_ik;
   double costhetajik_ij_jk;
   double costhetajik_ik_jk;
   
   double exp_ij_ik;

   double d_ij;
   double d_ik;
   
   double d_ij_2;
   double d_ik_2;
   
   double dd_ij;
   double dd_ik;
   
   double g1, g2, g3, g4, g5, powQ;
   
   double* dg1 = (double *) malloc(3*sizeof(double));
   
   c1 = *lambda * *epsilon;
   c2 = c1 / *sigma;
   c3 = c2 / *sigma;

   rij_cap = rij/(*sigma);
   rik_cap = rik/(*sigma);

   costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);

   diff_costhetajik = costhetajik - *costhetat;
  
   /* Derivatives of cosines w.r.t. r_ij, r_ik, r_jk */
   costhetajik_ij = (*sigma)*(pow(rij,2) - pow(rik,2) + pow(rjk,2))/(2*rij*rij*rik);
   costhetajik_ik = (*sigma)*(pow(rik,2) - pow(rij,2) + pow(rjk,2))/(2*rij*rik*rik);
   costhetajik_jk = -(*sigma)*rjk/(rij*rik);

   /* Hessian matrix of cosine */
   costhetajik_ij_ij = (*sigma * *sigma)*(pow(rik,2) - pow(rjk,2))/(rij*rij*rij*rik);
   costhetajik_ik_ik = (*sigma * *sigma)*(pow(rij,2) - pow(rjk,2))/(rij*rik*rik*rik);
   costhetajik_jk_jk = -(*sigma * *sigma)/(rij*rik);
   costhetajik_ij_ik = -(*sigma * *sigma)*(pow(rij,2) + pow(rik,2) + pow(rjk,2))/(2*rij*rij*rik*rik);
   costhetajik_ij_jk = (*sigma * *sigma)*rjk/(rij*rij*rik); 
   costhetajik_ik_jk = (*sigma * *sigma)*rjk/(rik*rik*rij); 
   
   /* Variables for simplifying terms */
   exp_ij_ik = exp((*gamma) * (1/(rij_cap - *a) + 1/(rik_cap - *a)));

   d_ij = -(*gamma)*pow((rij_cap - *a),-2);
   d_ik = -(*gamma)*pow((rik_cap - *a),-2);

   d_ij_2 = d_ij * d_ij;
   d_ik_2 = d_ik * d_ik;

   dd_ij = (2 * *gamma)*pow((rij_cap - *a),-3);
   dd_ik = (2 * *gamma)*pow((rik_cap - *a),-3);

   powQ = 4 * *Q * *Q;

   dg1[0]  =  diff_costhetajik * costhetajik_ij;
   dg1[1]  =  diff_costhetajik * costhetajik_ik;
   dg1[2]  =  diff_costhetajik * costhetajik_jk;
   
   if ((rij_cap < *a) && (rik_cap < *a))
   {
     g1  =  pow(diff_costhetajik,2);
     g2  =  exp(- *Q * g1);
     g3  =  1 - g2;
     g4  =  2 * *Q * g2;
     g5  =  powQ * g2;
     
     *phi  = c1 * exp_ij_ik * g3;
     
     /* 1st derivative */
     dphi[0] = c2 * exp_ij_ik *  (g3 * d_ij + g4 * (diff_costhetajik * costhetajik_ij)); /* w.r.t. ij */
     dphi[1] = c2 * exp_ij_ik *  (g3 * d_ik + g4 * (diff_costhetajik * costhetajik_ik)); /* w.r.t. ik */
     dphi[2] = c2 * exp_ij_ik *  (g4 * (diff_costhetajik * costhetajik_jk));             /* w.r.t. jk */
     
     /* Hessian */
     d2phi[0] = c3 * exp_ij_ik * ((g3 * (d_ij_2 + dd_ij)) + (2 * g4 * d_ij * dg1[0]) - (g5 * pow(dg1[0],2)) + (g4 * (diff_costhetajik * costhetajik_ij_ij +  pow(costhetajik_ij,2))));          
     d2phi[1] = c3 * exp_ij_ik * ((g3 * (d_ik_2 + dd_ik)) + (2 * g4 * d_ik * dg1[1]) - (g5 * pow(dg1[1],2)) + (g4 * (diff_costhetajik * costhetajik_ik_ik + pow(costhetajik_ik,2))));          
     d2phi[2] = c3 * exp_ij_ik * (- (g5 * pow(dg1[2],2)) + (g4 * (diff_costhetajik * costhetajik_jk_jk +  pow(costhetajik_jk,2))));                                                            
     d2phi[3] = c3 * exp_ij_ik * ((g3 * d_ij * d_ik) + g4 * (d_ik * dg1[0] + d_ij * dg1[1]) - (g5 * dg1[0] * dg1[1]) + (g4 * (diff_costhetajik * costhetajik_ij_ik + costhetajik_ij * costhetajik_ik))); 
     d2phi[4] = c3 * exp_ij_ik * ((g4 * d_ik * dg1[2]) - (g5 * dg1[1] * dg1[2]) + (g4 * (diff_costhetajik * costhetajik_ik_jk + costhetajik_ik * costhetajik_jk)));                                      
     d2phi[5] = c3 * exp_ij_ik * ((g4 * d_ij * dg1[2]) - (g5 * dg1[0] * dg1[2]) + (g4 * (diff_costhetajik * costhetajik_ij_jk + costhetajik_ij * costhetajik_jk)));                                      
   }

   /* d2phi[0] derivative is w.r.t. rij, rij */
   /* d2phi[1] derivative is w.r.t. rik, rik */
   /* d2phi[2] derivative is w.r.t. rjk, rjk */
   /* d2phi[3] derivative is w.r.t. rij, rik */
   /* d2phi[4] derivative is w.r.t. rij, rjk */
   /* d2phi[5] derivative is w.r.t. rik, rjk */
   
   free(dg1);

   return;
}

/* Calculate  four-body term phi_four(r_ij, r_ik, r_il, r_jk, r_jl, r_kl) */
static void calc_phi_four(double* a, double* lambda_2, double* gamma, double* sigma, double* epsilon, double* Q, double* costhetat,
                          double rij, double rik, double ril, double rjk, double rjl, double rkl, double* phi)
{
   /* local variables */
   double c1;
    
   double rij_cap;
   double rik_cap;
   double ril_cap;

   double costhetajik, costhetajil, costhetakil;

   double diff_costhetajik, diff_costhetajil, diff_costhetakil;

   double exp_ij_ik_il;

   c1 = *lambda_2 * *epsilon;

   rij_cap = rij/(*sigma);
   rik_cap = rik/(*sigma);
   ril_cap = ril/(*sigma);

   costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);
   costhetajil = (pow(rij,2) + pow(ril,2) - pow(rjl,2))/(2*rij*ril);
   costhetakil = (pow(ril,2) + pow(rik,2) - pow(rkl,2))/(2*ril*rik);

   /* Difference of two cosines */
   diff_costhetajik = costhetajik - *costhetat; diff_costhetajil = costhetajil - *costhetat; diff_costhetakil = costhetakil - *costhetat;   

   /* Variables for simplifying terms */
   exp_ij_ik_il = exp((*gamma) * (1/(rij_cap - *a) + 1/(rik_cap - *a) + 1/(ril_cap - *a)));

   if ((rij_cap < *a) && (rik_cap < *a) && (ril_cap < *a)) 
   {
     *phi  = c1 * exp_ij_ik_il * (1 - exp(- *Q *(pow(diff_costhetajik,2) + pow(diff_costhetajil,2) + pow(diff_costhetakil,2)))); 
   }
     
   return;
}

/* Calculate four-body term phi_four(r_ij, r_ik, r_il, r_jk, r_jl, r_kl) and its 1st derivative 
  dphi_four(r_ij, r_ik, r_il, r_jk, r_jl, r_kl)  
 
  dphi has six components as derivatives of phi w.r.t. r_ij, r_ik, r_il, r_jk, r_jl, r_kl 
*/
static void calc_phi_dphi_four(double* a, double* lambda_2, double* gamma, double* sigma, double* epsilon, double* Q, double* costhetat, 
                               double rij, double rik, double ril, double rjk, double rjl, double rkl, double* phi, double* dphi)
{
   /* local variables */
   double c1;
   double c2;
    
   double rij_cap;
   double rik_cap;
   double ril_cap;

   double costhetajik, costhetajil, costhetakil;

   double diff_costhetajik, diff_costhetajil, diff_costhetakil;

   /* 1st derivate of cosines w.r.t. ij, ik, il, jk, jl, kl in array order*/
   double* d_costhetajik = (double *) malloc(3*sizeof(double));
   double* d_costhetajil = (double *) malloc(3*sizeof(double));
   double* d_costhetakil = (double *) malloc(3*sizeof(double));

   /* Variables for simplifying terms */
   double exp_ij_ik_il;
   
   double d_ij;
   double d_ik;
   double d_il;

   double g1, g2, g3, g4;
   
   c1 = *lambda_2 * *epsilon;
   c2 = c1 / *sigma;

   rij_cap = rij/(*sigma);
   rik_cap = rik/(*sigma);
   ril_cap = ril/(*sigma);

   costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);
   costhetajil = (pow(rij,2) + pow(ril,2) - pow(rjl,2))/(2*rij*ril);
   costhetakil = (pow(ril,2) + pow(rik,2) - pow(rkl,2))/(2*ril*rik);

   /* Difference of two cosines */
   diff_costhetajik = costhetajik - *costhetat; diff_costhetajil = costhetajil - *costhetat; diff_costhetakil = costhetakil - *costhetat;   

   /* 1st derivative of cosines */
   d_costhetajik[0] = (*sigma)*(pow(rij,2) - pow(rik,2) + pow(rjk,2))/(2*rij*rij*rik); /* w.r.t. ij */
   d_costhetajik[1] = (*sigma)*(pow(rik,2) - pow(rij,2) + pow(rjk,2))/(2*rij*rik*rik); /* w.r.t. ik */
   d_costhetajik[2] = -(*sigma)*rjk/(rij*rik); /* w.r.t. jk */

   d_costhetajil[0] = (*sigma)*(pow(rij,2) - pow(ril,2) + pow(rjl,2))/(2*rij*rij*ril); /* w.r.t. ij */
   d_costhetajil[1] = (*sigma)*(pow(ril,2) - pow(rij,2) + pow(rjl,2))/(2*rij*ril*ril); /* w.r.t. il */
   d_costhetajil[2] = -(*sigma)*rjl/(rij*ril); /* w.r.t. jl */
   
   d_costhetakil[0] = (*sigma)*(pow(rik,2) - pow(ril,2) + pow(rkl,2))/(2*rik*rik*ril); /* w.r.t. ik */ 
   d_costhetakil[1] = (*sigma)*(pow(ril,2) - pow(rik,2) + pow(rkl,2))/(2*rik*ril*ril); /* w.r.t. il */
   d_costhetakil[2] = -(*sigma)*rkl/(rik*ril); /* w.r.t. kl */

   /* Variables for simplifying terms */
   exp_ij_ik_il = exp((*gamma) * (1/(rij_cap - *a) + 1/(rik_cap - *a) + 1/(ril_cap - *a)));

   d_ij = -(*gamma)*pow((rij_cap - *a),-2);
   d_ik = -(*gamma)*pow((rik_cap - *a),-2);
   d_il = -(*gamma)*pow((ril_cap - *a),-2);

   
   if ((rij_cap < *a) && (rik_cap < *a) && (ril_cap < *a)) 
   {
     g1  =  pow(diff_costhetajik,2) + pow(diff_costhetajil,2) + pow(diff_costhetakil,2);
     g2  =  exp(- *Q * g1);
     g3  =  1 - g2;
     g4  =  2 * *Q * g2;        

     *phi  = c1 * exp_ij_ik_il * g3;
     
     dphi[0] = c2 * exp_ij_ik_il *  (g3 * d_ij + g4 * (diff_costhetajik * d_costhetajik[0] 
                                                 + diff_costhetajil * d_costhetajil[0])); /* w.r.t. ij */
     dphi[1] = c2 * exp_ij_ik_il *  (g3 * d_ik + g4 * (diff_costhetajik * d_costhetajik[1]                 
                                                 + diff_costhetakil * d_costhetakil[0])); /* w.r.t. ik */
     dphi[2] = c2 * exp_ij_ik_il *  (g3 * d_il + g4 * (diff_costhetajil * d_costhetajil[1]                 
                                                 + diff_costhetakil * d_costhetakil[1])); /* w.r.t. il */
     dphi[3] = c2 * exp_ij_ik_il *  (g4 * (diff_costhetajik * d_costhetajik[2]));             /* w.r.t. jk */
     dphi[4] = c2 * exp_ij_ik_il *  (g4 * (diff_costhetajil * d_costhetajil[2]));             /* w.r.t. jl */
     dphi[5] = c2 * exp_ij_ik_il *  (g4 * (diff_costhetakil * d_costhetakil[2]));             /* w.r.t. kl */
   }
     
    free(d_costhetajik);
    free(d_costhetajil);
    free(d_costhetakil);
    
   return;
}

/* Calculate four-body term phi_four(r_ij, r_ik, r_il, r_jk, r_jl, r_kl) and its 1st & 2nd derivatives 
  dphi_four(r_ij, r_ik, r_il, r_jk, r_jl, r_kl), d2phi_four(r_ij, r_ik, r_il, r_jk, r_jl, r_kl) 
 
  dphi has six components as derivatives of phi w.r.t. r_ij, r_ik, r_il, r_jk, r_jl, r_kl 

  d2phi as symmetric Hessian matrix of phi has 21 components: [0]=(ij,ij), [6]=(ij,ik), [11]=(ij,il), [15]=(ij,jk), [18]=(ij,jl), [20]=(ij,kl)
                                                                           [1]=(ik,ik), [7] =(ik,il), [12]=(ik,jk), [16]=(ik,jl), [19]=(ik,kl) 
                                                                                        [2] =(il,il), [8] =(il,jk), [13]=(il,jl), [17]=(il,kl) 
                                                                                                      [3] =(jk,jk), [9] =(jk,jl), [14]=(jk,kl) 
                                                                                                                    [4] =(jl,jl), [10]=(jl,kl) 
                                                                                                                                  [5] =(kl,kl) 
*/
static void calc_phi_d2phi_four(double* a, double* lambda_2, double* gamma, double* sigma, double* epsilon, double* Q, double* costhetat, 
                                double rij, double rik, double ril, double rjk, double rjl, double rkl, double* phi, double* dphi, double* d2phi)
{
    /* local variables */
   double c1;
   double c2;
   double c3;
    
   double rij_cap;
   double rik_cap;
   double ril_cap;

   double costhetajik, costhetajil, costhetakil;

   double diff_costhetajik, diff_costhetajil, diff_costhetakil;

   /* 1st derivate of cosines w.r.t. ij, ik, il, jk, jl, kl in array order*/
   double* d_costhetajik = (double *) malloc(3*sizeof(double));
   double* d_costhetajil = (double *) malloc(3*sizeof(double));
   double* d_costhetakil = (double *) malloc(3*sizeof(double));

   /* Hessian of cosines in the given array order*/
   double* d2_costhetajik = (double *) malloc(6*sizeof(double));  /* [0] = (ij,ij), [1] = (ik,ik), [2] = (jk,jk), [3] = (ij,ik), [4] = (ij,jk), [5] = (ik,jk) */
   double* d2_costhetajil = (double *) malloc(6*sizeof(double));
   double* d2_costhetakil = (double *) malloc(6*sizeof(double));

   double exp_ij_ik_il;
   
   double d_ij;
   double d_ik;
   double d_il;
   
   double d_ij_2;
   double d_ik_2;
   double d_il_2;
   
   double dd_ij;
   double dd_ik;
   double dd_il;
   
   double g1, g2, g3, g4, g5, powQ;
   
   double* dg1 = (double *) malloc(6*sizeof(double));
   double* dg2 = (double *) malloc(6*sizeof(double));
   
   int i;
   
   c1 = *lambda_2 * *epsilon;
   c2 = c1 / *sigma;
   c3 = c2 / *sigma;

   rij_cap = rij/(*sigma);
   rik_cap = rik/(*sigma);
   ril_cap = ril/(*sigma);

   costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);
   costhetajil = (pow(rij,2) + pow(ril,2) - pow(rjl,2))/(2*rij*ril);
   costhetakil = (pow(ril,2) + pow(rik,2) - pow(rkl,2))/(2*ril*rik);

   /* Difference of two cosines */
   diff_costhetajik = costhetajik - *costhetat; diff_costhetajil = costhetajil - *costhetat; diff_costhetakil = costhetakil - *costhetat;   

   /* 1st derivative of cosines */
   d_costhetajik[0] = (*sigma)*(pow(rij,2) - pow(rik,2) + pow(rjk,2))/(2*rij*rij*rik); /* w.r.t. ij */
   d_costhetajik[1] = (*sigma)*(pow(rik,2) - pow(rij,2) + pow(rjk,2))/(2*rij*rik*rik); /* w.r.t. ik */
   d_costhetajik[2] = -(*sigma)*rjk/(rij*rik); /* w.r.t. jk */

   d_costhetajil[0] = (*sigma)*(pow(rij,2) - pow(ril,2) + pow(rjl,2))/(2*rij*rij*ril); /* w.r.t. ij */
   d_costhetajil[1] = (*sigma)*(pow(ril,2) - pow(rij,2) + pow(rjl,2))/(2*rij*ril*ril); /* w.r.t. il */
   d_costhetajil[2] = -(*sigma)*rjl/(rij*ril); /* w.r.t. jl */
   
   d_costhetakil[0] = (*sigma)*(pow(rik,2) - pow(ril,2) + pow(rkl,2))/(2*rik*rik*ril); /* w.r.t. ik */ 
   d_costhetakil[1] = (*sigma)*(pow(ril,2) - pow(rik,2) + pow(rkl,2))/(2*rik*ril*ril); /* w.r.t. il */
   d_costhetakil[2] = -(*sigma)*rkl/(rik*ril); /* w.r.t. kl */

   /* Hessian matrix of cosines */
   d2_costhetajik[0] = (*sigma * *sigma)*(pow(rik,2) - pow(rjk,2))/(pow(rij,3)*rik);                         /* w.r.t. ij,ij */
   d2_costhetajik[1] = -(*sigma * *sigma)*(pow(rij,2) + pow(rik,2) + pow(rjk,2))/(2*pow(rij,2)*pow(rik,2));  /* w.r.t. ij,ik */
   d2_costhetajik[2] = (*sigma * *sigma)*rjk/(pow(rij,2)*rik);                                               /* w.r.t. ij,jk */
   d2_costhetajik[3] = (*sigma * *sigma)*(pow(rij,2) - pow(rjk,2))/(rij*pow(rik,3));                         /* w.r.t. ik,ik */
   d2_costhetajik[4] = (*sigma * *sigma)*rjk/(pow(rik,2)*rij);                                               /* w.r.t. ik,jk */
   d2_costhetajik[5] = -(*sigma * *sigma)/(rij*rik);                                                         /* w.r.t. jk,jk */

   d2_costhetajil[0] = (*sigma * *sigma)*(pow(ril,2) - pow(rjl,2))/(pow(rij,3)*ril);                         /* w.r.t. ij,ij */
   d2_costhetajil[1] = -(*sigma * *sigma)*(pow(rij,2) + pow(ril,2) + pow(rjl,2))/(2*pow(rij,2)*pow(ril,2));  /* w.r.t. ij,il */
   d2_costhetajil[2] = (*sigma * *sigma)*rjl/(pow(rij,2)*ril);                                               /* w.r.t. ij,jl */
   d2_costhetajil[3] = (*sigma * *sigma)*(pow(rij,2) - pow(rjl,2))/(rij*pow(ril,3));                         /* w.r.t. il,il */
   d2_costhetajil[4] = (*sigma * *sigma)*rjl/(pow(ril,2)*rij);                                               /* w.r.t. il,jl */
   d2_costhetajil[5] = -(*sigma * *sigma)/(rij*ril);                                                         /* w.r.t. jl,jl */

   d2_costhetakil[0] = (*sigma * *sigma)*(pow(ril,2) - pow(rkl,2))/(pow(rik,3)*ril);                         /* w.r.t. ik,ik */
   d2_costhetakil[1] = -(*sigma * *sigma)*(pow(rik,2) + pow(ril,2) + pow(rkl,2))/(2*pow(rik,2)*pow(ril,2));  /* w.r.t. ik,il */
   d2_costhetakil[2] = (*sigma * *sigma)*rkl/(pow(rik,2)*ril);                                               /* w.r.t. ik,kl */
   d2_costhetakil[3] = (*sigma * *sigma)*(pow(rik,2) - pow(rkl,2))/(rik*pow(ril,3));                         /* w.r.t. il,il */
   d2_costhetakil[4] = (*sigma * *sigma)*rkl/(pow(ril,2)*rik);                                               /* w.r.t. il,kl */
   d2_costhetakil[5] = -(*sigma * *sigma)/(rik*ril);                                                         /* w.r.t. kl,kl */

   /* Variables for simplifying terms */
   exp_ij_ik_il = exp((*gamma) * (1/(rij_cap - *a) + 1/(rik_cap - *a) + 1/(ril_cap - *a)));

   d_ij = -(*gamma)*pow((rij_cap - *a),-2);
   d_ik = -(*gamma)*pow((rik_cap - *a),-2);
   d_il = -(*gamma)*pow((ril_cap - *a),-2);

   d_ij_2 = d_ij * d_ij;
   d_ik_2 = d_ik * d_ik;
   d_il_2 = d_il * d_il;

   dd_ij = (2 * *gamma)*pow((rij_cap - *a),-3);
   dd_ik = (2 * *gamma)*pow((rik_cap - *a),-3);
   dd_il = (2 * *gamma)*pow((ril_cap - *a),-3);

   *phi = 0.0;

   for(i = 0; i < 6; i++) dphi[i] = 0.0;
   for(i = 0; i < 21; i++) d2phi[i] = 0.0;

   powQ = 4 * *Q * *Q;

   if ((rij_cap < *a) && (rik_cap < *a) && (ril_cap < *a)) 
   {
     g1  =  pow(diff_costhetajik,2) + pow(diff_costhetajil,2) + pow(diff_costhetakil,2);
     g2  =  exp(- *Q * g1);
     g3  =  1 - g2;
     g4  =  2 * *Q * g2;
     g5  =  powQ * g2;

     dg1[0]  =  diff_costhetajik * d_costhetajik[0] + diff_costhetajil * d_costhetajil[0];
     dg1[1]  =  diff_costhetajik * d_costhetajik[1] + diff_costhetakil * d_costhetakil[0];
     dg1[2]  =  diff_costhetajil * d_costhetajil[1] + diff_costhetakil * d_costhetakil[1];
     dg1[3]  =  diff_costhetajik * d_costhetajik[2];
     dg1[4]  =  diff_costhetajil * d_costhetajil[2];
     dg1[5]  =  diff_costhetakil * d_costhetakil[2];

     *phi  += exp_ij_ik_il * g3;
     
     /* 1st derivative */
     dphi[0] += exp_ij_ik_il *  (g3 * d_ij + g4 * dg1[0]); /* w.r.t. ij */
     dphi[1] += exp_ij_ik_il *  (g3 * d_ik + g4 * dg1[1]); /* w.r.t. ik */                
     dphi[2] += exp_ij_ik_il *  (g3 * d_il + g4 * dg1[2]); /* w.r.t. il */                 
     dphi[3] += exp_ij_ik_il *  g4 * dg1[3];               /* w.r.t. jk */
     dphi[4] += exp_ij_ik_il *  g4 * dg1[4];               /* w.r.t. jl */
     dphi[5] += exp_ij_ik_il *  g4 * dg1[5];               /* w.r.t. kl */
    
     /* Hessian */
     dg2[0]  =  diff_costhetajik * d2_costhetajik[0] + diff_costhetajil * d2_costhetajil[0] + pow(d_costhetajik[0],2) + pow(d_costhetajil[0],2); 
     dg2[1]  =  diff_costhetajik * d2_costhetajik[1] + d_costhetajik[0] * d_costhetajik[1];
     dg2[2]  =  diff_costhetajil * d2_costhetajil[1] + d_costhetajil[0] * d_costhetajil[1];
     dg2[3]  =  diff_costhetajik * d2_costhetajik[2] + d_costhetajik[0] * d_costhetajik[2];
     dg2[4]  =  diff_costhetajil * d2_costhetajil[2] + d_costhetajil[0] * d_costhetajil[2];

     d2phi[0] += exp_ij_ik_il * ((g3 * (d_ij_2 + dd_ij)) + (2 * g4 * d_ij * dg1[0]) - (g5 * pow(dg1[0],2)) + (g4 * dg2[0]));          /* w.r.t. (ij,ij) */
     d2phi[6] += exp_ij_ik_il * ((g3 * d_ij * d_ik) + g4 * (d_ik * dg1[0] + d_ij * dg1[1]) - (g5 * dg1[0] * dg1[1]) + (g4 * dg2[1])); /* w.r.t. (ij,ik) */
     d2phi[11] += exp_ij_ik_il * ((g3 * d_ij * d_il) + g4 * (d_il * dg1[0] + d_ij * dg1[2]) - (g5 * dg1[0] * dg1[2]) + (g4 * dg2[2])); /* w.r.t. (ij,il) */
     d2phi[15] += exp_ij_ik_il * ((g4 * d_ij * dg1[3]) - (g5 * dg1[0] * dg1[3]) + (g4 * dg2[3]));                                      /* w.r.t. (ij,jk) */
     d2phi[18] += exp_ij_ik_il * ((g4 * d_ij * dg1[4]) - (g5 * dg1[0] * dg1[4]) + (g4 * dg2[4]));                                      /* w.r.t. (ij,jl) */
     d2phi[20] += exp_ij_ik_il * ((g4 * d_ij * dg1[5]) - (g5 * dg1[0] * dg1[5]));                                                      /* w.r.t. (ij,kl) */
     
     dg2[0]  =  diff_costhetajik * d2_costhetajik[3] + diff_costhetakil * d2_costhetakil[0] + pow(d_costhetajik[1],2) + pow(d_costhetakil[0],2);
     dg2[1]  =  diff_costhetakil * d2_costhetakil[1] + d_costhetakil[0] * d_costhetakil[1];
     dg2[2]  =  diff_costhetajik * d2_costhetajik[4] + d_costhetajik[1] * d_costhetajik[2];
     dg2[3]  =  diff_costhetakil * d2_costhetakil[2] + d_costhetakil[0] * d_costhetakil[2];

     d2phi[1] += exp_ij_ik_il * ((g3 * (d_ik_2 + dd_ik)) + (2 * g4 * d_ik * dg1[1]) - (g5 * pow(dg1[1],2)) + (g4 * dg2[0]));          /* w.r.t. (ik,ik) */
     d2phi[7] += exp_ij_ik_il * ((g3 * d_ik * d_il) + g4 * (d_il * dg1[1] + d_ik * dg1[2]) - (g5 * dg1[1] * dg1[2]) + (g4 * dg2[1])); /* w.r.t. (ik,il) */
     d2phi[12] += exp_ij_ik_il * ((g4 * d_ik * dg1[3]) - (g5 * dg1[1] * dg1[3]) + (g4 * dg2[2]));                                      /* w.r.t. (ik,jk) */
     d2phi[16] += exp_ij_ik_il * ((g4 * d_ik * dg1[4]) - (g5 * dg1[1] * dg1[4]));                                                      /* w.r.t. (ik,jl) */
     d2phi[19] += exp_ij_ik_il * ((g4 * d_ik * dg1[5]) - (g5 * dg1[1] * dg1[5]) + (g4 * dg2[3]));                                     /* w.r.t. (ik,kl) */

     dg2[0]  =  diff_costhetajil * d2_costhetajil[3] + diff_costhetakil * d2_costhetakil[3] + pow(d_costhetajil[1],2) + pow(d_costhetakil[1],2);
     dg2[1]  =  diff_costhetajil * d2_costhetajil[4] + d_costhetajil[1] * d_costhetajil[2];
     dg2[2]  =  diff_costhetakil * d2_costhetakil[4] + d_costhetakil[1] * d_costhetakil[2];

     d2phi[2] += exp_ij_ik_il * ((g3 * (d_il_2 + dd_il)) + (2 * g4 * d_il * dg1[2]) - (g5 * pow(dg1[2],2)) + (g4 * dg2[0]));         /* w.r.t. (il,il) */
     d2phi[8] += exp_ij_ik_il * ((g4 * d_il * dg1[3]) - (g5 * dg1[2] * dg1[3]));                                                     /* w.r.t. (il,jk) */
     d2phi[13] += exp_ij_ik_il * ((g4 * d_il * dg1[4]) - (g5 * dg1[2] * dg1[4]) + (g4 * dg2[1]));                                     /* w.r.t. (il,jl) */
     d2phi[17] += exp_ij_ik_il * ((g4 * d_il * dg1[5]) - (g5 * dg1[2] * dg1[5]) + (g4 * dg2[2]));                                     /* w.r.t. (il,kl) */
     
     dg2[0]  =  diff_costhetajik * d2_costhetajik[5] +  pow(d_costhetajik[2],2);
     
     d2phi[3] += exp_ij_ik_il * (- (g5 * pow(dg1[3],2)) + (g4 * dg2[0]));                                                            /* w.r.t. (jk,jk) */
     d2phi[9] += exp_ij_ik_il * (- g5 * dg1[3] * dg1[4]);                                                                            /* w.r.t. (jk,jl) */
     d2phi[14] += exp_ij_ik_il * (- g5 * dg1[3] * dg1[5]);                                                                            /* w.r.t. (jk,kl) */
     
     dg2[0]  =  diff_costhetajil * d2_costhetajil[5] +  pow(d_costhetajil[2],2);
     
     d2phi[4] += exp_ij_ik_il * (- (g5 * pow(dg1[4],2)) + (g4 * dg2[0]));                                                            /* w.r.t. (jl,jl) */
     d2phi[10] += exp_ij_ik_il * (- g5 * dg1[4] * dg1[5]);                                                                            /* w.r.t. (jl,kl) */
     
     dg2[0]  =  diff_costhetakil * d2_costhetakil[5] +  pow(d_costhetakil[2],2);
     
     d2phi[5] += exp_ij_ik_il * (- (g5 * pow(dg1[5],2)) + (g4 * dg2[0]));                                                            /* w.r.t. (kl,kl) */
   }

   *phi  *= c1;

   for(i = 0; i < 6; i++) dphi[i] *= c2;
   for(i = 0; i < 21; i++) d2phi[i] *= c3;
   
    free(d_costhetajik);
    free(d_costhetajil);
    free(d_costhetakil);
    
    free(d2_costhetajik);
    free(d2_costhetajil);
    free(d2_costhetakil);
    
    free(dg1);
    free(dg2);

   return;
}

/* compute function */
static int compute(void* km) 
{
   /* local variables */
   intptr_t* pkim = *((intptr_t**) km);
   double R1, R2, R3, R4, R5, R6;
   double R_pairs[2];
   double *pR_pairs = &(R_pairs[0]);
   double Rsqij, Rsqik, Rsqil, Rsqjk, Rsqjl, Rsqkl;
   double phi_two;
   double dphi_two;
   double d2phi_two;
   double dEidr_two;
   double d2Eidr_two;
   double phi_three;
   double* dphi_three;
   double* d2phi_three;
   double* dEidr_three;
   double* d2Eidr_three;
   double phi_four;
   double* dphi_four;
   double* d2phi_four;
   double* dEidr_four;
   double* d2Eidr_four;
   double Rij[DIM];
   double Rik[DIM];
   double Ril[DIM];
   double Rjk[DIM];
   double Rjl[DIM];
   double Rkl[DIM];
   double *pRij = &(Rij[0]);
   double *pRik = &(Rik[0]);
   double *pRil = &(Ril[0]);
   double *pRjk = &(Rjk[0]);
   double *pRjl = &(Rjl[0]);
   double *pRkl = &(Rkl[0]);
   double Rij_pairs[2][3];
   double *pRij_pairs = &(Rij_pairs[0][0]);
   int ier;
   int i;
   int i_pairs[2];
   int *pi_pairs = &(i_pairs[0]);
   int j;
   int j_pairs[2];
   int *pj_pairs = &(j_pairs[0]);
   int jj;
   int k;
   int kk;
   int l;
   int ll;
   int kdim;
   int currentAtom;
   int* neighListOfCurrentAtom;
   struct model_buffer* buffer;
   int comp_energy;
   int comp_force;
   int comp_particleEnergy;
   int comp_process_dEdr;
   int comp_process_d2Edr2;
   int NBC;
   int IterOrLoca;
   int model_index_shift;
   int zero = 0;
   int one = 1;
   int request;
   
   int* nAtoms;
   int* nSpecies;
   int* particleSpecies;
   double* cutoff;
   double* cutsq;
   double* A;
   double* B;
   double* p;
   double* q;
   double* a;
   double* lambda;
   double* lambda_2;
   double* gamma;
   double* sigma;
   double* epsilon;
   double* Q;
   double* costhetat;
   double* Rij_list;
   double* coords;
   double* energy;
   double* force;
   double* particleEnergy;
   double* boxSideLengths;
   int numOfAtomNeigh;
   int iSpecies;
   int jSpecies;
   int kSpecies;
   int lSpecies;
   int interaction_index;
   typedef int (*get_neigh_ptr)(void *,int *,int *,int *, int *, int **, double **);
   get_neigh_ptr get_neigh;

   dphi_three = (double *) malloc(3*sizeof(double));
   d2phi_three = (double *) malloc(6*sizeof(double));
   dEidr_three = (double *) malloc(3*sizeof(double));
   d2Eidr_three = (double *) malloc(6*sizeof(double));

   dphi_four = (double *) malloc(6*sizeof(double));
   d2phi_four = (double *) malloc(21*sizeof(double));
   dEidr_four = (double *) malloc(6*sizeof(double));
   d2Eidr_four = (double *) malloc(21*sizeof(double));

   /* get buffer from KIM object */
   buffer = (struct model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
      return ier;
   }

   /* unpack info from the buffer */
   NBC = buffer->NBC;
   IterOrLoca = buffer->IterOrLoca;
   model_index_shift = buffer->model_index_shift;
   
   /* check to see if we have been asked to compute the forces, particleEnergy, and dEdr */
   KIM_API_getm_compute_by_index(pkim, &ier, 5*3,
                                 buffer->energy_ind,         &comp_energy,         1,
                                 buffer->forces_ind,         &comp_force,          1,
                                 buffer->particleEnergy_ind, &comp_particleEnergy, 1,
                                 buffer->process_dEdr_ind,   &comp_process_dEdr,   1,
                                 buffer->process_d2Edr2_ind, &comp_process_d2Edr2, 1);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute_by_index", ier);
      return ier;
   }

   KIM_API_getm_data_by_index(pkim, &ier, 9*3,
                              buffer->cutoff_ind,                      &cutoff,         1,
                              buffer->numberOfParticles_ind,           &nAtoms,         1,
                              buffer->numberOfSpecies_ind,             &nSpecies,       1,
                              buffer->particleSpecies_ind,             &particleSpecies,1,
                              buffer->coordinates_ind,                 &coords,         1,
                              buffer->boxSideLengths_ind,              &boxSideLengths, (NBC==2),
                              buffer->energy_ind,                      &energy,         comp_energy,
                              buffer->forces_ind,                      &force,          comp_force,
                              buffer->particleEnergy_ind,              &particleEnergy, comp_particleEnergy);
     
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data_by_index", ier);
      return ier;
   }
   /* unpack the Model's parameters stored in the KIM API object */
   if(*nSpecies == 1)
   {
        interaction_index = 0; 
        A        = &(buffer->A)[interaction_index];
        B        = &(buffer->B)[interaction_index];
        p        = &(buffer->p)[interaction_index];
        q        = &(buffer->q)[interaction_index];
        lambda   = &(buffer->lambda)[interaction_index];
        lambda_2 = &(buffer->lambda_2)[interaction_index];
        gamma    = &(buffer->gamma)[interaction_index];
        sigma    = &(buffer->sigma)[interaction_index];
        epsilon  = &(buffer->epsilon)[interaction_index];
        Q        = &(buffer->Q)[interaction_index];
        costhetat= &(buffer->costhetat)[interaction_index];
        a   = &(buffer->a)[interaction_index];
        cutsq    = &(buffer->cutsq)[interaction_index];
   }
   
   if (NBC!=3)
   {
      get_neigh = (get_neigh_ptr) KIM_API_get_method_by_index(pkim, buffer->get_neigh_ind, &ier);
      if (KIM_STATUS_OK > ier)
      {
         KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_method_by_index", ier);
         return ier;
      }
   }

   /* Check to be sure that the atom types are correct */
   /**/
   ier = KIM_STATUS_FAIL; /* assume an error */
   for (i = 0; i < *nAtoms; ++i)
   {
      if ( particleSpecies[i] > 2)
      {
         KIM_API_report_error(__LINE__, __FILE__, "Unexpected species type detected", ier);
         return ier;
      }
   }
   ier = KIM_STATUS_OK; /* everything is ok */

   /* initialize potential energies, forces, and virial term */
   if (comp_particleEnergy)
   {
      for (i = 0; i < *nAtoms; ++i)
      {
         particleEnergy[i] = 0.0;
      }
   }
   else if (comp_energy)
   {
      *energy = 0.0;
   }

   if (comp_force)
   {
      for (i = 0; i < *nAtoms; ++i)
      {
         for (kdim = 0; kdim < DIM; ++kdim)
         {
            force[i*DIM + kdim] = 0.0;
         }
      }
   }

   /* Initialize neighbor handling for CLUSTER NBC */
   if (3 == NBC) /* CLUSTER */
   {
      neighListOfCurrentAtom = (int *) malloc((*nAtoms)*sizeof(int));
   }

   /* Initialize neighbor handling for Iterator mode */

   if (1 == IterOrLoca)
   {
      ier = (*get_neigh)(&pkim, &zero, &zero, &currentAtom, &numOfAtomNeigh,
                          &neighListOfCurrentAtom, &Rij_list);
      /* check for successful initialization */
      if (KIM_STATUS_NEIGH_ITER_INIT_OK != ier)
      {
         KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", ier);
         ier = KIM_STATUS_FAIL;
         return ier;
      }
   }

   /* Compute enery and forces */

   /* loop over particles and compute enregy and forces */
   i = -1;
   while( 1 )
   {

      /* Set up neighbor list for next atom for all NBC methods */
      if (1 == IterOrLoca) /* ITERATOR mode */
      {
         ier = (*get_neigh)(&pkim, &zero, &one, &currentAtom, &numOfAtomNeigh,
                             &neighListOfCurrentAtom, &Rij_list);
         if (KIM_STATUS_NEIGH_ITER_PAST_END == ier) /* the end of the list, terminate loop */
         {
            break;
         }
         if (KIM_STATUS_OK > ier) /* some sort of problem, exit */
         {
            KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", ier);
            return ier;
         }

         i = currentAtom + model_index_shift;
      }
      else
      {
         i++;
         if (*nAtoms <= i) /* incremented past end of list, terminate loop */
         {
            break;
         }

         if (3 == NBC) /* CLUSTER NBC method */
         {
            numOfAtomNeigh = *nAtoms - 1;
            for (kdim = 0; kdim < *nAtoms; ++kdim)
            {
               if (kdim < i)
                 neighListOfCurrentAtom[kdim] = kdim - model_index_shift;
               if (kdim > i)
                 neighListOfCurrentAtom[kdim - 1] = kdim - model_index_shift;
            }
            ier = KIM_STATUS_OK;
         }
         else
         {
            request = i - model_index_shift;
            ier = (*get_neigh)(&pkim, &one, &request,
                                &currentAtom, &numOfAtomNeigh,
                                &neighListOfCurrentAtom, &Rij_list);
            if (KIM_STATUS_OK != ier) /* some sort of problem, exit */
            {
               KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
               ier = KIM_STATUS_FAIL;
               return ier;
            }
         }
      }
      iSpecies = particleSpecies[i];
            
      /* loop over the neighbors of atom i */
      for (jj = 0; jj < numOfAtomNeigh; ++ jj)
      {

         j = neighListOfCurrentAtom[jj] + model_index_shift; /* get neighbor ID */
         jSpecies = particleSpecies[j];

         /* get corresponding parameters if *nSpecies == 2*/
         if (*nSpecies == 2)
         {
             if (iSpecies == SPEC1 && jSpecies == SPEC1) interaction_index = 0;
             else if (iSpecies == SPEC2 && jSpecies == SPEC2) interaction_index = 15;
             else interaction_index = 7;
             
             A        = &(buffer->A)[interaction_index];
             B        = &(buffer->B)[interaction_index];
             p        = &(buffer->p)[interaction_index];
             q        = &(buffer->q)[interaction_index];
             a        = &(buffer->a)[interaction_index];
             sigma    = &(buffer->sigma)[interaction_index];
             epsilon  = &(buffer->epsilon)[interaction_index];
             cutsq    = &(buffer->cutsq)[interaction_index];
         }

         /* compute relative position vector and squared distance */
         Rsqij = 0.0;
         for (kdim = 0; kdim < DIM; ++kdim)
         {
            if (0 != NBC) /* all methods except NEIGH_RVEC */
            {
               Rij[kdim] = coords[j*DIM + kdim] - coords[i*DIM + kdim];
            }
            else          /* NEIGH_RVEC method */
            {
               Rij[kdim] = Rij_list[jj*DIM + kdim];
            }

            /* apply periodic boundary conditions if required */
            if (2 == NBC)
            {
               if (fabs(Rij[kdim]) > 0.5*boxSideLengths[kdim])
               {
                  Rij[kdim] -= (Rij[kdim]/fabs(Rij[kdim]))*boxSideLengths[kdim];
               }
            }
            
            /* compute squared distance */
            Rsqij += Rij[kdim]*Rij[kdim];
         }
         
         /* compute energy and force */
         if (Rsqij > *cutsq) continue; /* particles are not interacting  */
         R1 = sqrt(Rsqij);
         if (comp_process_d2Edr2)
         {
            /* compute pair potential and its derivatives */
            calc_phi_d2phi_two(A, B, p, q, a, sigma, epsilon, 
                               R1, &phi_two, &dphi_two, &d2phi_two);
               
            /* compute dEidr */
                 /* Full mode -- regular contribution */
            dEidr_two  = 0.5*dphi_two;
            d2Eidr_two = 0.5*d2phi_two;
        }
        else if (comp_force || comp_process_dEdr)
        {
            /* compute pair potential and its derivative */
            calc_phi_dphi_two(A, B, p, q, a, sigma, epsilon, 
                              R1, &phi_two, &dphi_two);

            /* compute dEidr */
                 /* Full mode -- regular contribution */
            dEidr_two = 0.5*dphi_two;
        }
        else
        {
           /* compute just pair potential */
            calc_phi_two(A, B, p, q, a, sigma, epsilon, 
                         R1, &phi_two);
        }
          
         /* contribution to energy */
        if (comp_particleEnergy)
        {
           particleEnergy[i] += 0.5*phi_two;
        }
        if (comp_energy)
        {
            /* Full mode -- add half v to total energy */
            *energy += 0.5*phi_two;
        }
          
        /* contribution to process_dEdr */
        if (comp_process_dEdr)
        {
           ier = KIM_API_process_dEdr(km, &dEidr_two, &R1, &pRij, &i, &j);
        }
            
        /* contribution to process_d2Edr2 */
        if (comp_process_d2Edr2)
        {
            R_pairs[0] = R_pairs[1] = R1;
            Rij_pairs[0][0] = Rij_pairs[1][0] = Rij[0];
            Rij_pairs[0][1] = Rij_pairs[1][1] = Rij[1];
            Rij_pairs[0][2] = Rij_pairs[1][2] = Rij[2];
            i_pairs[0] = i_pairs[1] = i;
            j_pairs[0] = j_pairs[1] = j;

            ier = KIM_API_process_d2Edr2(km, &d2Eidr_two, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
        } 
          
        /* contribution to forces */
        if (comp_force)
        {
           for (kdim = 0; kdim < DIM; ++kdim)
           {
               force[i*DIM + kdim] += dEidr_two*Rij[kdim]/R1; /* accumulate force on atom i */
               force[j*DIM + kdim] -= dEidr_two*Rij[kdim]/R1; /* accumulate force on atom j */
           }
        }
        
        /* Start adding three body terms */
        /*********************************/
        if(jj == numOfAtomNeigh-1) continue;

        for (kk = jj+1; kk < numOfAtomNeigh; ++kk)
        {
            k = neighListOfCurrentAtom[kk] + model_index_shift; /* get neighbor ID */
            kSpecies = particleSpecies[k];
            
            if (*nSpecies == 2)
            {
                if (iSpecies == SPEC1 && jSpecies == SPEC1 && kSpecies == SPEC1) interaction_index = 0;
                else if (iSpecies == SPEC1 && jSpecies == SPEC1 && kSpecies == SPEC2) interaction_index = 2;
                else if (iSpecies == SPEC1 && jSpecies == SPEC2 && kSpecies == SPEC1) interaction_index = 4;
                else if (iSpecies == SPEC1 && jSpecies == SPEC2 && kSpecies == SPEC2) interaction_index = 6;
                else if (iSpecies == SPEC2 && jSpecies == SPEC1 && kSpecies == SPEC1) interaction_index = 8;
                else if (iSpecies == SPEC2 && jSpecies == SPEC1 && kSpecies == SPEC2) interaction_index = 10;
                else if (iSpecies == SPEC2 && jSpecies == SPEC2 && kSpecies == SPEC1) interaction_index = 12;
                else interaction_index = 14;
                
                a        = &(buffer->a)[interaction_index];
                lambda   = &(buffer->lambda)[interaction_index];
                gamma    = &(buffer->gamma)[interaction_index];
                sigma    = &(buffer->sigma)[interaction_index];
                epsilon  = &(buffer->epsilon)[interaction_index];
                Q        = &(buffer->Q)[interaction_index];
                costhetat= &(buffer->costhetat)[interaction_index];
                cutsq    = &(buffer->cutsq)[interaction_index];
            }
        
            /* compute relative position vector and squared distance */
            Rsqik = 0.0;
            Rsqjk = 0.0;
            for (kdim = 0; kdim < DIM; ++kdim)
            {
               if (0 != NBC) /* all methods except NEIGH_RVEC */
               {
                  Rik[kdim] = coords[k*DIM + kdim] - coords[i*DIM + kdim];
                  Rjk[kdim] = Rik[kdim] - Rij[kdim];
               }
               else          /* NEIGH_RVEC method */
               {
                  Rik[kdim] = Rij_list[kk*DIM + kdim];
                  Rjk[kdim] = Rik[kdim] - Rij[kdim];
               }

               /* apply periodic boundary conditions if required */
               if (2 == NBC)
               {
                  if (fabs(Rik[kdim]) > 0.5*boxSideLengths[kdim])
                  {
                     Rik[kdim] -= (Rik[kdim]/fabs(Rik[kdim]))*boxSideLengths[kdim];
                     Rjk[kdim] = Rik[kdim] - Rij[kdim];
                  } 
                  /*if (fabs(Rjk[kdim]) > 0.5*boxSideLengths[kdim])
                  {
                     Rjk[kdim] -= (Rjk[kdim]/fabs(Rjk[kdim]))*boxSideLengths[kdim];
                  }*/
               }
               
               /* compute squared distance */
               Rsqik += Rik[kdim]*Rik[kdim];
               Rsqjk += Rjk[kdim]*Rjk[kdim];
            }

            /* compute energy and force */
            if (Rsqik > *cutsq) continue; /* particles are interacting ? */
        
            R2 = sqrt(Rsqik);
            R4 = sqrt(Rsqjk);

            if (comp_process_d2Edr2)
            {
               /* compute three-body potential and its derivatives */
               calc_phi_d2phi_three(a, lambda, gamma, sigma, epsilon, Q, costhetat, 
                                    R1, R2, R4, &phi_three, dphi_three, d2phi_three);
               
               /* compute dEidr */
                  /* Full mode -- regular contribution */
               dEidr_three[0]  = dphi_three[0];
               dEidr_three[1]  = dphi_three[1];
               dEidr_three[2]  = dphi_three[2];
               
               d2Eidr_three[0] = d2phi_three[0];
               d2Eidr_three[1] = d2phi_three[1];
               d2Eidr_three[2] = d2phi_three[2];
               d2Eidr_three[3] = d2phi_three[3];
               d2Eidr_three[4] = d2phi_three[4];
               d2Eidr_three[5] = d2phi_three[5];
            }
            else if (comp_force || comp_process_dEdr)
            {
               /* compute three-body potential and its derivative */
               calc_phi_dphi_three(a, lambda, gamma, sigma, epsilon, Q, costhetat, 
                                   R1, R2, R4, &phi_three, dphi_three);

               /* compute dEidr */
                  /* Full mode -- regular contribution */
               dEidr_three[0]  =  dphi_three[0];
               dEidr_three[1]  =  dphi_three[1];
               dEidr_three[2]  =  dphi_three[2];
            }
            else
            {
               /* compute just three-body potential */
               calc_phi_three(a, lambda, gamma, sigma, epsilon, Q, costhetat, 
                             R1, R2, R4, &phi_three);
            }
            
            /* contribution to energy */
            if (comp_particleEnergy)
            {
               particleEnergy[i] += phi_three;
            }
            if (comp_energy)
            {
                  /* Full mode -- add v to total energy */
                *energy += phi_three;
            }
            
            /* contribution to process_dEdr */
            if (comp_process_dEdr)
            {
               ier = KIM_API_process_dEdr(km, dEidr_three, &R1, &pRij, &i, &j);
               ier = KIM_API_process_dEdr(km, dEidr_three+1, &R2, &pRik, &i, &k);
               ier = KIM_API_process_dEdr(km, dEidr_three+2, &R4, &pRjk, &j, &k);
            }
            
            /* contribution to process_d2Edr2 */
            if (comp_process_d2Edr2)
            {
               R_pairs[0] = R_pairs[1] = R1;
               Rij_pairs[0][0] = Rij_pairs[1][0] = Rij[0];
               Rij_pairs[0][1] = Rij_pairs[1][1] = Rij[1];
               Rij_pairs[0][2] = Rij_pairs[1][2] = Rij[2];
               i_pairs[0] = i_pairs[1] = i;
               j_pairs[0] = j_pairs[1] = j;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);

               R_pairs[0] = R_pairs[1] = R2;
               Rij_pairs[0][0] = Rij_pairs[1][0] = Rik[0];
               Rij_pairs[0][1] = Rij_pairs[1][1] = Rik[1];
               Rij_pairs[0][2] = Rij_pairs[1][2] = Rik[2];
               i_pairs[0] = i_pairs[1] = i;
               j_pairs[0] = j_pairs[1] = k;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three+1, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);

               R_pairs[0] = R_pairs[1] = R4;
               Rij_pairs[0][0] = Rij_pairs[1][0] = Rjk[0];
               Rij_pairs[0][1] = Rij_pairs[1][1] = Rjk[1];
               Rij_pairs[0][2] = Rij_pairs[1][2] = Rjk[2];
               i_pairs[0] = i_pairs[1] = j;
               j_pairs[0] = j_pairs[1] = k;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three+2, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);


               R_pairs[0] = R1;
               R_pairs[1] = R2;
               Rij_pairs[0][0] =  Rij[0];
               Rij_pairs[0][1] =  Rij[1];
               Rij_pairs[0][2] =  Rij[2];
               Rij_pairs[1][0] =  Rik[0];
               Rij_pairs[1][1] =  Rik[1];
               Rij_pairs[1][2] =  Rik[2];
               i_pairs[0]  = i;
               j_pairs[0]  = j;
               i_pairs[1]  = i;
               j_pairs[1]  = k;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three+3, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);

               R_pairs[0] = R2;
               R_pairs[1] = R1;
               Rij_pairs[0][0] =  Rik[0];
               Rij_pairs[0][1] =  Rik[1];
               Rij_pairs[0][2] =  Rik[2];
               Rij_pairs[1][0] =  Rij[0];
               Rij_pairs[1][1] =  Rij[1];
               Rij_pairs[1][2] =  Rij[2];
               i_pairs[0]  = i;
               j_pairs[0]  = k;
               i_pairs[1]  = i;
               j_pairs[1]  = j;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three+3, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
               
               R_pairs[0] = R1;
               R_pairs[1] = R4;
               Rij_pairs[0][0] =  Rij[0];
               Rij_pairs[0][1] =  Rij[1];
               Rij_pairs[0][2] =  Rij[2];
               Rij_pairs[1][0] =  Rjk[0];
               Rij_pairs[1][1] =  Rjk[1];
               Rij_pairs[1][2] =  Rjk[2];
               i_pairs[0]  = i;
               j_pairs[0]  = j;
               i_pairs[1]  = j;
               j_pairs[1]  = k;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three+4, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);

               R_pairs[0] = R4;
               R_pairs[1] = R1;
               Rij_pairs[0][0] =  Rjk[0];
               Rij_pairs[0][1] =  Rjk[1];
               Rij_pairs[0][2] =  Rjk[2];
               Rij_pairs[1][0] =  Rij[0];
               Rij_pairs[1][1] =  Rij[1];
               Rij_pairs[1][2] =  Rij[2];
               i_pairs[0]  = j;
               j_pairs[0]  = k;
               i_pairs[1]  = i;
               j_pairs[1]  = j;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three+4, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);

               R_pairs[0] = R2;
               R_pairs[1] = R4;
               Rij_pairs[0][0] =  Rik[0];
               Rij_pairs[0][1] =  Rik[1];
               Rij_pairs[0][2] =  Rik[2];
               Rij_pairs[1][0] =  Rjk[0];
               Rij_pairs[1][1] =  Rjk[1];
               Rij_pairs[1][2] =  Rjk[2];
               i_pairs[0]  = i;
               j_pairs[0]  = k;
               i_pairs[1]  = j;
               j_pairs[1]  = k;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three+5, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);

               R_pairs[0] = R4;
               R_pairs[1] = R2;
               Rij_pairs[0][0] =  Rjk[0];
               Rij_pairs[0][1] =  Rjk[1];
               Rij_pairs[0][2] =  Rjk[2];
               Rij_pairs[1][0] =  Rik[0];
               Rij_pairs[1][1] =  Rik[1];
               Rij_pairs[1][2] =  Rik[2];
               i_pairs[0]  = j;
               j_pairs[0]  = k;
               i_pairs[1]  = i;
               j_pairs[1]  = k;

               ier = KIM_API_process_d2Edr2(km, d2Eidr_three+5, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
            }
            /* contribution to forces */
            if (comp_force)
            {
               for (kdim = 0; kdim < DIM; ++kdim)
               {
                  force[i*DIM + kdim] += dEidr_three[0]*Rij[kdim]/R1 + dEidr_three[1]*Rik[kdim]/R2; /* accumulate force on atom i */
                  force[j*DIM + kdim] -= dEidr_three[0]*Rij[kdim]/R1 - dEidr_three[2]*Rjk[kdim]/R4; /* accumulate force on atom j */
                  force[k*DIM + kdim] -= dEidr_three[2]*Rjk[kdim]/R4 + dEidr_three[1]*Rik[kdim]/R2; /* accumulate force on atom k */
               }
            }

            /* Start adding four body terms */
            /*********************************/
            for (ll = kk+1; ll < numOfAtomNeigh; ++ll)
            {
                l = neighListOfCurrentAtom[ll] + model_index_shift; /* get neighbor ID */
                lSpecies = particleSpecies[l];
                
                if (*nSpecies == 2)
                {
                    if (iSpecies == SPEC1 && jSpecies == SPEC1 && kSpecies == SPEC1 && lSpecies == SPEC1) interaction_index = 0;
                    else if (iSpecies == SPEC1 && jSpecies == SPEC1 && kSpecies == SPEC1 && lSpecies == SPEC2) interaction_index = 1;
                    else if (iSpecies == SPEC1 && jSpecies == SPEC1 && kSpecies == SPEC2 && lSpecies == SPEC1) interaction_index = 2;
                    else if (iSpecies == SPEC1 && jSpecies == SPEC1 && kSpecies == SPEC2 && lSpecies == SPEC2) interaction_index = 3;
                    else if (iSpecies == SPEC1 && jSpecies == SPEC2 && kSpecies == SPEC1 && lSpecies == SPEC1) interaction_index = 4;
                    else if (iSpecies == SPEC1 && jSpecies == SPEC2 && kSpecies == SPEC1 && lSpecies == SPEC2) interaction_index = 5;
                    else if (iSpecies == SPEC1 && jSpecies == SPEC2 && kSpecies == SPEC2 && lSpecies == SPEC1) interaction_index = 6;
                    else if (iSpecies == SPEC1 && jSpecies == SPEC2 && kSpecies == SPEC2 && lSpecies == SPEC2) interaction_index = 7;
                    else if (iSpecies == SPEC2 && jSpecies == SPEC1 && kSpecies == SPEC1 && lSpecies == SPEC1) interaction_index = 8;
                    else if (iSpecies == SPEC2 && jSpecies == SPEC1 && kSpecies == SPEC1 && lSpecies == SPEC2) interaction_index = 9;
                    else if (iSpecies == SPEC2 && jSpecies == SPEC1 && kSpecies == SPEC2 && lSpecies == SPEC1) interaction_index = 10;
                    else if (iSpecies == SPEC2 && jSpecies == SPEC1 && kSpecies == SPEC2 && lSpecies == SPEC2) interaction_index = 11;
                    else if (iSpecies == SPEC2 && jSpecies == SPEC2 && kSpecies == SPEC1 && lSpecies == SPEC1) interaction_index = 12;
                    else if (iSpecies == SPEC2 && jSpecies == SPEC2 && kSpecies == SPEC1 && lSpecies == SPEC2) interaction_index = 13;
                    else if (iSpecies == SPEC2 && jSpecies == SPEC2 && kSpecies == SPEC2 && lSpecies == SPEC1) interaction_index = 14;
                    else interaction_index = 15;
                    
                    a        = &(buffer->a)[interaction_index];
                    lambda_2 = &(buffer->lambda_2)[interaction_index];
                    gamma    = &(buffer->gamma)[interaction_index];
                    sigma    = &(buffer->sigma)[interaction_index];
                    epsilon  = &(buffer->epsilon)[interaction_index];
                    Q        = &(buffer->Q)[interaction_index];
                    costhetat= &(buffer->costhetat)[interaction_index];
                    cutsq    = &(buffer->cutsq)[interaction_index];
                }
        
                /* compute relative position vector and squared distance */
                Rsqil = 0.0;
                Rsqjl = 0.0;
                Rsqkl = 0.0;
                for (kdim = 0; kdim < DIM; ++kdim)
                {
                    if (0 != NBC) /* all methods except NEIGH_RVEC */
                    {
                        Ril[kdim] = coords[l*DIM + kdim] - coords[i*DIM + kdim];
                        Rjl[kdim] = Ril[kdim] - Rij[kdim];
                        Rkl[kdim] = Ril[kdim] - Rik[kdim];
                    }
                    else          /* NEIGH_RVEC method */
                    {
                        Ril[kdim] = Rij_list[ll*DIM + kdim];
                        Rjl[kdim] = Ril[kdim] - Rij[kdim];
                        Rkl[kdim] = Ril[kdim] - Rik[kdim];
                    }

                    /* apply periodic boundary conditions if required */
                    if (2 == NBC)
                    {
                        if (fabs(Ril[kdim]) > 0.5*boxSideLengths[kdim])
                        {
                            Ril[kdim] -= (Ril[kdim]/fabs(Ril[kdim]))*boxSideLengths[kdim];
                            Rjl[kdim] = Ril[kdim] - Rij[kdim];
                            Rkl[kdim] = Ril[kdim] - Rik[kdim];
                        } 
                        /*if (fabs(Rjl[kdim]) > 0.5*boxSideLengths[kdim])
                        {
                            Rjl[kdim] -= (Rjl[kdim]/fabs(Rjl[kdim]))*boxSideLengths[kdim];
                        } 
                        if (fabs(Rkl[kdim]) > 0.5*boxSideLengths[kdim])
                        {
                            Rkl[kdim] -= (Rkl[kdim]/fabs(Rkl[kdim]))*boxSideLengths[kdim];
                        }*/ 
                    }
                    
                    /* compute squared distance */
                    Rsqil += Ril[kdim]*Ril[kdim];
                    Rsqjl += Rjl[kdim]*Rjl[kdim];
                    Rsqkl += Rkl[kdim]*Rkl[kdim];
                }

                /* compute energy and force */
                /* if (Rsqil > *cutsq || Rsqjl > *cutsq || Rsqkl > *cutsq) continue; */ /* particles are interacting ? */
                if (Rsqil > *cutsq) continue; /* particles are interacting ? */
        
                R3 = sqrt(Rsqil);
                R5 = sqrt(Rsqjl);
                R6 = sqrt(Rsqkl);

                if (comp_process_d2Edr2)
                {
                    /* compute four-body potential and its derivatives */
                    calc_phi_d2phi_four(a, lambda_2, gamma, sigma, epsilon, Q, costhetat, 
                                        R1, R2, R3, R4, R5, R6, &phi_four, dphi_four, d2phi_four);
               
                    /* compute dEidr */
                        /* Full mode -- regular contribution */
                    dEidr_four[0]  =  dphi_four[0];
                    dEidr_four[1]  =  dphi_four[1];
                    dEidr_four[2]  =  dphi_four[2];
                    dEidr_four[3]  =  dphi_four[3];
                    dEidr_four[4]  =  dphi_four[4];
                    dEidr_four[5]  =  dphi_four[5];

                    
                    d2Eidr_four[0]  =  d2phi_four[0];
                    d2Eidr_four[1]  =  d2phi_four[1];
                    d2Eidr_four[2]  =  d2phi_four[2];
                    d2Eidr_four[3]  =  d2phi_four[3];
                    d2Eidr_four[4]  =  d2phi_four[4];
                    d2Eidr_four[5]  =  d2phi_four[5];
                    d2Eidr_four[6]  =  d2phi_four[6];
                    d2Eidr_four[7]  =  d2phi_four[7];
                    d2Eidr_four[8]  =  d2phi_four[8];
                    d2Eidr_four[9]  =  d2phi_four[9];
                    d2Eidr_four[10] =  d2phi_four[10];
                    d2Eidr_four[11] =  d2phi_four[11];
                    d2Eidr_four[12] =  d2phi_four[12];
                    d2Eidr_four[13] =  d2phi_four[13];
                    d2Eidr_four[14] =  d2phi_four[14];
                    d2Eidr_four[15] =  d2phi_four[15];
                    d2Eidr_four[16] =  d2phi_four[16];
                    d2Eidr_four[17] =  d2phi_four[17];
                    d2Eidr_four[18] =  d2phi_four[18];
                    d2Eidr_four[19] =  d2phi_four[19];
                    d2Eidr_four[20] =  d2phi_four[20];
                }
                else if (comp_force || comp_process_dEdr)
                {
                    /* compute four-body potential and its derivative */
                    calc_phi_dphi_four(a, lambda_2, gamma, sigma, epsilon, Q, costhetat, 
                                       R1, R2, R3, R4, R5, R6, &phi_four, dphi_four);

                    /* compute dEidr */
                        /* Full mode -- regular contribution */
                    dEidr_four[0]  =  dphi_four[0];
                    dEidr_four[1]  =  dphi_four[1];
                    dEidr_four[2]  =  dphi_four[2];
                    dEidr_four[3]  =  dphi_four[3];
                    dEidr_four[4]  =  dphi_four[4];
                    dEidr_four[5]  =  dphi_four[5];
                }
                else
                {
                    /* compute just four-body potential */
                    calc_phi_four(a, lambda_2, gamma, sigma, epsilon, Q, costhetat, 
                                  R1, R2, R3, R4, R5, R6, &phi_four);
                }
            
                /* contribution to energy */
                if (comp_particleEnergy)
                {
                    particleEnergy[i] += phi_four;
                }
                if (comp_energy)
                {
                        /* Full mode -- add v to total energy */
                    *energy += phi_four;
                }
            
                /* contribution to process_dEdr */
                if (comp_process_dEdr)
                {
                    ier = KIM_API_process_dEdr(km, dEidr_four, &R1, &pRij, &i, &j);
                    ier = KIM_API_process_dEdr(km, dEidr_four+1, &R2, &pRik, &i, &k);
                    ier = KIM_API_process_dEdr(km, dEidr_four+2, &R3, &pRil, &i, &l);
                    ier = KIM_API_process_dEdr(km, dEidr_four+3, &R4, &pRjk, &j, &k);
                    ier = KIM_API_process_dEdr(km, dEidr_four+4, &R5, &pRjl, &j, &l);
                    ier = KIM_API_process_dEdr(km, dEidr_four+5, &R6, &pRkl, &k, &l);
                }
            
                /* contribution to process_d2Edr2 */
                if (comp_process_d2Edr2)
                {
                    R_pairs[0] = R_pairs[1] = R1;                               
                    Rij_pairs[0][0] = Rij_pairs[1][0] = Rij[0];
                    Rij_pairs[0][1] = Rij_pairs[1][1] = Rij[1];
                    Rij_pairs[0][2] = Rij_pairs[1][2] = Rij[2];
                    i_pairs[0] = i_pairs[1] = i;
                    j_pairs[0] = j_pairs[1] = j;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
               
                    R_pairs[0] = R1;
                    R_pairs[1] = R2;                               
                    Rij_pairs[0][0] =  Rij[0];
                    Rij_pairs[0][1] =  Rij[1];
                    Rij_pairs[0][2] =  Rij[2];
                    Rij_pairs[1][0] =  Rik[0];
                    Rij_pairs[1][1] =  Rik[1];
                    Rij_pairs[1][2] =  Rik[2];
                    i_pairs[0]  = i;
                    j_pairs[0]  = j;
                    i_pairs[1]  = i;
                    j_pairs[1]  = k;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+6, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R2;
                    R_pairs[1] = R1;                               
                    Rij_pairs[0][0] =  Rik[0];
                    Rij_pairs[0][1] =  Rik[1];
                    Rij_pairs[0][2] =  Rik[2];
                    Rij_pairs[1][0] =  Rij[0];
                    Rij_pairs[1][1] =  Rij[1];
                    Rij_pairs[1][2] =  Rij[2];
                    i_pairs[0]  = i;
                    j_pairs[0]  = k;
                    i_pairs[1]  = i;
                    j_pairs[1]  = j;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+6, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R1;
                    R_pairs[1] = R3;                               
                    Rij_pairs[0][0] =  Rij[0];
                    Rij_pairs[0][1] =  Rij[1];
                    Rij_pairs[0][2] =  Rij[2];
                    Rij_pairs[1][0] =  Ril[0];
                    Rij_pairs[1][1] =  Ril[1];
                    Rij_pairs[1][2] =  Ril[2];
                    i_pairs[0]  = i;
                    j_pairs[0]  = j;
                    i_pairs[1]  = i;
                    j_pairs[1]  = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+11, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R3;
                    R_pairs[1] = R1;                               
                    Rij_pairs[0][0] =  Ril[0];
                    Rij_pairs[0][1] =  Ril[1];
                    Rij_pairs[0][2] =  Ril[2];
                    Rij_pairs[1][0] =  Rij[0];
                    Rij_pairs[1][1] =  Rij[1];
                    Rij_pairs[1][2] =  Rij[2];
                    i_pairs[0]  = i;
                    j_pairs[0]  = l;
                    i_pairs[1]  = i;
                    j_pairs[1]  = j;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+11, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R1;
                    R_pairs[1] = R4;                               
                    Rij_pairs[0][0] =  Rij[0];
                    Rij_pairs[0][1] =  Rij[1];
                    Rij_pairs[0][2] =  Rij[2];
                    Rij_pairs[1][0] =  Rjk[0];
                    Rij_pairs[1][1] =  Rjk[1];
                    Rij_pairs[1][2] =  Rjk[2];
                    i_pairs[0]  = i;
                    j_pairs[0]  = j;
                    i_pairs[1]  = j;
                    j_pairs[1]  = k;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+15, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R4;
                    R_pairs[1] = R1;                               
                    Rij_pairs[0][0] =  Rjk[0];
                    Rij_pairs[0][1] =  Rjk[1];
                    Rij_pairs[0][2] =  Rjk[2];
                    Rij_pairs[1][0] =  Rij[0];
                    Rij_pairs[1][1] =  Rij[1];
                    Rij_pairs[1][2] =  Rij[2];
                    i_pairs[0]  = j;
                    j_pairs[0]  = k;
                    i_pairs[1]  = i;
                    j_pairs[1]  = j;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+15, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R1;
                    R_pairs[1] = R5;                               
                    Rij_pairs[0][0] =  Rij[0];
                    Rij_pairs[0][1] =  Rij[1];
                    Rij_pairs[0][2] =  Rij[2];
                    Rij_pairs[1][0] =  Rjl[0];
                    Rij_pairs[1][1] =  Rjl[1];
                    Rij_pairs[1][2] =  Rjl[2];
                    i_pairs[0]  = i;
                    j_pairs[0]  = j;
                    i_pairs[1]  = j;
                    j_pairs[1]  = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+18, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R5;
                    R_pairs[1] = R1;                               
                    Rij_pairs[0][0] =  Rjl[0];
                    Rij_pairs[0][1] =  Rjl[1];
                    Rij_pairs[0][2] =  Rjl[2];
                    Rij_pairs[1][0] =  Rij[0];
                    Rij_pairs[1][1] =  Rij[1];
                    Rij_pairs[1][2] =  Rij[2];
                    i_pairs[0]  = j;
                    j_pairs[0]  = l;
                    i_pairs[1]  = i;
                    j_pairs[1]  = j;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+18, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R1;
                    R_pairs[1] = R6;                               
                    Rij_pairs[0][0] =  Rij[0];
                    Rij_pairs[0][1] =  Rij[1];
                    Rij_pairs[0][2] =  Rij[2];
                    Rij_pairs[1][0] =  Rkl[0];
                    Rij_pairs[1][1] =  Rkl[1];
                    Rij_pairs[1][2] =  Rkl[2];
                    i_pairs[0]  = i;
                    j_pairs[0]  = j;
                    i_pairs[1]  = k;
                    j_pairs[1]  = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+20, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R6;
                    R_pairs[1] = R1;                               
                    Rij_pairs[0][0] =  Rkl[0];
                    Rij_pairs[0][1] =  Rkl[1];
                    Rij_pairs[0][2] =  Rkl[2];
                    Rij_pairs[1][0] =  Rij[0];
                    Rij_pairs[1][1] =  Rij[1];
                    Rij_pairs[1][2] =  Rij[2];
                    i_pairs[0]  = k;
                    j_pairs[0]  = l;
                    i_pairs[1]  = i;
                    j_pairs[1]  = j;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+20, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R_pairs[1] = R2;                               
                    Rij_pairs[0][0] = Rij_pairs[1][0] = Rik[0];
                    Rij_pairs[0][1] = Rij_pairs[1][1] = Rik[1];
                    Rij_pairs[0][2] = Rij_pairs[1][2] = Rik[2];
                    i_pairs[0] = i_pairs[1] = i;
                    j_pairs[0] = j_pairs[1] = k;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+1, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
               
                    R_pairs[0] = R2;
                    R_pairs[1] = R3;                               
                    Rij_pairs[0][0] =  Rik[0];
                    Rij_pairs[0][1] =  Rik[1];
                    Rij_pairs[0][2] =  Rik[2];
                    Rij_pairs[1][0] =  Ril[0];
                    Rij_pairs[1][1] =  Ril[1];
                    Rij_pairs[1][2] =  Ril[2];
                    i_pairs[0]  = i;
                    j_pairs[0]  = k;
                    i_pairs[1]  = i;
                    j_pairs[1]  = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+7, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R3;
                    R_pairs[1] = R2;                               
                    Rij_pairs[0][0] =  Ril[0];
                    Rij_pairs[0][1] =  Ril[1];
                    Rij_pairs[0][2] =  Ril[2];
                    Rij_pairs[1][0] =  Rik[0];
                    Rij_pairs[1][1] =  Rik[1];
                    Rij_pairs[1][2] =  Rik[2];
                    i_pairs[0]  = i;
                    j_pairs[0]  = l;
                    i_pairs[1]  = i;
                    j_pairs[1]  = k;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+7, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R2;
                    R_pairs[1] = R4;                               
                    Rij_pairs[0][0] =  Rik[0];
                    Rij_pairs[0][1] =  Rik[1];
                    Rij_pairs[0][2] =  Rik[2];
                    Rij_pairs[1][0] =  Rjk[0];
                    Rij_pairs[1][1] =  Rjk[1];
                    Rij_pairs[1][2] =  Rjk[2];
                    i_pairs[0]  = i;
                    j_pairs[0]  = k;
                    i_pairs[1]  = j;
                    j_pairs[1]  = k;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+12, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R4;
                    R_pairs[1] = R2;                               
                    Rij_pairs[0][0] =  Rjk[0];
                    Rij_pairs[0][1] =  Rjk[1];
                    Rij_pairs[0][2] =  Rjk[2];
                    Rij_pairs[1][0] =  Rik[0];
                    Rij_pairs[1][1] =  Rik[1];
                    Rij_pairs[1][2] =  Rik[2];
                    i_pairs[0]  = j;
                    j_pairs[0]  = k;
                    i_pairs[1]  = i;
                    j_pairs[1]  = k;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+12, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R2;
                    R_pairs[1] = R5;                               
                    Rij_pairs[0][0] =  Rik[0];
                    Rij_pairs[0][1] =  Rik[1];
                    Rij_pairs[0][2] =  Rik[2];
                    Rij_pairs[1][0] =  Rjl[0];
                    Rij_pairs[1][1] =  Rjl[1];
                    Rij_pairs[1][2] =  Rjl[2];
                    i_pairs[0]  = i;
                    j_pairs[0]  = k;
                    i_pairs[1]  = j;
                    j_pairs[1]  = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+16, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R5;
                    R_pairs[1] = R2;                               
                    Rij_pairs[0][0] =  Rjl[0];
                    Rij_pairs[0][1] =  Rjl[1];
                    Rij_pairs[0][2] =  Rjl[2];
                    Rij_pairs[1][0] =  Rik[0];
                    Rij_pairs[1][1] =  Rik[1];
                    Rij_pairs[1][2] =  Rik[2];
                    i_pairs[0]  = j;
                    j_pairs[0]  = l;
                    i_pairs[1]  = i;
                    j_pairs[1]  = k;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+16, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R2;
                    R_pairs[1] = R6;                               
                    Rij_pairs[0][0] =  Rik[0];
                    Rij_pairs[0][1] =  Rik[1];
                    Rij_pairs[0][2] =  Rik[2];
                    Rij_pairs[1][0] =  Rkl[0];
                    Rij_pairs[1][1] =  Rkl[1];
                    Rij_pairs[1][2] =  Rkl[2];
                    i_pairs[0]  = i;
                    j_pairs[0]  = k;
                    i_pairs[1]  = k;
                    j_pairs[1]  = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+19, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R6;
                    R_pairs[1] = R2;                               
                    Rij_pairs[0][0] =  Rkl[0];
                    Rij_pairs[0][1] =  Rkl[1];
                    Rij_pairs[0][2] =  Rkl[2];
                    Rij_pairs[1][0] =  Rik[0];
                    Rij_pairs[1][1] =  Rik[1];
                    Rij_pairs[1][2] =  Rik[2];
                    i_pairs[0]  = k;
                    j_pairs[0]  = l;
                    i_pairs[1]  = i;
                    j_pairs[1]  = k;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+19, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R_pairs[1] = R3;                               
                    Rij_pairs[0][0] = Rij_pairs[1][0] = Ril[0];
                    Rij_pairs[0][1] = Rij_pairs[1][1] = Ril[1];
                    Rij_pairs[0][2] = Rij_pairs[1][2] = Ril[2];
                    i_pairs[0] = i_pairs[1] = i;
                    j_pairs[0] = j_pairs[1] = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+2, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
               
                    R_pairs[0] = R3;
                    R_pairs[1] = R4;                               
                    Rij_pairs[0][0] =  Ril[0];
                    Rij_pairs[0][1] =  Ril[1];
                    Rij_pairs[0][2] =  Ril[2];
                    Rij_pairs[1][0] =  Rjk[0];
                    Rij_pairs[1][1] =  Rjk[1];
                    Rij_pairs[1][2] =  Rjk[2];
                    i_pairs[0]  = i;
                    j_pairs[0]  = l;
                    i_pairs[1]  = j;
                    j_pairs[1]  = k;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+8, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R4;
                    R_pairs[1] = R3;                               
                    Rij_pairs[0][0] =  Rjk[0];
                    Rij_pairs[0][1] =  Rjk[1];
                    Rij_pairs[0][2] =  Rjk[2];
                    Rij_pairs[1][0] =  Ril[0];
                    Rij_pairs[1][1] =  Ril[1];
                    Rij_pairs[1][2] =  Ril[2];
                    i_pairs[0]  = j;
                    j_pairs[0]  = k;
                    i_pairs[1]  = i;
                    j_pairs[1]  = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+8, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R3;
                    R_pairs[1] = R5;                               
                    Rij_pairs[0][0] =  Ril[0];
                    Rij_pairs[0][1] =  Ril[1];
                    Rij_pairs[0][2] =  Ril[2];
                    Rij_pairs[1][0] =  Rjl[0];
                    Rij_pairs[1][1] =  Rjl[1];
                    Rij_pairs[1][2] =  Rjl[2];
                    i_pairs[0]  = i;
                    j_pairs[0]  = l;
                    i_pairs[1]  = j;
                    j_pairs[1]  = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+13, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R5;
                    R_pairs[1] = R3;                               
                    Rij_pairs[0][0] =  Rjl[0];
                    Rij_pairs[0][1] =  Rjl[1];
                    Rij_pairs[0][2] =  Rjl[2];
                    Rij_pairs[1][0] =  Ril[0];
                    Rij_pairs[1][1] =  Ril[1];
                    Rij_pairs[1][2] =  Ril[2];
                    i_pairs[0]  = j;
                    j_pairs[0]  = l;
                    i_pairs[1]  = i;
                    j_pairs[1]  = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+13, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R3;
                    R_pairs[1] = R6;                               
                    Rij_pairs[0][0] =  Ril[0];
                    Rij_pairs[0][1] =  Ril[1];
                    Rij_pairs[0][2] =  Ril[2];
                    Rij_pairs[1][0] =  Rkl[0];
                    Rij_pairs[1][1] =  Rkl[1];
                    Rij_pairs[1][2] =  Rkl[2];
                    i_pairs[0]  = i;
                    j_pairs[0]  = l;
                    i_pairs[1]  = k;
                    j_pairs[1]  = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+17, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R6;
                    R_pairs[1] = R3;                               
                    Rij_pairs[0][0] =  Rkl[0];
                    Rij_pairs[0][1] =  Rkl[1];
                    Rij_pairs[0][2] =  Rkl[2];
                    Rij_pairs[1][0] =  Ril[0];
                    Rij_pairs[1][1] =  Ril[1];
                    Rij_pairs[1][2] =  Ril[2];
                    i_pairs[0]  = k;
                    j_pairs[0]  = l;
                    i_pairs[1]  = i;
                    j_pairs[1]  = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+17, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R_pairs[1] = R4;                               
                    Rij_pairs[0][0] = Rij_pairs[1][0] = Rjk[0];
                    Rij_pairs[0][1] = Rij_pairs[1][1] = Rjk[1];
                    Rij_pairs[0][2] = Rij_pairs[1][2] = Rjk[2];
                    i_pairs[0] = i_pairs[1] = j;
                    j_pairs[0] = j_pairs[1] = k;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+3, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
               
                    R_pairs[0] = R4;
                    R_pairs[1] = R5;                               
                    Rij_pairs[0][0] =  Rjk[0];
                    Rij_pairs[0][1] =  Rjk[1];
                    Rij_pairs[0][2] =  Rjk[2];
                    Rij_pairs[1][0] =  Rjl[0];
                    Rij_pairs[1][1] =  Rjl[1];
                    Rij_pairs[1][2] =  Rjl[2];
                    i_pairs[0]  = j;
                    j_pairs[0]  = k;
                    i_pairs[1]  = j;
                    j_pairs[1]  = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+9, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R5;
                    R_pairs[1] = R4;                               
                    Rij_pairs[0][0] =  Rjl[0];
                    Rij_pairs[0][1] =  Rjl[1];
                    Rij_pairs[0][2] =  Rjl[2];
                    Rij_pairs[1][0] =  Rjk[0];
                    Rij_pairs[1][1] =  Rjk[1];
                    Rij_pairs[1][2] =  Rjk[2];
                    i_pairs[0]  = j;
                    j_pairs[0]  = l;
                    i_pairs[1]  = j;
                    j_pairs[1]  = k;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+9, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R4;
                    R_pairs[1] = R6;                               
                    Rij_pairs[0][0] =  Rjk[0];
                    Rij_pairs[0][1] =  Rjk[1];
                    Rij_pairs[0][2] =  Rjk[2];
                    Rij_pairs[1][0] =  Rkl[0];
                    Rij_pairs[1][1] =  Rkl[1];
                    Rij_pairs[1][2] =  Rkl[2];
                    i_pairs[0]  = j;
                    j_pairs[0]  = k;
                    i_pairs[1]  = k;
                    j_pairs[1]  = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+14, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R6;
                    R_pairs[1] = R4;                               
                    Rij_pairs[0][0] =  Rkl[0];
                    Rij_pairs[0][1] =  Rkl[1];
                    Rij_pairs[0][2] =  Rkl[2];
                    Rij_pairs[1][0] =  Rjk[0];
                    Rij_pairs[1][1] =  Rjk[1];
                    Rij_pairs[1][2] =  Rjk[2];
                    i_pairs[0]  = k;
                    j_pairs[0]  = l;
                    i_pairs[1]  = j;
                    j_pairs[1]  = k;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+14, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R_pairs[1] = R5;                               
                    Rij_pairs[0][0] = Rij_pairs[1][0] = Rjl[0];
                    Rij_pairs[0][1] = Rij_pairs[1][1] = Rjl[1];
                    Rij_pairs[0][2] = Rij_pairs[1][2] = Rjl[2];
                    i_pairs[0] = i_pairs[1] = j;
                    j_pairs[0] = j_pairs[1] = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+4, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
               
                    R_pairs[0] = R5;
                    R_pairs[1] = R6;                               
                    Rij_pairs[0][0] =  Rjl[0];
                    Rij_pairs[0][1] =  Rjl[1];
                    Rij_pairs[0][2] =  Rjl[2];
                    Rij_pairs[1][0] =  Rkl[0];
                    Rij_pairs[1][1] =  Rkl[1];
                    Rij_pairs[1][2] =  Rkl[2];
                    i_pairs[0]  = j;
                    j_pairs[0]  = l;
                    i_pairs[1]  = k;
                    j_pairs[1]  = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+10, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R6;
                    R_pairs[1] = R5;                               
                    Rij_pairs[0][0] =  Rkl[0];
                    Rij_pairs[0][1] =  Rkl[1];
                    Rij_pairs[0][2] =  Rkl[2];
                    Rij_pairs[1][0] =  Rjl[0];
                    Rij_pairs[1][1] =  Rjl[1];
                    Rij_pairs[1][2] =  Rjl[2];
                    i_pairs[0]  = k;
                    j_pairs[0]  = l;
                    i_pairs[1]  = j;
                    j_pairs[1]  = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+10, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
              
                    R_pairs[0] = R_pairs[1] = R6;                               
                    Rij_pairs[0][0] = Rij_pairs[1][0] = Rkl[0];
                    Rij_pairs[0][1] = Rij_pairs[1][1] = Rkl[1];
                    Rij_pairs[0][2] = Rij_pairs[1][2] = Rkl[2];
                    i_pairs[0] = i_pairs[1] = k;
                    j_pairs[0] = j_pairs[1] = l;

                    ier = KIM_API_process_d2Edr2(km, d2Eidr_four+5, &pR_pairs, &pRij_pairs, &pi_pairs, &pj_pairs);
                } 
                /* contribution to forces */
                if (comp_force)
                {
                    for (kdim = 0; kdim < DIM; ++kdim)
                    {
                        force[i*DIM + kdim] +=  dEidr_four[0]*Rij[kdim]/R1 + dEidr_four[1]*Rik[kdim]/R2 + dEidr_four[2]*Ril[kdim]/R3; /* accumulate force on atom i */
                        force[j*DIM + kdim] += -dEidr_four[0]*Rij[kdim]/R1 + dEidr_four[3]*Rjk[kdim]/R4 + dEidr_four[4]*Rjl[kdim]/R5; /* accumulate force on atom j */
                        force[k*DIM + kdim] += -dEidr_four[1]*Rik[kdim]/R2 - dEidr_four[3]*Rjk[kdim]/R4 + dEidr_four[5]*Rkl[kdim]/R6; /* accumulate force on atom k */
                        force[l*DIM + kdim] += -dEidr_four[2]*Ril[kdim]/R3 - dEidr_four[4]*Rjl[kdim]/R5 - dEidr_four[5]*Rkl[kdim]/R6; /* accumulate force on atom l */
                    }    
                }
            } /* loop on ll */
            /* End adding four body terms */
            /******************************/
        } /* loop on kk */

        /* End adding three body terms */
        /*******************************/

      } /* loop on jj */

   }    /* infinite while loop terminated by break statements above */

   /* perform final tasks */
   
   if (comp_particleEnergy && comp_energy)
   {
      *energy = 0.0;
      for (k = 0; k < *nAtoms; ++k)
      {
         *energy += particleEnergy[k];
      }
   }

   /* Free temporary storage */
   if (3 == NBC) 
   {
      free(neighListOfCurrentAtom);
   }

   free(dphi_three);
   free(d2phi_three);
   free(dEidr_three);
   free(d2Eidr_three);
   
   free(dphi_four);
   free(d2phi_four);
   free(dEidr_four);
   free(d2Eidr_four);
   
   /* everything is great */
   ier = KIM_STATUS_OK;
   
   return ier;
}

/* Initialization function */
int model_driver_init(void *km, char* paramfile_names, int* nmstrlen, int* numparamfiles)
{
   /* KIM variables */
   intptr_t* pkim = *((intptr_t**) km);
   char* paramfile1name;

   /* Local variables */
   FILE* fid;
   double* A;
   double* B;
   double* p;
   double* q;
   double* a;
   double* lambda;
   double* lambda_2;
   double* gamma;
   double* sigma;
   double* epsilon;
   double* Q;
   double* costhetat;
   /* double <FILL parameter 1>; */
   double* model_Pcutoff;
   double* model_cutoff;
   double* model_cutsq;
   double* model_A;
   double* model_B;
   double* model_p;
   double* model_q;
   double* model_a;
   double* model_lambda;
   double* model_lambda_2;
   double* model_gamma;
   double* model_sigma;
   double* model_epsilon;
   double* model_Q;
   double* model_costhetat;
   /* double* model_<FILL parameter 1>; */
   /* double* model_<FILL parameter 2>; */
   /* FILL as many parameters as needed */
   int ier;
   struct model_buffer* buffer;
   const char* NBCstr;

   int num_species;
   int num_interactions;
   fpos_t filepos;
   char dummy[255];
   int i;
   
   /* set paramfile1name */
   if (*numparamfiles != 1)
   {
       ier = KIM_STATUS_FAIL;
       KIM_API_report_error(__LINE__, __FILE__, "Incorrect number of parameter files.", ier);
       return ier;
   }
   paramfile1name = paramfile_names;

   /* store pointer to functions in KIM object */
   KIM_API_setm_method(pkim, &ier, 3*4,
                     "compute", 1, &compute, 1,
                     "reinit",  1, &reinit,  1,
                     "destroy", 1, &destroy, 1);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_setm_method", ier);
      return ier;
   }
   
   /* Read in model parameters from parameter file */
   fid = fopen(paramfile1name, "r");
   if (fid == NULL)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "Unable to open parameter file for Four_Body_Mistriotis_Flytzanis_Farantos parameters", ier);
      return ier;
   }
   
   /* get rid of initial comments in the parameter file */
   fgetpos(fid, &filepos);
   fgets(dummy, 255, fid);
   while (dummy[0] == '#' || isspace(dummy[0])) {
     fgetpos(fid, &filepos);
     fgets(dummy, 255, fid);
   }
   fsetpos(fid, &filepos);
   /* get rid of comments end */
   
   /* read number of species */  
   ier = fscanf(fid, "%d\n", &num_species);
   /* read number of species end */
   if(!(num_species == 1 || num_species == 2))
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "the parameter file must have either 1 or 2 number of species", ier);
      return ier;
   }
   if (ier != 1)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "error reading first line of parameter file", ier);
      return ier;
   }
   
   /* For two species in four body system, we have num_interactions = \sum_{i=0}^{i=4} 4Ci * (4-i)C(4-i) = 2^4 */
   if (num_species == 1) num_interactions = 1;
   else num_interactions = 16;
   
   /* Allocate memory for local parameters */
   A         = (double*) malloc(num_interactions*sizeof(double));  
   B         = (double*) malloc(num_interactions*sizeof(double));  
   p         = (double*) malloc(num_interactions*sizeof(double));  
   q         = (double*) malloc(num_interactions*sizeof(double));  
   a         = (double*) malloc(num_interactions*sizeof(double));  
   lambda    = (double*) malloc(num_interactions*sizeof(double));  
   lambda_2  = (double*) malloc(num_interactions*sizeof(double));  
   gamma     = (double*) malloc(num_interactions*sizeof(double));  
   sigma     = (double*) malloc(num_interactions*sizeof(double));  
   epsilon   = (double*) malloc(num_interactions*sizeof(double));  
   Q         = (double*) malloc(num_interactions*sizeof(double));  
   costhetat = (double*) malloc(num_interactions*sizeof(double));  
   
   if(A==NULL || B==NULL || p==NULL || q==NULL || a==NULL 
           || lambda==NULL || lambda_2==NULL || gamma==NULL 
           || sigma==NULL || epsilon==NULL || Q ==NULL || costhetat==NULL)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
      return ier;
   }
   
   /* read parameters */
   for (i=0; i< num_interactions; ++i) {
     /* get rid of comments begin */
     fgetpos(fid, &filepos);
     fgets(dummy, 255, fid);
     while (dummy[0] == '#' || isspace(dummy[0])) {
       fgetpos(fid, &filepos);
       fgets(dummy, 255, fid);
     }
     fsetpos(fid, &filepos);
     /* get rid of comments end */ 

     ier = fscanf(fid, "%lf\n %lf\n %lf\n %lf\n %lf\n %lf\n %lf\n %lf\n %lf\n %lf\n %lf\n %lf\n",
                  &A[i],
                  &B[i],
                  &p[i],
                  &q[i],
                  &a[i],
                  &lambda[i],
                  &lambda_2[i],
                  &gamma[i],
                  &sigma[i],
                  &epsilon[i],
                  &Q[i],
                  &costhetat[i]);

     /* check that we read the right number of parameters */
     if (12 != ier)
     {
         ier = KIM_STATUS_FAIL;
         KIM_API_report_error(__LINE__, __FILE__, "Unable to read all Four_Body_Mistriotis_Flytzanis_Farantos parameters", ier);
         return ier;
     }
   }
   fclose(fid);
   
   /* convert to appropriate units */
   for (i=0; i< num_interactions; ++i) {
       a[i] *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "ps",
                                               0.0, 0.0,  0.0, 0.0, 0.0, &ier);
       if (KIM_STATUS_OK > ier)
       {
           KIM_API_report_error(__LINE__, __FILE__,"KIM_API_convert_to_act_unit", ier);
      return ier;
       }
       sigma[i] *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "ps",
                                               1.0, 0.0,  0.0, 0.0, 0.0, &ier);
       if (KIM_STATUS_OK > ier)
       {
           KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
           return ier;
       }
       epsilon[i] *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "ps",
                                               0.0, 1.0,  0.0, 0.0, 0.0, &ier);
       if (KIM_STATUS_OK > ier)
       {
           KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
           return ier;
       }
   }
   
   /* store model cutoff in KIM object */
   model_cutoff = (double*) KIM_API_get_data(pkim, "cutoff", &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", ier);
      return ier;
   }
   
   /* allocate memory for parameters */
   model_Pcutoff   = (double*) malloc(1*sizeof(double));
   model_cutsq     = (double*) malloc(num_interactions*sizeof(double));  
   model_A         = (double*) malloc(num_interactions*sizeof(double));  
   model_B         = (double*) malloc(num_interactions*sizeof(double));  
   model_p         = (double*) malloc(num_interactions*sizeof(double));  
   model_q         = (double*) malloc(num_interactions*sizeof(double));  
   model_a         = (double*) malloc(num_interactions*sizeof(double));  
   model_lambda    = (double*) malloc(num_interactions*sizeof(double));  
   model_lambda_2  = (double*) malloc(num_interactions*sizeof(double));  
   model_gamma     = (double*) malloc(num_interactions*sizeof(double));  
   model_sigma     = (double*) malloc(num_interactions*sizeof(double));  
   model_epsilon   = (double*) malloc(num_interactions*sizeof(double));  
   model_Q         = (double*) malloc(num_interactions*sizeof(double));  
   model_costhetat = (double*) malloc(num_interactions*sizeof(double));  
   
   if( model_Pcutoff   ==NULL
     || model_cutsq    ==NULL
     || model_A        ==NULL
     || model_B        ==NULL
     || model_p        ==NULL
     || model_q        ==NULL
     || model_a        ==NULL
     || model_lambda   ==NULL
     || model_lambda_2 ==NULL
     || model_gamma    ==NULL
     || model_sigma    ==NULL
     || model_epsilon  ==NULL
     || model_Q        ==NULL
     || model_costhetat==NULL )
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
      return ier;
   }

   /* store parameters in KIM object */
   KIM_API_setm_data(pkim, &ier, 14*4,
                             "PARAM_FREE_cutoff",   1,  model_cutoff,              1,
                             "PARAM_FIXED_cutsq",  num_interactions,  model_cutsq,                1,
                             "PARAM_FREE_A",         num_interactions,  model_A,                    1,
                             "PARAM_FREE_B",         num_interactions,  model_B,                    1,
                             "PARAM_FREE_p",         num_interactions,  model_p,                    1,
                             "PARAM_FREE_q",         num_interactions,  model_q,                    1,
                             "PARAM_FREE_a",         num_interactions,  model_a,                    1,
                             "PARAM_FREE_lambda",    num_interactions,  model_lambda,               1,
                             "PARAM_FREE_lambda_2",  num_interactions,  model_lambda_2,             1,
                             "PARAM_FREE_gamma",     num_interactions,  model_gamma,                1,
                             "PARAM_FREE_sigma",     num_interactions,  model_sigma,                1,
                             "PARAM_FREE_epsilon",   num_interactions,  model_epsilon,              1,
                             "PARAM_FREE_Q",         num_interactions,  model_Q,                    1,
                             "PARAM_FREE_costhetat",  num_interactions,  model_costhetat,           1
                            );
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_setm_data", ier);
      return ier;
   }
   
   /* set value of parameters */
   for (i=0; i< num_interactions; ++i) {
       model_A[i] = A[i];
       model_B[i] = B[i];
       model_p[i] = p[i];
       model_q[i] = q[i];
       model_a[i] = a[i];
       model_lambda[i] = lambda[i];
       model_lambda_2[i] = lambda_2[i];
       model_gamma[i] = gamma[i];
       model_sigma[i] = sigma[i];
       model_epsilon[i] = epsilon[i];
       model_Q[i] = Q[i];
       model_costhetat[i] = costhetat[i];
       *model_cutoff = a[i]*sigma[i];
       model_cutsq[i] = (*model_cutoff)*(*model_cutoff);
   }
   *model_Pcutoff = *model_cutoff;
   
   /* allocate buffer */
   buffer = (struct model_buffer*) malloc(sizeof(struct model_buffer));
   if (NULL == buffer)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
      return ier;
   }

   /* setup buffer */
   /* Determine neighbor list boundary condition (NBC) */
   ier = KIM_API_get_NBC_method(pkim, &NBCstr);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_NBC_method", ier);
      return ier;
   }
   if (!strcmp("NEIGH_RVEC_F",NBCstr))
   {
      buffer->NBC = 0;
   }
   else if (!strcmp("NEIGH_PURE_F",NBCstr))
   {
      buffer->NBC = 1;
   }
   else if (!strcmp("MI_OPBC_F",NBCstr))
   {
      buffer->NBC = 2;
   }
   else if (!strcmp("CLUSTER",NBCstr))
   {
      buffer->NBC = 3;
   }
   else
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "Unknown NBC method", ier);
      return ier;
   }

   /* determine neighbor list handling mode */
   if (buffer->NBC != 3)
   {
      /*****************************
       * IterOrLoca = 1 -- Iterator
       *            = 2 -- Locator
       *****************************/
      buffer->IterOrLoca = KIM_API_get_neigh_mode(pkim, &ier);
      if (KIM_STATUS_OK > ier)
      {
         KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh_mode", ier);
         return ier;
      }
      if ((buffer->IterOrLoca != 1) && (buffer->IterOrLoca != 2))
      {
         printf("* ERROR: Unsupported IterOrLoca mode = %i\n", buffer->IterOrLoca);
         return ier;
      }
   }
   else
   {
      buffer->IterOrLoca = 2;   /* for CLUSTER NBC */
   }

   buffer->model_index_shift = KIM_API_get_model_index_shift(pkim);

   KIM_API_getm_index(pkim, &ier, 12*3,
                      "cutoff",                      &(buffer->cutoff_ind),                      1,
                      "numberOfParticles",           &(buffer->numberOfParticles_ind),           1,
                      "numberOfSpecies",             &(buffer->numberOfSpecies_ind),             1,
                      "particleSpecies",             &(buffer->particleSpecies_ind),               1,
                      "coordinates",                 &(buffer->coordinates_ind),                 1,
                      "get_neigh",                   &(buffer->get_neigh_ind),                   1,
                      "boxSideLengths",              &(buffer->boxSideLengths_ind),              1,
                      "energy",                      &(buffer->energy_ind),                      1,
                      "forces",                      &(buffer->forces_ind),                      1,
                      "particleEnergy",              &(buffer->particleEnergy_ind),              1,
                      "process_dEdr",                &(buffer->process_dEdr_ind),                1,
                      "process_d2Edr2",              &(buffer->process_d2Edr2_ind),              1
                     );
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_index", ier);
      return ier;
   }

   /* store parameters in buffer */
   buffer->Pcutoff    = model_Pcutoff;
   buffer->cutsq      = model_cutsq;
   buffer->A          = model_A;
   buffer->B          = model_B;
   buffer->p          = model_p;
   buffer->q          = model_q;
   buffer->a          = model_a;
   buffer->lambda     = model_lambda;
   buffer->lambda_2   = model_lambda_2;
   buffer->gamma      = model_gamma;
   buffer->sigma      = model_sigma;
   buffer->epsilon    = model_epsilon;
   buffer->Q          = model_Q;
   buffer->costhetat  = model_costhetat;
     
   /* store in model buffer */
   KIM_API_set_model_buffer(pkim, (void*) buffer, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_set_model_buffer", ier);
      return ier;
   }
   
   free(A);        
   free(B);        
   free(p);        
   free(q);        
   free(a);        
   free(lambda);   
   free(lambda_2); 
   free(gamma);    
   free(epsilon);    
   free(Q);    
   free(costhetat);    

   ier = KIM_STATUS_OK;
   return ier;
}

/* Reinitialization function */
static int reinit(void *km)
{
    /* Local variables */
   intptr_t* pkim = *((intptr_t**) km);
   int ier;
   double *cutoff;
   struct model_buffer* buffer;
   int* nSpecies;
   int num_interactions;
   int i;

   /* get buffer from KIM object */
   buffer = (struct model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
      return ier;
   }
   
   /* For two species in four body system, we have  num_interactions = \sum_{i=0}^{i=4} 4Ci * (4-i)C(4-i) = 2^4*/
   nSpecies = KIM_API_get_data_by_index(pkim, buffer->numberOfSpecies_ind, &ier);

   if (*nSpecies == 1) num_interactions = 1;
   else if (*nSpecies == 2) num_interactions = 16;
   else printf("The number of species should be either 1 or 2.\n"); 

   /* set new values in KIM object     */
   /*                                  */
   /* update cutoff in KIM API and also cutsq */
   cutoff = KIM_API_get_data(pkim, "cutoff", &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", ier);
      return ier;
   }
   *cutoff = *buffer->Pcutoff;
   
   /* set value of parameter cutsq */
   for (i=0; i< num_interactions; ++i) 
       buffer->cutsq[i] = (*cutoff)*(*cutoff);

   /* FILL: store any other FIXED parameters whose values depend on FREE parameters */

   ier = KIM_STATUS_OK;
   return ier;
}

/* destroy function */
static int destroy(void *km)
{
   /* Local variables */
   intptr_t* pkim = *((intptr_t**) km);
   struct model_buffer* buffer;
   int ier;

   /* get model buffer from KIM object */
   buffer = (struct model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
      return ier;
   }

   /* free parameters */
   free(buffer->Pcutoff);
   free(buffer->cutsq);
   free(buffer->A);
   free(buffer->B);
   free(buffer->p);
   free(buffer->q);
   free(buffer->a);
   free(buffer->lambda);
   free(buffer->lambda_2);
   free(buffer->gamma);
   free(buffer->sigma);
   free(buffer->epsilon);
   free(buffer->Q);
   free(buffer->costhetat);
   /* FILL: repeat above statements as many times as necessary for all FREE and FIXED parameters. */

   /* destroy the buffer */
   free(buffer);

   ier = KIM_STATUS_OK;
   return ier;
}
