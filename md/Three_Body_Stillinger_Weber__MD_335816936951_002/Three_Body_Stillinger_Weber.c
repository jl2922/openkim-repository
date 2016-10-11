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
*    Ryan S. Elliott
*    Ellad B. Tadmor
*    Valeriu Smirichinski
*    Amit Singh
*    Mingjian Wen
*/

/*******************************************************************************
*
*  Three_Body_Stillinger_Weber
*
*  Stillinger-Weber type Three-Body potential KIM Model Driver
*
*  Language: C
*
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "KIM_API_C.h"
#include "KIM_API_status.h"

#define DIM 3       /* dimensionality of space */

/******************************************************************************
* Below are the definitions and values of all Model parameters
*******************************************************************************/

/* Define prototypes for Model Driver init */
/* must be all lowercase to be compatible with the KIM API (to support Fortran Tests) */
int model_driver_init(void* km, char* paramfile_names, int* nmstrlen, int* numparamfiles);

/* Define prototypes for Model (Driver) reinit, compute, and destroy */
/* defined as static to avoid namespace clashes with other Models    */
static int reinit(void* km);
static int destroy(void* km);
static int compute(void* km);

/* local function prototypes */
static int get_param_index(int i, int j, int num_species);

static void calc_phi_two(double A, double B, double p, double q, double cutoff,
                         double sigma, double r, double* phi);

static void calc_phi_dphi_two(double A, double B, double p, double q,
                              double cutoff, double sigma, double r,
                              double* phi, double* dphi);

static void calc_phi_d2phi_two(double A, double B, double p, double q,
                               double cutoff, double sigma, double r,
                               double* phi, double* dphi, double* d2phi);

static void calc_phi_three(double cutoff_ij,double cutoff_ik, double lambda,
                           double gamma_ij, double gamma_ik, double costheta, 
                           double rij, double rik, double rjk, double* phi);

static void calc_phi_dphi_three(double cutoff_ij,double cutoff_ik, double lambda,
                                double gamma_ij, double gamma_ik, double costheta, 
                                double rij, double rik, double rjk,
                                double* phi, double* dphi); 

static void calc_phi_d2phi_three(double cutoff_ij, double cutoff_ik, double lambda,
                                 double gamma_ij, double gamma_ik, double costheta,
                                 double rij, double rik, double rjk, double* phi,
                                 double* dphi, double* d2phi);

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
   int num_model_species;
   
   int* sim_species_code; 

   double* cutoff;
   double* cutsq;
   double* A;
   double* B;
   double* p;
   double* q;
   double* lambda;
   double* gamma;
   double* sigma;
   double* costheta;
};

/* Calculate pair potential phi_two(r) */
static void calc_phi_two(double A, double B, double p, double q, double cutoff,
                         double sigma, double r, double* phi)
{
  /* Local variables */
   double r_cap;

   r_cap = r/sigma;

   if (r >= cutoff)
   {
      *phi = 0.0;
   }
   else
   {
      *phi = A * (B*pow(r_cap,-p) - pow(r_cap,-q)) * exp(sigma/(r - cutoff));
   }

   return;
}

/* Calculate pair potential phi_two(r) and its derivative dphi_two(r) */
static void calc_phi_dphi_two(double A, double B, double p, double q,
                              double cutoff, double sigma, double r,
                              double* phi, double* dphi)
{
   /* Local variables */
   double r_cap;

   r_cap = r/sigma;

   if (r >= cutoff)
   {
      *phi = 0.0;
      *dphi = 0.0;
   }
   else
   {
      *phi = A * (B*pow(r_cap,-p) - pow(r_cap,-q)) * exp(sigma/(r - cutoff));

      *dphi = (q*pow(r_cap,-(q+1)) - p*B*pow(r_cap,-(p+1)))
            - (B*pow(r_cap,-p) - pow(r_cap,-q)) * pow((r - cutoff)/sigma, -2);
      *dphi *= (1/sigma) * A * exp(sigma/(r - cutoff));
   }

   return;
}

/* Calculate pair potential phi_two(r) and its 1st & 2nd derivatives dphi_two(r), d2phi_two(r) */
static void calc_phi_d2phi_two(double A, double B, double p, double q,
                               double cutoff, double sigma, double r,
                               double* phi, double* dphi, double* d2phi)
{
   /* Local variables */
   double r_cap;

   r_cap = r/sigma;

   if (r >= cutoff)
   {
      *phi = 0.0;
      *dphi = 0.0;
      *d2phi = 0.0;
   }
   else
   {
      *phi = A * (B*pow(r_cap,-p) - pow(r_cap,-q)) * exp(sigma/(r - cutoff));

      *dphi = (q*pow(r_cap,-(q+1)) - p*B*pow(r_cap,-(p+1)))
              - (B*pow(r_cap,-p) - pow(r_cap,-q)) * pow((r - cutoff)/sigma, -2);
      *dphi *= (1/sigma) * A * exp(sigma/(r - cutoff));

      *d2phi = (B*pow(r_cap,-p) - pow(r_cap,-q)) * 
               (pow((r - cutoff)/sigma, -4) + 2*pow((r - cutoff)/sigma, -3))
               + 2*(p*B*pow(r_cap,-(p+1)) - q*pow(r_cap,-(q+1)))
               * pow((r - cutoff)/sigma, -2)
               + (p*(p+1)*B*pow(r_cap, -(p+2)) - q*(q+1)*pow(r_cap, -(q+2)));
      *d2phi *= (1 / (sigma*sigma)) * A * exp(sigma/(r-cutoff));
   }

   return;

}

/* Calculate pair potential phi_three(rij, rik, rjk) */
static void calc_phi_three(double cutoff_ij,double cutoff_ik, double lambda,
                           double gamma_ij, double gamma_ik, double costheta, 
                           double rij, double rik, double rjk, double* phi)
{
   /* local variables */
   double costhetajik;
   double diff_costhetajik;
   double exp_ij_ik;

   if ((rij < cutoff_ij) && (rik < cutoff_ik))
   {
     costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);
     diff_costhetajik = costhetajik - costheta;
     exp_ij_ik = exp(gamma_ij/(rij - cutoff_ij) + gamma_ik/(rik - cutoff_ik));

     *phi  = lambda * exp_ij_ik * diff_costhetajik *  diff_costhetajik;
   } else
   {
     *phi = 0.0;

   }
   return;
}

/* Calculate pair potential phi_three(rij, rik, rjk) and its 1st derivative
  dphi_three(rij, rik, rjk)

  dphi has three components as derivatives of phi w.r.t. rij, rik, rjk
*/
static void calc_phi_dphi_three(double cutoff_ij, double cutoff_ik, double lambda,
                                double gamma_ij, double gamma_ik, double costheta,
                                double rij, double rik, double rjk, double* phi, double* dphi)
{
   double costhetajik;
   double diff_costhetajik;
   double costhetajik_ij;
   double costhetajik_ik;
   double costhetajik_jk;
   double exp_ij_ik;
   double d_ij;
   double d_ik;

   if ((rij < cutoff_ij) && (rik < cutoff_ik))
   {
     costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);

     diff_costhetajik = costhetajik - costheta;

     /* Derivatives of cosines w.r.t rij, rik, rjk */
     costhetajik_ij = (pow(rij,2) - pow(rik,2) + pow(rjk,2))/(2*rij*rij*rik);
     costhetajik_ik = (pow(rik,2) - pow(rij,2) + pow(rjk,2))/(2*rij*rik*rik);
     costhetajik_jk = -rjk/(rij*rik);

     /* Variables for simplifying terms */
     exp_ij_ik = exp(gamma_ij/(rij - cutoff_ij) + gamma_ik/(rik - cutoff_ik));

     d_ij = -gamma_ij*pow(rij - cutoff_ij, -2);
     d_ik = -gamma_ik*pow(rik - cutoff_ik, -2);

     *phi    = lambda * exp_ij_ik * diff_costhetajik *  diff_costhetajik;
     dphi[0] = lambda * diff_costhetajik * exp_ij_ik * (d_ij * diff_costhetajik + 2*costhetajik_ij);
     dphi[1] = lambda * diff_costhetajik * exp_ij_ik * (d_ik * diff_costhetajik + 2*costhetajik_ik);
     dphi[2] = lambda * diff_costhetajik * exp_ij_ik * 2 * costhetajik_jk;
   } else
   {
     *phi    = 0.0;
     dphi[0] = 0.0;
     dphi[1] = 0.0;
     dphi[2] = 0.0;
   }
   return;
}

/* Calculate pair potential phi_three(rij, rik, rjk) and its 1st & 2nd derivatives
  dphi_three(rij, rik, rjk), d2phi_three(rij, rik, rjk)

  dphi has three components as derivatives of phi w.r.t. rij, rik, rjk

  d2phi as symmetric Hessian matrix of phi has six components: [0]=(ij,ij), [3]=(ij,ik), [4]=(ij,jk)
                                                                            [1]=(ik,ik), [5]=(ik,jk)
                                                                                         [2]=(jk,jk)
*/
static void calc_phi_d2phi_three(double cutoff_ij, double cutoff_ik, double lambda,
                                 double gamma_ij, double gamma_ik, double costheta,
                                 double rij, double rik, double rjk,
                                 double* phi, double* dphi, double* d2phi)
{
    /* local variables */
   double costhetajik;
   double diff_costhetajik;
   double diff_costhetajik_2;
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

   if ((rij < cutoff_ij) && (rik < cutoff_ik))
   {
     costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);
     diff_costhetajik = costhetajik - costheta;
     diff_costhetajik_2 = diff_costhetajik * diff_costhetajik;

     /* Derivatives of cosines w.r.t. r_ij, r_ik, r_jk */
     costhetajik_ij = (pow(rij,2) - pow(rik,2) + pow(rjk,2))/(2*rij*rij*rik);
     costhetajik_ik = (pow(rik,2) - pow(rij,2) + pow(rjk,2))/(2*rij*rik*rik);
     costhetajik_jk = -rjk/(rij*rik);

     /* Hessian matrix of cosine */
     costhetajik_ij_ij = (pow(rik,2) - pow(rjk,2))/(rij*rij*rij*rik);
     costhetajik_ik_ik = (pow(rij,2) - pow(rjk,2))/(rij*rik*rik*rik);
     costhetajik_jk_jk = -1/(rij*rik);
     costhetajik_ij_ik = -(pow(rij,2) + pow(rik,2) + pow(rjk,2))/(2*rij*rij*rik*rik);
     costhetajik_ij_jk = rjk/(rij*rij*rik);
     costhetajik_ik_jk = rjk/(rik*rik*rij);

     /* Variables for simplifying terms */
     exp_ij_ik = exp(gamma_ij/(rij - cutoff_ij) + gamma_ik/(rik - cutoff_ik));

     d_ij = -gamma_ij*pow(rij - cutoff_ij, -2);
     d_ik = -gamma_ik*pow(rik - cutoff_ik, -2);
     d_ij_2 = d_ij * d_ij;
     d_ik_2 = d_ik * d_ik;
     dd_ij = 2*gamma_ij*pow(rij - cutoff_ij, -3);
     dd_ik = 2*gamma_ik*pow(rik - cutoff_ik, -3);


     *phi    = lambda * exp_ij_ik * diff_costhetajik *  diff_costhetajik;
     dphi[0] = lambda * diff_costhetajik * exp_ij_ik * (d_ij * diff_costhetajik + 2*costhetajik_ij);
     dphi[1] = lambda * diff_costhetajik * exp_ij_ik * (d_ik * diff_costhetajik + 2*costhetajik_ik);
     dphi[2] = lambda * diff_costhetajik * exp_ij_ik * 2 * costhetajik_jk;

     d2phi[0] = exp_ij_ik * ((d_ij_2 + dd_ij) * diff_costhetajik_2 +
                  (4 * d_ij * costhetajik_ij + 2 * costhetajik_ij_ij) *
                  diff_costhetajik + 2 * costhetajik_ij * costhetajik_ij);
     d2phi[1] = exp_ij_ik * ((d_ik_2 + dd_ik) * diff_costhetajik_2 +
                  (4 * d_ik * costhetajik_ik + 2 * costhetajik_ik_ik) *
                  diff_costhetajik + 2 * costhetajik_ik * costhetajik_ik);
     d2phi[2] = 2 * exp_ij_ik * ( costhetajik_jk_jk * diff_costhetajik + costhetajik_jk * costhetajik_jk);
     d2phi[3] = exp_ij_ik * (d_ij * d_ik * diff_costhetajik_2 +
                 (d_ij * costhetajik_ik + d_ik * costhetajik_ij + costhetajik_ij_ik) *
                 2 * diff_costhetajik + 2 * costhetajik_ij * costhetajik_ik);
     d2phi[4] = exp_ij_ik * ((d_ij * costhetajik_jk + costhetajik_ij_jk) * 
                 2 * diff_costhetajik + 2 * costhetajik_ij * costhetajik_jk);
     d2phi[5] = exp_ij_ik * ((d_ik * costhetajik_jk + costhetajik_ik_jk) *
                 2 * diff_costhetajik + 2 * costhetajik_ik * costhetajik_jk);

     d2phi[0]  *= lambda;   /*derivative is w.r.t. rij, rij*/
     d2phi[1]  *= lambda;   /*derivative is w.r.t. rik, rik*/
     d2phi[2]  *= lambda;   /*derivative is w.r.t. rjk, rjk*/ 
     d2phi[3]  *= lambda;   /*derivative is w.r.t. rij, rik*/ 
     d2phi[4]  *= lambda;   /*derivative is w.r.t. rij, rjk*/
     d2phi[5]  *= lambda;   /*derivative is w.r.t. rik, rjk*/ 
   } else
   {
     *phi = 0.0;
     dphi[0]  = dphi[1] = dphi[2] = 0.0;
     d2phi[0] = d2phi[1] = d2phi[2] = 0.0;
     d2phi[3] = d2phi[4] = d2phi[5] = 0.0;
   }

   return;
}

/* compute function */
static int compute(void* km)
{
   /* local variables */
   intptr_t* pkim = *((intptr_t**) km);
   double R1;
   double R2;
   double R3;
   double R_pairs[2];
   double *pR_pairs = &(R_pairs[0]);
   double Rsqij;
   double Rsqjk;
   double Rsqik;
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
   double Rij[DIM];
   double Rjk[DIM];
   double Rik[DIM];
   double *pRij = &(Rij[0]);
   double *pRjk = &(Rjk[0]);
   double *pRik = &(Rik[0]);
   double Rij_pairs[2][3];
   double *pRij_pairs = &(Rij_pairs[0][0]);
   int ier;
   int i;
   int i_pairs[2];
   int *pi_pairs = &(i_pairs[0]);
   int j;
   int j_pairs[2];
   int *pj_pairs = &(j_pairs[0]);
   int k;
   int jj;
   int kk;
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
   int* particleSpecies;
   double* cutoff_ij;
   double* cutsq_ij;
   double* cutoff_ik;
   double* cutsq_ik;
   double* A;
   double* B;
   double* p;
   double* q;
   double* sigma;
   double* lambda_ij;
   double* lambda_ik;
   double lambda;
   double* gamma_ij;
   double* gamma_ik;
   double* costheta;
   double* Rij_list;
   double* coords;
   double* energy;
   double* force;
   double* particleEnergy;
   double* boxSideLengths;
   int numOfAtomNeigh;
   int num_model_species;
   int iSpecies;
   int jSpecies;
   int kSpecies;
   int index_ij;   /* index of parameters for interaction between species i and j */
   int index_ik;
   
   typedef int (*get_neigh_ptr)(void *,int *,int *,int *, int *, int **, double **);
   get_neigh_ptr get_neigh;

   dphi_three = (double *) malloc(3*sizeof(double));
   d2phi_three = (double *) malloc(6*sizeof(double));
   dEidr_three = (double *) malloc(3*sizeof(double));
   d2Eidr_three = (double *) malloc(6*sizeof(double));

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
   num_model_species = buffer->num_model_species;

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

   KIM_API_getm_data_by_index(pkim, &ier, 7*3,
                              buffer->numberOfParticles_ind,           &nAtoms,         1,
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
   ier = KIM_STATUS_FAIL; /* assume an error */
   for (i = 0; i < *nAtoms; ++i)
   {
      if (particleSpecies[i] > num_model_species)
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

        /* get corresponding parameters */
        index_ij    = get_param_index(iSpecies, jSpecies, num_model_species);  
        cutoff_ij   = &buffer->cutoff[index_ij];
        cutsq_ij    = &buffer->cutsq[index_ij];
        A           = &buffer->A[index_ij];
        B           = &buffer->B[index_ij];
        p           = &buffer->p[index_ij];
        q           = &buffer->q[index_ij];
        sigma       = &buffer->sigma[index_ij];
        lambda_ij   = &buffer->lambda[index_ij];
        gamma_ij    = &buffer->gamma[index_ij];
        costheta    = &buffer->costheta[index_ij];

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
         if (Rsqij > *cutsq_ij) continue; /* particles are not interacting  */
         R1 = sqrt(Rsqij);
         if (comp_process_d2Edr2)
         {
             /* compute pair potential and its derivatives */
             calc_phi_d2phi_two(*A, *B, *p, *q, *cutoff_ij, *sigma,
                                R1, &phi_two, &dphi_two, &d2phi_two);

             /* compute dEidr */
             dEidr_two  = 0.5*dphi_two;
             d2Eidr_two = 0.5*d2phi_two;
         }
         else if (comp_force || comp_process_dEdr)
         {
             /* compute pair potential and its derivative */
             calc_phi_dphi_two(*A, *B, *p, *q, *cutoff_ij, *sigma,
                               R1, &phi_two, &dphi_two);

             dEidr_two = 0.5*dphi_two;
         }
         else
         {
            /* compute just pair potential */
             calc_phi_two(*A, *B, *p, *q, *cutoff_ij, *sigma,
                          R1, &phi_two);
         }

         /* contribution to energy */
         if (comp_particleEnergy)
         {
            particleEnergy[i] += 0.5*phi_two;
         }
         if (comp_energy)
         {
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
         for (kk = jj+1; kk < numOfAtomNeigh; ++kk)
         {
           k = neighListOfCurrentAtom[kk] + model_index_shift; /* get neighbor ID */
           kSpecies = particleSpecies[k];

           /* get corresponding parameters */
           index_ik    = get_param_index(iSpecies, kSpecies, num_model_species);  
           cutoff_ik   = &buffer->cutoff[index_ik];
           cutsq_ik    = &buffer->cutsq[index_ik];
           lambda_ik   = &buffer->lambda[index_ik];
           gamma_ik    = &buffer->gamma[index_ik];

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
              }

              /* compute squared distance */
              Rsqik += Rik[kdim]*Rik[kdim];
              Rsqjk += Rjk[kdim]*Rjk[kdim];
           }

           /* compute energy and force */
           if (Rsqik > *cutsq_ik) continue; /* particles are interacting ? */
            R2 = sqrt(Rsqik);
            R3 = sqrt(Rsqjk);

            if (comp_process_d2Edr2)
            {
               /* compute three-body potential and its derivatives */
               lambda = sqrt(fabs(*lambda_ij)*fabs(*lambda_ik));
               calc_phi_d2phi_three(*cutoff_ij, *cutoff_ik, lambda,
                                    *gamma_ij, *gamma_ik, *costheta,
                                    R1, R2, R3, &phi_three, dphi_three, d2phi_three);

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
               lambda = sqrt(fabs(*lambda_ij)*fabs(*lambda_ik));
               calc_phi_dphi_three(*cutoff_ij, *cutoff_ik, lambda,
                                   *gamma_ij, *gamma_ik, *costheta,
                                   R1, R2, R3, &phi_three, dphi_three);

               dEidr_three[0]  =  dphi_three[0];
               dEidr_three[1]  =  dphi_three[1];
               dEidr_three[2]  =  dphi_three[2];
            }
            else
            {
               /* compute just three-body potential */
               lambda = sqrt(fabs(*lambda_ij)*fabs(*lambda_ik));
               calc_phi_three(*cutoff_ij, *cutoff_ik, lambda,
                              *gamma_ij, *gamma_ik, *costheta,
                              R1, R2, R3, &phi_three);
            }

            /* contribution to energy */
            if (comp_particleEnergy)
            {
               particleEnergy[i] += phi_three;
            }
            if (comp_energy)
            {
              *energy += phi_three;
            }

            /* contribution to process_dEdr */
            if (comp_process_dEdr)
            {
               ier = KIM_API_process_dEdr(km, dEidr_three, &R1, &pRij, &i, &j);
               ier = KIM_API_process_dEdr(km, dEidr_three+1, &R2, &pRik, &i, &k);
               ier = KIM_API_process_dEdr(km, dEidr_three+2, &R3, &pRjk, &j, &k);
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

               R_pairs[0] = R_pairs[1] = R3;
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
               R_pairs[1] = R3;
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

               R_pairs[0] = R3;
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
               R_pairs[1] = R3;
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

               R_pairs[0] = R3;
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
                  force[i*DIM + kdim] += dEidr_three[0]*Rij[kdim]/R1 + dEidr_three[1]*Rik[kdim]/R2;
                  force[j*DIM + kdim] -= dEidr_three[0]*Rij[kdim]/R1 - dEidr_three[2]*Rjk[kdim]/R3;
                  force[k*DIM + kdim] -= dEidr_three[2]*Rjk[kdim]/R3 + dEidr_three[1]*Rik[kdim]/R2;
               }
            }
         } /* loop on kk */

      /* End adding three body terms */
      /*******************************/

      } /* loop on jj */
   }    /* infinite while loop terminated by break statements above */

   /* Free temporary storage */
   if (3 == NBC)
   {
      free(neighListOfCurrentAtom);
   }

   free(dphi_three);
   free(d2phi_three);
   free(dEidr_three);
   free(d2Eidr_three);

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
   double* model_cutoff;
   double* model_cutsq;
   double* model_A;
   double* model_B;
   double* model_p;
   double* model_q;
   double* model_lambda;
   double* model_gamma;
   double* model_sigma;
   double* model_costheta;

   int ier;
   struct model_buffer* buffer;
   const char* NBCstr;
   double* cutoff;
   int num_model_species;
   int num_interactions;
   int num_sim_species;
   int* sim_species_code; 
   const char* species;
   int index; 
   int dummy;
   int i, j;

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
      KIM_API_report_error(__LINE__, __FILE__, "Unable to open parameter file", ier);
      return ier;
   }

   /* read number of species */
   ier = fscanf(fid, "%d\n", &num_model_species);
   if (ier != 1)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "error reading first line of parameter file", ier);
      return ier;
   }

   num_interactions = (num_model_species + 1)*num_model_species/2;  


   /* allocate memory for parameters */
   model_A        = (double*) malloc(num_interactions*sizeof(double));  
   model_B        = (double*) malloc(num_interactions*sizeof(double));  
   model_p        = (double*) malloc(num_interactions*sizeof(double));  
   model_q        = (double*) malloc(num_interactions*sizeof(double));  
   model_sigma    = (double*) malloc(num_interactions*sizeof(double));  
   model_lambda   = (double*) malloc(num_interactions*sizeof(double));  
   model_gamma    = (double*) malloc(num_interactions*sizeof(double));  
   model_costheta = (double*) malloc(num_interactions*sizeof(double));  
   model_cutoff   = (double*) malloc(num_interactions*sizeof(double));  
   model_cutsq    = (double*) malloc(num_interactions*sizeof(double));  

   if( model_cutoff==NULL
     || model_cutsq==NULL
     || model_A==NULL
     || model_B==NULL
     || model_p==NULL
     || model_q==NULL
     || model_sigma==NULL
     || model_lambda==NULL
     || model_gamma==NULL
     || model_costheta==NULL )
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
      return ier;
   }

   /* read parameters for two body interactions */
   for (i=0; i< num_interactions; ++i) {
     ier = fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                  &model_A[i],
                  &model_B[i],
                  &model_p[i],
                  &model_q[i],
                  &model_sigma[i],
                  &model_lambda[i],
                  &model_gamma[i],
                  &model_costheta[i],
                  &model_cutoff[i]);
     /* check that we read the right number of parameters */
     if (9 != ier)
     {
       ier = KIM_STATUS_FAIL;
       KIM_API_report_error(__LINE__, __FILE__, "corrupted parameter file", ier);
       return ier;
     }
   }
   /* close param file */
   fclose(fid);

   /* convert units */
   for (i=0; i< num_interactions; ++i) {
     model_A[i] *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
         0.0, 1.0, 0.0, 0.0, 0.0, &ier);
     if (KIM_STATUS_OK > ier)
     {
       KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
       return ier;
     }

     model_sigma[i] *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
         1.0, 0.0, 0.0, 0.0, 0.0, &ier);
     if (KIM_STATUS_OK > ier)
     {
       KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
       return ier;
     }

     model_lambda[i] *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
         0.0, 1.0, 0.0, 0.0, 0.0, &ier);
     if (KIM_STATUS_OK > ier)
     {
       KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
       return ier;
     }

     model_gamma[i] *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
         1.0, 0.0, 0.0, 0.0, 0.0, &ier);
     if (KIM_STATUS_OK > ier)
     {
       KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
       return ier;
     }
   }

   /* store model cutoff in KIM object */
   cutoff = (double*) KIM_API_get_data(pkim, "cutoff", &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", ier);
      return ier;
   }
  
   /* get the number of species in simulator and the corresponding species code */
   ier = KIM_API_get_num_sim_species(pkim, &num_sim_species, &dummy);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_num_sim_species", ier);
      return ier;
   }
   
   sim_species_code = (int*)malloc(num_sim_species*sizeof(int));
   for (i=0; i<num_sim_species; i++)
   {
     ier = KIM_API_get_sim_species(pkim, i, &species);
     if (KIM_STATUS_OK > ier)
     {
       KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_num_sim_species", ier);
       return ier;
     }

     sim_species_code[i] = KIM_API_get_species_code(pkim, species, &ier);
     if (KIM_STATUS_OK > ier)
     {
       KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_num_sim_species", ier);
       return ier;
     }
   }

   for (i=0; i< num_interactions; ++i) {
     model_cutsq[i] = model_cutoff[i]*model_cutoff[i];
   }
   
   /* let KIM cutoff be the largest among the species interaction in simulator */
   *cutoff = 0.0;
   for (i=0; i<num_sim_species; ++i) {
     for (j=i; j<num_sim_species; ++j) {
       index = get_param_index(sim_species_code[i], sim_species_code[j], num_model_species);
       if (model_cutoff[index] > *cutoff) {
         *cutoff = model_cutoff[index];
       }
     }
   }
   
   /* store parameters in KIM object */
   KIM_API_setm_data(pkim, &ier, 10*4,
                             "PARAM_FREE_cutoff",   num_interactions,  model_cutoff,       1,
                             "PARAM_FIXED_cutsq",   num_interactions,  model_cutsq,        1,
                             "PARAM_FREE_A",        num_interactions,  model_A,            1,
                             "PARAM_FREE_B",        num_interactions,  model_B,            1,
                             "PARAM_FREE_p",        num_interactions,  model_p,            1,
                             "PARAM_FREE_q",        num_interactions,  model_q,            1,
                             "PARAM_FREE_sigma",    num_interactions,  model_sigma,        1,
                             "PARAM_FREE_lambda",   num_interactions,  model_lambda,       1,
                             "PARAM_FREE_gamma",    num_interactions,  model_gamma,        1,
                             "PARAM_FREE_costheta", num_interactions,  model_costheta,     1);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_setm_data", ier);
      return ier;
   }

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
                      "particleSpecies",             &(buffer->particleSpecies_ind),             1,
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
    
   /* store simulator species code and number of model species */
   buffer->sim_species_code = sim_species_code;
   buffer->num_model_species = num_model_species;
   
   /* store parameters in buffer */
   buffer->cutoff     = model_cutoff;
   buffer->cutsq      = model_cutsq;
   buffer->A          = model_A;
   buffer->B          = model_B;
   buffer->p          = model_p;
   buffer->q          = model_q;
   buffer->sigma      = model_sigma;
   buffer->lambda     = model_lambda;
   buffer->gamma      = model_gamma;
   buffer->costheta   = model_costheta;

   /* store in model buffer */
   KIM_API_set_model_buffer(pkim, (void*) buffer, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_set_model_buffer", ier);
      return ier;
   }

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
   int* sim_species_code;
   int num_model_species;
   int index;
   int i, j;

   /* get buffer from KIM object */
   buffer = (struct model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
      return ier;
   }
   sim_species_code = buffer->sim_species_code;
   num_model_species = buffer->num_model_species;
   
   /* update cutoff in KIM API and also cutsq */
   cutoff = KIM_API_get_data(pkim, "cutoff", &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", ier);
      return ier;
   }

   nSpecies = KIM_API_get_data_by_index(pkim, buffer->numberOfSpecies_ind, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data_by_index", ier);
      return ier;
   }


   /* let KIM cutoff be the largest among the species interaction in simulator */
   *cutoff = 0.0;
   for (i=0; i < *nSpecies; ++i) {
     for (j=i; j < *nSpecies; ++j) {
       index = get_param_index(sim_species_code[i], sim_species_code[j], num_model_species);
       buffer->cutsq[index] = buffer->cutoff[index]*buffer->cutoff[index];
       if (buffer->cutoff[index] > *cutoff) {
         *cutoff = buffer->cutoff[index];
       }
     }
   }
   
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
    
   /* free  */
   free(buffer->sim_species_code);
   /* free parameters */
   free(buffer->cutoff);
   free(buffer->cutsq);
   free(buffer->A);
   free(buffer->B);
   free(buffer->p);
   free(buffer->q);
   free(buffer->lambda);
   free(buffer->gamma);
   free(buffer->sigma);
   free(buffer->costheta);

   /* destroy the buffer */
   free(buffer);

   ier = KIM_STATUS_OK;
   return ier;
}


/* species number starting from 1 */
static int get_param_index(int i, int j, int num_species)
{
  int index;
  int k;

  /* swap i j if  i>j */
  if (i > j) {
    k = i;
    i = j;
    j = k;
  }

  index = (i-1)*(num_species-1) - (i-1)*(i-2)/2 + (j-1);

  return index;
}


