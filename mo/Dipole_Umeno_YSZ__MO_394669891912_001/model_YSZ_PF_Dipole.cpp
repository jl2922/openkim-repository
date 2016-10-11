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

* Copyright (c) 2012, Institute of Industrial Science, the University of Tokyo.
* All rights reserved.
*
* Contributors: Yoshitaka Umeno
* First written Oct 2012
* Modified Jan 2013:
*  - Reallocation of buffer arrays when numberOfParticles changes
*  - Missing check of comp_energy and comp_force is added
*  - model_destroy funtion is added
* Modified Oct 2013:
*  - Bug fix: sign in E_ind calculation
*/

/*******************************************************************************
*
*  model_YSZ_PF_Dipole
*
*  Dipole model potential for YSZ (Yttria Stabilized Zirconia)
*
*  Reference: Parameters for YSZ are not yet published (and are likely to be modified).
*             For the dipole model, see:
*              Tangney and Scandolo, J.Chem.Phys.117,8898-8904(2002) and
*              Beck et al., J.Chem.Phys.135,234512(2011) (9pages)
*
*  Language: CPP
*
*  Release:
*
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
#define MODEL_CUTOFF 12.0 /* cutoff radius in angstroms */
#define MODEL_CUTSQ  (MODEL_CUTOFF * MODEL_CUTOFF)
#define DP_EPS 14.40
#define DP_CUT 12.0
#define DP_TOL 1.e-7
#define DP_MIX 0.2
//#define ntype 3
//#define nptype 3

/* Define prototypes for model init */
/* must be all lowercase to be compatible with the KIM API (to support Fortran Tests) */
/**/
extern "C" void model_init(void* km);

/* Define prototypes for model reinit, compute, and destroy */
/* defined as static to avoid namespace clashes with other Models */
/**/
//static void compute(void* km, int* ier);
static void compute(void* km);
static int model_destroy(void* km);
/**/
static void buckingham(double r, double &e, double &f, double a, double sig, double c);
static void elstat_shift(double r, double dp_kappa, double &value_tail, double &grad_tail,
		  double &ggrad_tail);
static void elstat_value(double r, double dp_kappa, double &ftail, double &gtail,
		  double &ggtail);
static double shortrange_value(double r, double a, double b, double c);
static void shortrange_term(double r, double b, double c, double &srval_tail,
		     double &srgrad_tail);
static int pairtyp(int i, int j);
static double dsquare(double d);

// Define model_buffer structure
struct model_buffer {
  double *dp_alpha, *dp_b, *dp_c;
  double *E_statx, *E_staty, *E_statz;
  double *E_indx, *E_indy, *E_indz;
  double *E_oldx, *E_oldy, *E_oldz;
  double *E_totx, *E_toty, *E_totz;
  double *p_srx, *p_sry, *p_srz;
  double *p_indx, *p_indy, *p_indz;
  //double dp_eps, dp_cut, dp_tol, dp_mix;
  double *ratio, *charge, last_charge, dp_kappa;
  double *buck_a, *buck_s, *buck_c;
  int sw_kappa;
  bool initialize;
  int ntype, nptype;
  int nparticle;
};

static void buckingham(double r, double &e, double &f, double a, double sig, double c)
{
  //  double rang = r;
  double x = sig/r; double x2 = x * x; double x3 = x2 * x;
  double x6 = x3 * x3; double x5 = x3 * x2;
  double exprs = exp(-r/sig);
  e = a * exprs - c * x6;
  f = -1.0/sig * a * exprs + c * 6.0 * x6 / r;
}

static double shortrange_value(double r, double a, double b, double c)
{
  static double x[5];

  x[0] = b * r;
  x[1] = x[0] * x[0];
  x[2] = x[1] * x[0];
  x[3] = x[1] * x[1];
  x[4] = 1 + x[0] + x[1] / 2 + x[2] / 6 + x[3] / 24;

  return a * c * x[4] * exp(-x[0]) / DP_EPS;
}
void shortrange_term(double r, double b, double c, double &srval_tail,
  double &srgrad_tail)
{
  static double x[6];

  x[0] = b * r;
  x[1] = x[0] * x[0];
  x[2] = x[1] * x[0];
  x[3] = x[1] * x[1];
  x[4] = 1 + x[0] + x[1] / 2 + x[2] / 6 + x[3] / 24;
  x[5] = exp(-x[0]);

  srval_tail = c * x[4] * x[5] / DP_EPS;
  //srgrad_tail = -c * b * x[3] * x[5] / (24 * DP_EPS * r);
  srgrad_tail = -c * b * x[3] * x[5] / (24 * DP_EPS);
}

static void elstat_shift(double r, double dp_kappa, double &value_tail, double &grad_tail,
  double &ggrad_tail)
{
  static double ftail, gtail, ggtail, ftail_cut, gtail_cut, ggtail_cut;
  static double x[3];

  x[0] = r * r;
  x[1] = DP_CUT * DP_CUT;
  x[2] = x[0] - x[1];

  elstat_value(r, dp_kappa, ftail, gtail, ggtail);
  elstat_value(DP_CUT, dp_kappa, ftail_cut, gtail_cut, ggtail_cut);

  value_tail = ftail - ftail_cut - x[2] * gtail_cut / 2;
  grad_tail = gtail - gtail_cut;
  ggrad_tail = 0.;

  value_tail -= x[2] * x[2] * ggtail_cut / 8;
  grad_tail -= x[2] * ggtail_cut / 2;
  ggrad_tail = ggtail - ggtail_cut;
}

static void elstat_value(double r, double dp_kappa, double &ftail, double &gtail, double &ggtail)
{
  static double x[4];

  x[0] = r * r;
  x[1] = dp_kappa * dp_kappa;
  x[2] = 2 * DP_EPS * dp_kappa / sqrt(M_PI);
  x[3] = exp(-x[0] * x[1]);

  ftail = DP_EPS * erfc(dp_kappa * r) / r;
  gtail = -(ftail + x[2] * x[3]) / x[0];
  ggtail = (2 * x[1] * x[2] * x[3] - gtail * 3) / x[0];
}

static double dsquare(double d)
{
  return d * d;
}

// pairtyp is to give index for pair type
// 1:Zr-O, 2:O-O 3:Y-O (Caution: 1-based numbering here) (OLD)
// 1:Zr-Zr, 2:Zr-O 3:Zr-Y, 4:O-O, 5:O-Y, 6:Y-Y (NEW)
static int pairtyp(int i, int j)
{
  if ((i==1)&&(j==1)) { return 1;
  } else if (((i==1)&&(j==2))||((i==2)&&(j==1))) { return 2;
  } else if (((i==1)&&(j==3))||((i==3)&&(j==1))) { return 3;
  } else if ((i==2)&&(j==2)) { return 4;
  } else if (((i==2)&&(j==3))||((i==3)&&(j==2))) { return 5;
  } else if ((i==3)&&(j==3)) { return 6;
  } else { return 0; }
  /*
  if (((i==1)&&(j==2))||((i==2)&&(j==1))) { return 1;
  } else if ((i==2)&&(j==2)) { return 2;
  } else if (((i==2)&&(j==3))||((i==3)&&(j==2))) { return 3;
  } else { return 0; }
  */
}

// destroy function
static int model_destroy(void* km)
{
  intptr_t* pkim = *((intptr_t**) km);
  struct model_buffer* buffer;
  int ier;
  buffer = (struct model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
  delete[] buffer->E_statx, buffer->E_staty, buffer->E_statz;
  delete[] buffer->E_indx,  buffer->E_indy,  buffer->E_indz;
  delete[] buffer->E_oldx,  buffer->E_oldy,  buffer->E_oldz;
  delete[] buffer->E_totx,  buffer->E_toty,  buffer->E_totz;
  delete[] buffer->p_srx,  buffer->p_sry,  buffer->p_srz;
  delete[] buffer->p_indx, buffer->p_indy, buffer->p_indz;
  delete buffer->dp_alpha, buffer->dp_b, buffer->dp_c;
  delete buffer->charge, buffer->buck_a, buffer->buck_s, buffer->buck_c;
  free(buffer);
  ier = KIM_STATUS_OK;
  return ier;
}

/* compute function */
//static void compute(void* km, int* ier)
static void compute(void* km)
{
   /* local static parameters */
   const static double cutsq = MODEL_CUTOFF * MODEL_CUTOFF;
   /* local variables */
   intptr_t* pkim = *((intptr_t**) km);
   double R;
   double Rsqij;
   double phi;
   double dphi;
   double dEidr;
   double Rij[DIM];
   int i;
   int j;
   int jj;
   int k;
   int numOfAtomNeigh;
   int currentAtom;
   int comp_energy;
   int comp_force;
   int comp_particleEnergy;
   int comp_virial;
   int IterOrLoca;
   int HalfOrFull;
   int NBC;
   const char* NBCstr;
   int numberContrib;

   int* nAtoms;
   int* particleSpecies;
   double* Rij_list;
   double* coords;
   double* energy;
   double* force;
   double* particleEnergy;
   double* virial;
   int* neighListOfCurrentAtom;
   double* boxSideLengths;
   int* numContrib;

   // For Dipole
   struct model_buffer* buffer;
   double rr, rr2, drx, dry, drz, vp0, v0;
   int ix, iy, iz;
   double rcut = DP_CUT; double rcut2 = rcut * rcut;
   int    self;
   int    h, l, typ1, typ2, uf, us, stresses, ipair;
   double value, grad, value_tail, grad_tail, grad_i, grad_j, p_sr_tail;
   double value_el, grad_el, ggrad_el; 
   // buffered...
   double *dp_alpha, *dp_b, *dp_c;
   double *E_statx, *E_staty, *E_statz;
   double *E_indx, *E_indy, *E_indz;
   double *E_oldx, *E_oldy, *E_oldz;
   double *E_totx, *E_toty, *E_totz;
   double *p_srx, *p_sry, *p_srz;
   double *p_indx, *p_indy, *p_indz;
   double *ratio, *charge, last_charge, dp_kappa;
   double *buck_a, *buck_s, *buck_c;
   int sw_kappa;
   bool initialize;
   int ntype, nptype;
   int nparticle;
   int status;

   //printf("### Dipole ###\n");

   // get buffer from KIM object
   //buffer = (struct model_buffer*) KIM_API_get_model_buffer(pkim, ier);
   buffer = (struct model_buffer*) KIM_API_get_model_buffer(pkim, &status);

   //KIM_API_print(pkim,&status);
   //printf("unpack\n");
   if (KIM_STATUS_OK > status)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", status);
      return;
   }
   int *natom;
   natom = (int*) KIM_API_get_data(pkim, "numberOfParticles", &status);
   //printf("check %d %d\n",buffer->nparticle,*natom);
   if (buffer->nparticle != *natom) {
     printf("-- reallocate buffer arrays --\n");
     buffer->nparticle = *natom;
     delete[] buffer->E_statx, buffer->E_staty, buffer->E_statz;
     delete[] buffer->E_indx,  buffer->E_indy,  buffer->E_indz;
     delete[] buffer->E_oldx,  buffer->E_oldy,  buffer->E_oldz;
     delete[] buffer->E_totx,  buffer->E_toty,  buffer->E_totz;
     delete[] buffer->p_srx,  buffer->p_sry,  buffer->p_srz;
     delete[] buffer->p_indx, buffer->p_indy, buffer->p_indz;

     E_statx  = new double[*natom];
     E_staty  = new double[*natom];
     E_statz  = new double[*natom];
     E_indx   = new double[*natom];
     E_indy   = new double[*natom];
     E_indz   = new double[*natom];
     E_oldx   = new double[*natom];
     E_oldy   = new double[*natom];
     E_oldz   = new double[*natom];
     E_totx   = new double[*natom];
     E_toty   = new double[*natom];
     E_totz   = new double[*natom];
     p_srx    = new double[*natom];
     p_sry    = new double[*natom];
     p_srz    = new double[*natom];
     p_indx   = new double[*natom];
     p_indy   = new double[*natom];
     p_indz   = new double[*natom];
     buffer->E_statx = E_statx;
     buffer->E_staty = E_staty;
     buffer->E_statz = E_statz;
     buffer->E_indx = E_indx;
     buffer->E_indy = E_indy;
     buffer->E_indz = E_indz;
     buffer->E_oldx = E_oldx;
     buffer->E_oldy = E_oldy;
     buffer->E_oldz = E_oldz;
     buffer->E_totx = E_totx;
     buffer->E_toty = E_toty;
     buffer->E_totz = E_totz;
     buffer->p_srx = p_srx;
     buffer->p_sry = p_sry;
     buffer->p_srz = p_srz;
     buffer->p_indx = p_indx;
     buffer->p_indy = p_indy;
     buffer->p_indz = p_indz;
   }

   // unpack info from the buffer
   dp_alpha = buffer->dp_alpha;
   dp_b = buffer->dp_b;
   dp_c = buffer->dp_c;
   E_statx = buffer->E_statx;
   E_staty = buffer->E_staty;
   E_statz = buffer->E_statz;
   E_indx = buffer->E_indx;
   E_indy = buffer->E_indy;
   E_indz = buffer->E_indz;
   E_oldx = buffer->E_oldx;
   E_oldy = buffer->E_oldy;
   E_oldz = buffer->E_oldz;
   E_totx = buffer->E_totx;
   E_toty = buffer->E_toty;
   E_totz = buffer->E_totz;
   p_srx = buffer->p_srx;
   p_sry = buffer->p_sry;
   p_srz = buffer->p_srz;
   p_indx = buffer->p_indx;
   p_indy = buffer->p_indy;
   p_indz = buffer->p_indz;
   ratio = buffer->ratio;
   charge = buffer->charge;
   last_charge = buffer->last_charge;
   dp_kappa = buffer->dp_kappa;
   buck_a = buffer->buck_a;
   buck_s = buffer->buck_s;
   buck_c = buffer->buck_c;
   sw_kappa = buffer->sw_kappa;
   initialize = buffer->initialize;
   ntype = buffer->ntype;
   nptype = buffer->nptype;
   nparticle = buffer->nparticle;
   // unpack end  
   //printf("Unpack end\n");

   /* Determine neighbor list boundary condition (NBC) */
   /* and half versus full mode: */
   /*****************************
    * HalfOrFull = 1 -- Half
    *            = 2 -- Full
    *****************************/
   //NBCstr = (char*) KIM_API_get_NBC_method(pkim, ier);
   status = KIM_API_get_NBC_method(pkim, &NBCstr);
   //if (KIM_STATUS_OK > *ier)
   if (KIM_STATUS_OK > status)
   {
     //KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_NBC_method", *ier);
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_NBC_method", status);
      return;
   }
   /*
   if (!strcmp("CLUSTER",NBCstr))
   {
      NBC = 0;
      HalfOrFull = 1;
   }
   else if (!strcmp("MI_OPBC_H",NBCstr))
   {
      NBC = 1;
      HalfOrFull = 1;
   }
   */
   if (!strcmp("MI_OPBC_F",NBCstr))
   {
      NBC = 1;
      HalfOrFull = 2;
   }

   //else if (!strcmp("NEIGH_PURE_H",NBCstr))
   //{
   //   NBC = 2;
   //   HalfOrFull = 1;
   //}
   else if (!strcmp("NEIGH_PURE_F",NBCstr))
   {
      NBC = 2;
      HalfOrFull = 2;
   }
   else if (!strcmp("NEIGH_RVEC_F",NBCstr))
   {
      NBC = 3;
      HalfOrFull = 2;
   }
   else
   {
      status = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "Unknown NBC method", status);
      return;
   }

   /* determine neighbor list handling mode */
   if (NBC != 0)
   {
      /*****************************
       * IterOrLoca = 1 -- Iterator
       *            = 2 -- Locator
       *****************************/
      IterOrLoca = KIM_API_get_neigh_mode(pkim, &status);
      if (KIM_STATUS_OK > status)
      {
         KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh_mode", status);
         return;
      }
      //if ((IterOrLoca != 1) && (IterOrLoca != 2))
	//if (IterOrLoca != 2)
      //{
      //   printf("* ERROR: Unsupported IterOrLoca mode = %i\n", IterOrLoca);
      //   exit(-1);
      //}
   }
   else
   {
     //IterOrLoca = 2;   /* for CLUSTER NBC */
      printf("* ERROR: Cluster is unsupported\n");
      exit(-1);
   }

   /* check to see if we have been asked to compute the forces, particleEnergy, energy and virial */
   KIM_API_getm_compute(pkim, &status, 3*3,
                        "energy",         &comp_energy,         1,
                        "forces",         &comp_force,          1,
                        "particleEnergy", &comp_particleEnergy, 1);
   //                        "virial",         &comp_virial,         1);
   if (KIM_STATUS_OK > status)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute", status);
      return;
   }

   /* unpack data from KIM object */
   KIM_API_getm_data(pkim, &status, 8*3,
                     "numberOfParticles",           &nAtoms,         1,
                     "particleSpecies",               &particleSpecies,  1,
                     "coordinates",                 &coords,         1,
                     "numberContributingParticles", &numContrib,     (HalfOrFull==1),
                     "boxSideLengths",              &boxSideLengths, (NBC==1),
                     "energy",                      &energy,         (comp_energy==1),
                     "forces",                      &force,          (comp_force==1),
                     "particleEnergy",              &particleEnergy, (comp_particleEnergy==1));
   //                     "virial",                      &virial,         (comp_virial==1));
   if (KIM_STATUS_OK > status)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", status);
      return;
   }
   /*
   if (HalfOrFull == 1)
   {
      if (0 != NBC) // non-CLUSTER cases
      {
         numberContrib = *numContrib;
      }
      else
      {
         numberContrib = *nAtoms;
      }
   }
   */

   /* Check to be sure that the atom types are correct */
   /**/
   status = KIM_STATUS_FAIL; /* assume an error */
   for (i = 0; i < *nAtoms; ++i)
   {
     if ((particleSpecies[i] < 1)||(particleSpecies[i] > 3))
      {
         KIM_API_report_error(__LINE__, __FILE__, "Unexpected species type detected", status);
         return;
      }
   }
   status = KIM_STATUS_OK; /* everything is ok */

   for (i=0; i<*nAtoms; i++) {
     E_totx[i]  = 0; E_toty[i]  = 0; E_totz[i]  = 0;
     E_indx[i]  = 0; E_indy[i]  = 0; E_indz[i]  = 0;
     p_indx[i]  = 0; p_indy[i]  = 0; p_indz[i]  = 0;
     E_statx[i] = 0; E_staty[i] = 0; E_statz[i] = 0;
     p_srx[i]   = 0; p_sry[i]   = 0; p_srz[i]   = 0;
   }
   
   // FIRST LOOP
   /* initialize potential energies, forces, and virial term */
   //printf("First loop\n");
   if (comp_particleEnergy)
   {
      for (i = 0; i < *nAtoms; ++i)
      {
         particleEnergy[i] = 0.0;
      }
   }
   if (comp_energy)
   {
      *energy = 0.0;
   }

   if (comp_force)
   {
      for (i = 0; i < *nAtoms; ++i)
      {
         for (k = 0; k < DIM; ++k)
         {
            force[i*DIM + k] = 0.0;
         }
      }
   }

   if (comp_virial)
   {
      for (i = 0; i < 6; ++i)
      {
         virial[i] = 0.0;
      }
   }

   /* Initialize neighbor handling for CLUSTER NBC */
   /*
   if (0 == NBC) // CLUSTER
   {
      neighListOfCurrentAtom = (int *) malloc((*nAtoms)*sizeof(int));
   }
   */

   // SECOND LOOP: calculate short-range and monopole forces,
   // calculate static field- and dipole-contributions
   //printf("Second loop\n");
   for (i = 0; i < *nAtoms; i++)
   {
     status = KIM_API_get_neigh(pkim, 1, i, &currentAtom, &numOfAtomNeigh,
			      &neighListOfCurrentAtom, &Rij_list);
     if (KIM_STATUS_OK != status) /* some sort of problem, exit */
       {
	 KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", status);
	 status = KIM_STATUS_FAIL;
	 return;
       }
     typ1 = particleSpecies[i];
     /* loop over the neighbors of atom i */
     for (jj = 0; jj < numOfAtomNeigh; ++ jj)
       {
         j = neighListOfCurrentAtom[jj]; /* get neighbor ID */
	 if (j<i) continue;
         /* compute relative position vector and squared distance */
         Rsqij = 0.0;
         for (k = 0; k < DIM; ++k)
	   {
	     if (3 != NBC) /* all methods except NEIGH_RVEC_F */
	       {
		 Rij[k] = coords[j*DIM + k] - coords[i*DIM + k];
	       }
	     else          /* NEIGH_RVEC_F method */
	       {
		 Rij[k] = Rij_list[jj*DIM + k];
	       }
	     
	     /* apply periodic boundary conditions if required */
	     if (1 == NBC)
	       {
		 if (abs(Rij[k]) > 0.5*boxSideLengths[k])
		   {
		     Rij[k] -= (Rij[k]/fabs(Rij[k]))*boxSideLengths[k];
		   }
	       }
	     
	     /* compute squared distance */
	     Rsqij += Rij[k]*Rij[k];
	   }

         // contribution to electron density
         if (Rsqij < MODEL_CUTSQ) /* particles are interacting ? */
         {
            R = sqrt(Rsqij);
	    rr2 = Rsqij; rr = R;
	    typ2 = particleSpecies[j];
	    ipair = pairtyp(typ1, typ2);
	    elstat_shift(rr, dp_kappa, value_el, grad_el, ggrad_el); //tail-functions
	    drx = Rij[0]; dry = Rij[1]; drz = Rij[2];
	    if ((rr < rcut)&&(ipair!=0)) { // calculate short-range forces
	      self = 0; if (i==j) { self = 1; }
	      buckingham(rr, value, grad, buck_a[ipair-1], buck_s[ipair-1]
			 , buck_c[ipair-1]); //Buckingham type
	      if (self) { value *= 0.5; grad *= 0.5; }
	      grad /= rr;
	      if (comp_particleEnergy) {
		particleEnergy[i] += value/2.0; particleEnergy[j] += value/2.0; }
	      if (comp_energy) {
		*energy +=value; }
	      if (comp_force) {
	      force[i*DIM] += drx*grad; force[i*DIM+1] += dry*grad; force[i*DIM+2] += drz*grad;
	      force[j*DIM] -= drx*grad; force[j*DIM+1] -= dry*grad; force[j*DIM+2] -= drz*grad;
	      }
	    }
	    if (rr < DP_CUT) { // calculate monopole forces
	      self = 0; if (i==j) { self = 1; }
	      value_tail = value_el;
	      grad_tail  = grad_el;
	      grad_i = charge[typ2-1] * grad_tail;
	      if (typ1 == typ2) { grad_j = grad_i;
	      } else { grad_j = charge[typ1-1] * grad_tail; }
	      value = charge[typ1-1] * charge[typ2-1] * value_tail;
	      grad  = charge[typ1-1] * grad_i;
	      if (self) { grad_i *= 0.5; grad_j *= 0.5; value *= 0.5; grad *= 0.5; }
	      if (comp_particleEnergy) {
		particleEnergy[i] += value/2.0; particleEnergy[j] += value/2.0; }
	      if (comp_energy) {
		*energy +=value; }
	      if (comp_force) {
	      force[i*DIM] += drx*grad; force[i*DIM+1] += dry*grad; force[i*DIM+2] += drz*grad;
	      force[j*DIM] -= drx*grad; force[j*DIM+1] -= dry*grad; force[j*DIM+2] -= drz*grad;
	      }
	      // calculate static field-contributions
	      E_statx[i] += drx*grad_i; E_statx[j] -= drx*grad_j;
	      E_staty[i] += dry*grad_i; E_staty[j] -= dry*grad_j;
	      E_statz[i] += drz*grad_i; E_statz[j] -= drz*grad_j;
	      // calculate short-range dipoles
	      if (ipair!=0) {
		p_sr_tail = grad_tail*rr
		  *shortrange_value(rr, dp_alpha[typ1-1], dp_b[ipair-1],
				    dp_c[ipair-1]);
		p_srx[i] += charge[typ2-1]*drx/rr*p_sr_tail;
		p_sry[i] += charge[typ2-1]*dry/rr*p_sr_tail;
		p_srz[i] += charge[typ2-1]*drz/rr*p_sr_tail;
		if (!self) {
		  p_sr_tail = grad_tail*rr
		    *shortrange_value(rr, dp_alpha[typ2-1], dp_b[ipair-1],
				      dp_c[ipair-1]);
		  p_srx[j] -= charge[typ1-1]*drx/rr*p_sr_tail;
		  p_sry[j] -= charge[typ1-1]*dry/rr*p_sr_tail;
		  p_srz[j] -= charge[typ1-1]*drz/rr*p_sr_tail; }
	      }
	    }
         }
       } // jj
   } // i
   // SECOND LOOP END

   // THIRD LOOP: calculate whole dipole moment for every atom
   //printf("Third loop\n");
   double rp, dp_sum;
   int dp_converged = 0, dp_it = 0;
   double max_diff = 10;
   while (dp_converged ==0) {
     dp_sum = 0;
     for (i = 0; i < *nAtoms; i++)
       {
	 status = KIM_API_get_neigh(pkim, 1, i, &currentAtom, &numOfAtomNeigh,
				  &neighListOfCurrentAtom, &Rij_list);
	 typ1 = particleSpecies[i];
	 if (dp_alpha[typ1-1]!=0) {
	   if (dp_it) {
	     E_totx[i] = (1-DP_MIX)*E_indx[i]
	       + DP_MIX*E_oldx[i] + E_statx[i];
	     E_toty[i] = (1-DP_MIX)*E_indy[i]
	       + DP_MIX*E_oldy[i] + E_staty[i];
	     E_totz[i] = (1-DP_MIX)*E_indz[i]
	       + DP_MIX*E_oldz[i] + E_statz[i];
	   } else {
	     E_totx[i] = E_indx[i] + E_statx[i];
	     E_toty[i] = E_indy[i] + E_staty[i];
	     E_totz[i] = E_indz[i] + E_statz[i];
	   } // if dp_it
	   p_indx[i] = dp_alpha[typ1-1] * E_totx[i] + p_srx[i];
	   p_indy[i] = dp_alpha[typ1-1] * E_toty[i] + p_sry[i];
	   p_indz[i] = dp_alpha[typ1-1] * E_totz[i] + p_srz[i];
	   E_oldx[i] = E_indx[i]; E_indx[i] = 0.;
	   E_oldy[i] = E_indy[i]; E_indy[i] = 0.;
	   E_oldz[i] = E_indz[i]; E_indz[i] = 0.;
	 }
       } // loop of i (atoms)
     for (i = 0; i < *nAtoms; i++)
       {
	 status = KIM_API_get_neigh(pkim, 1, i, &currentAtom, &numOfAtomNeigh,
				  &neighListOfCurrentAtom, &Rij_list);
	 typ1 = particleSpecies[i];
	 for (jj = 0; jj < numOfAtomNeigh; ++ jj)
	   {
	     j = neighListOfCurrentAtom[jj];
	     if (j<i) continue;
	     Rsqij = 0.0;
	     for (k = 0; k < DIM; ++k)
	       {
		 if (3 != NBC) /* all methods except NEIGH_RVEC_F */
		   { Rij[k] = coords[j*DIM + k] - coords[i*DIM + k]; }
		 else          /* NEIGH_RVEC_F method */
		   { Rij[k] = Rij_list[jj*DIM + k]; }
		 /* apply periodic boundary conditions if required */
		 if (1 == NBC)
		   {
		     if (abs(Rij[k]) > 0.5*boxSideLengths[k])
		       { Rij[k] -= (Rij[k]/fabs(Rij[k]))*boxSideLengths[k]; }
		   }
		 /* compute squared distance */
		 Rsqij += Rij[k]*Rij[k];
	       }
	     // contribution to electron density
	     if (Rsqij < MODEL_CUTSQ) /* particles are interacting ? */
	       {
		 R = sqrt(Rsqij);
		 rr2 = Rsqij; rr = R;
		 typ2 = particleSpecies[j];
		 ipair = pairtyp(typ1, typ2);
		 elstat_shift(rr, dp_kappa, value_el, grad_el, ggrad_el); //tail-functions
		 if (rr < rcut) { // calculate short-range forces
		   if ((dp_alpha[typ1-1]!=0)&&(dp_alpha[typ2-1]!=0)) {
		     self = 0; if (i==j) { self = 1; }
		     drx = Rij[0]; dry = Rij[1]; drz = Rij[2];
		     drx /=rr; dry /=rr; drz /=rr;
		     rp = p_indx[j]*drx + p_indy[j]*dry + p_indz[j]*drz;
		     E_indx[i] -= grad_el*(3*rp*drx-p_indx[j]);
		     E_indy[i] -= grad_el*(3*rp*dry-p_indy[j]);
		     E_indz[i] -= grad_el*(3*rp*drz-p_indz[j]);
		     if (!self) {
		       rp = p_indx[i]*drx + p_indy[i]*dry + p_indz[i]*drz;
		       E_indx[j] -= grad_el*(3*rp*drx-p_indx[i]);
		       E_indy[j] -= grad_el*(3*rp*dry-p_indy[i]);
		       E_indz[j] -= grad_el*(3*rp*drz-p_indz[i]);
		     }
		   }}
	       } } // j
       } // loop of i (atoms)  
     for (i = 0; i < *nAtoms; i++)
       {
	 typ1 = particleSpecies[i];
	 dp_sum += dsquare(dp_alpha[typ1-1]*(E_oldx[i]-E_indx[i]));
	 dp_sum += dsquare(dp_alpha[typ1-1]*(E_oldy[i]-E_indy[i]));
	 dp_sum += dsquare(dp_alpha[typ1-1]*(E_oldz[i]-E_indz[i]));
       } // loop of i (atoms)
     dp_sum /= 3*((double)(*nAtoms));
     dp_sum = sqrt(dp_sum);
     if (dp_it) {
       if ((dp_sum > max_diff) || (dp_it > 50)) {
	 dp_converged = 1;
	 //printf("dp (failure) dp_sum=%e  dp_it=%d\n",dp_sum,dp_it);
	 for (i = 0; i < *nAtoms; i++) {
	   typ1 = particleSpecies[i];
	   if (dp_alpha[typ1-1]!=0) {
	     p_indx[i] = dp_alpha[typ1-1]*E_statx[i]+p_srx[i];
	     p_indy[i] = dp_alpha[typ1-1]*E_staty[i]+p_sry[i];
	     p_indz[i] = dp_alpha[typ1-1]*E_statz[i]+p_srz[i];
	     E_indx[i] = E_statx[i];
	     E_indy[i] = E_staty[i];
	     E_indz[i] = E_statz[i]; }
	 }
       }
     } // if dp_it
     if (dp_sum < DP_TOL) {
       dp_converged = 1;
       //printf("dp (success) dp_it=%d\n",dp_it);
     }
     dp_it++;
   } // while dp_converged (THIRD LOOP END)
   //printf("dp iter %d\n",dp_it);
   //for (i = 0; i < *nAtoms; i++) {
   //  printf("%d %e\n",i,p_indx[i]); }
   
   // FOURTH LOOP: calculate monopole-dipole and dipole-dipole forces
   //printf("Fourth loop\n");
   double  rp_i, rp_j, pp_ij, tmp_1, tmp_2, tmpx, tmpy, tmpz;
   double  grad_1, grad_2, srval, srgrad, srval_tail, srgrad_tail, value_sum,
     grad_sum;
   for (i = 0; i < *nAtoms; i++)
     {
       status = KIM_API_get_neigh(pkim, 1, i, &currentAtom, &numOfAtomNeigh,
				&neighListOfCurrentAtom, &Rij_list);
       typ1 = particleSpecies[i];
       for (jj = 0; jj < numOfAtomNeigh; ++ jj)
	 {
	   j = neighListOfCurrentAtom[jj];
	   if (j<i) continue;
	   Rsqij = 0.0;
	   for (k = 0; k < DIM; ++k)
	     {
	       if (3 != NBC) /* all methods except NEIGH_RVEC_F */
		 { Rij[k] = coords[j*DIM + k] - coords[i*DIM + k]; }
	       else          /* NEIGH_RVEC_F method */
		 { Rij[k] = Rij_list[jj*DIM + k]; }
	       /* apply periodic boundary conditions if required */
	       if (1 == NBC)
		 {
		   if (abs(Rij[k]) > 0.5*boxSideLengths[k])
		     { Rij[k] -= (Rij[k]/fabs(Rij[k]))*boxSideLengths[k]; }
		 }
	       /* compute squared distance */
	       Rsqij += Rij[k]*Rij[k];
	     }
	   // contribution to electron density
	   if (Rsqij < MODEL_CUTSQ) /* particles are interacting ? */
	     {
	       R = sqrt(Rsqij);
	       rr2 = Rsqij; rr = R;
	       typ2 = particleSpecies[j];
	       ipair = pairtyp(typ1, typ2);
	       elstat_shift(rr, dp_kappa, value_el, grad_el, ggrad_el); //tail-functions
	       if ((rr < DP_CUT)&&((dp_alpha[typ1-1]!=0)||(dp_alpha[typ2-1]!=0))) {
		 self = 0; if (i==j) { self = 1; }
		 value_tail = -grad_el;
		 grad_tail  = -ggrad_el;
		 if (ipair!=0) {
		   shortrange_term(rr, dp_b[ipair-1], dp_c[ipair-1], srval_tail, srgrad_tail);
		   srval = value_tail * srval_tail;
		   srgrad = value_tail * srgrad_tail + grad_tail * srval_tail;
		 }
		 if (self) { value_tail *= 0.5; grad_tail *= 0.5; }
		 // monopole-dipole contributions //if (charge[typ1-1] && dp_alpha[typ2-1]) 
		 if (ipair!=0) {
		   value_sum = value_tail + srval;
		   grad_sum  = grad_tail + srgrad;
		 } else {
		   value_sum = value_tail;
		   grad_sum  = grad_tail;
		 }
		 drx = Rij[0]; dry = Rij[1]; drz = Rij[2];
		 drx /=rr; dry /=rr; drz /=rr;
		 rp_j = p_indx[j]*drx + p_indy[j]*dry + p_indz[j]*drz;
		 value  = charge[typ1-1] * rp_j * value_sum * rr;
		 grad_1 = charge[typ1-1] * rp_j * grad_sum * rr2;
		 grad_2 = charge[typ1-1] * value_sum;
		 if (comp_particleEnergy) {
		   particleEnergy[i] -= value/2.0; particleEnergy[j] -= value/2.0; }
		 if (comp_energy) {
		   *energy -=value; }
		 tmpx = drx * grad_1 + p_indx[j] * grad_2;
		 tmpy = dry * grad_1 + p_indy[j] * grad_2;
		 tmpz = drz * grad_1 + p_indz[j] * grad_2;
		 if (comp_force) {
		 force[i*DIM  ] -= tmpx; force[i*DIM+1] -= tmpy; force[i*DIM+2] -= tmpz;
		 force[j*DIM  ] += tmpx; force[j*DIM+1] += tmpy; force[j*DIM+2] += tmpz;
		 }
		 // dipole-monopole contributions //if (dp_alpha[typ2-1] && charge[typ2-1])
		 if (ipair!=0) {
		   value_sum = value_tail + srval;
		   grad_sum  = grad_tail + srgrad;
		 } else {
		   value_sum = value_tail;
		   grad_sum  = grad_tail;
		 }
		 rp_i = p_indx[i]*drx + p_indy[i]*dry + p_indz[i]*drz;
		 value  = charge[typ2-1] * rp_i * value_sum * rr;
		 grad_1 = charge[typ2-1] * rp_i * grad_sum * rr2;
		 grad_2 = charge[typ2-1] * value_sum;
		 if (comp_particleEnergy) {
		   particleEnergy[i] += value/2.0; particleEnergy[j] += value/2.0; }
		 if (comp_energy) {
		   *energy +=value; }
		 tmpx = drx * grad_1 + p_indx[i] * grad_2;
		 tmpy = dry * grad_1 + p_indy[i] * grad_2;
		 tmpz = drz * grad_1 + p_indz[i] * grad_2;
		 if (comp_force) {
		 force[i*DIM  ] += tmpx; force[i*DIM+1] += tmpy; force[i*DIM+2] += tmpz;
		 force[j*DIM  ] -= tmpx; force[j*DIM+1] -= tmpy; force[j*DIM+2] -= tmpz;
		 }
		 // dipole-dipole contributions //if (dp_alpha[typ1-1] && dp_alpha[typ2-1])
		 if ((dp_alpha[typ1-1]!=0)&&(dp_alpha[typ2-1]!=0)) {
		   pp_ij = p_indx[i]*p_indx[j]
		     + p_indy[i]*p_indy[j] + p_indz[i]*p_indz[j];
		   tmp_1 = 3 * rp_i * rp_j;
		   tmp_2 = 3 * value_tail / rr2;
		   value = -(tmp_1 - pp_ij) * value_tail;
		   grad_1 = (tmp_1 - pp_ij) * grad_tail;
		   grad_2 = 2 * rp_i * rp_j;
		   if (comp_particleEnergy) {
		     particleEnergy[i] += value/2.0; particleEnergy[j] += value/2.0; }
		   if (comp_energy) {
		     *energy +=value; }
		   tmpx = grad_1 * rr * drx - tmp_2 * 
		     (grad_2 * rr * drx - rp_i * rr * p_indx[j] 
		      - rp_j * rr * p_indx[i]);
		   tmpy = grad_1 * rr * dry - tmp_2 * 
		     (grad_2 * rr * dry - rp_i * rr * p_indy[j] 
		      - rp_j * rr * p_indy[i]);
		   tmpz = grad_1 * rr * drz - tmp_2 * 
		     (grad_2 * rr * drz - rp_i * rr * p_indz[j] 
		      - rp_j * rr * p_indz[i]);
		   if (comp_force) {
		   force[i*DIM  ] -= tmpx; force[j*DIM  ] += tmpx;
		   force[i*DIM+1] -= tmpy; force[j*DIM+1] += tmpy;
		   force[i*DIM+2] -= tmpz; force[j*DIM+2] += tmpz;
		   }
		 }
	       }
	     } }
     } // FOURTH LOOP END

   // FIFTH LOOP: self energy contributions and sum-up force contributions
   //printf("Fifth loop\n");
   double qq, pp;
   for (i = 0; i < *nAtoms; i++)
     {
       typ1 = particleSpecies[i];
       // self energy contributions
       qq = charge[typ1-1]*charge[typ1-1];
       value = DP_EPS * dp_kappa * qq / sqrt(M_PI);
       if (comp_particleEnergy) {
	 particleEnergy[i] -= value; }
       *energy -=value;
       if (dp_alpha[typ1-1] != 0.0) {
	 pp = p_indx[i]*p_indx[i]
	   + p_indy[i]*p_indy[i] + p_indz[i]*p_indz[i];
	 value = pp / (2 * dp_alpha[typ1-1]);
	 if (comp_particleEnergy) {
	   particleEnergy[i] += value; }
	 *energy +=value;
       }
     } // FIFTH LOOP END
   
   /* everything is great */
   status = KIM_STATUS_OK;
   return;
}

/* Initialization function */
void model_init(void *km)
{
   /* Local variables */
   intptr_t* pkim = *((intptr_t**) km);
   double* model_cutoff;
   int ier;

   int* natom;
   double *dp_alpha, *dp_b, *dp_c;
   double *E_statx, *E_staty, *E_statz;
   double *E_indx, *E_indy, *E_indz;
   double *E_oldx, *E_oldy, *E_oldz;
   double *E_totx, *E_toty, *E_totz;
   double *p_srx, *p_sry, *p_srz;
   double *p_indx, *p_indy, *p_indz;
   double *ratio, *charge, last_charge, dp_kappa;
   double *buck_a, *buck_s, *buck_c;
   int sw_kappa;
   bool initialize;
   int ntype, nptype;
   int nparticle;
   struct model_buffer* buffer;

   /* store pointer to compute function in KIM object */
   ier = KIM_API_set_method(pkim, "compute", 1, (func_ptr) &compute);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_set_method", ier);
      exit(1);
   }

   /* store model cutoff in KIM object */
   model_cutoff = (double*) KIM_API_get_data(pkim, "cutoff", &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", ier);
      exit(1);
   }
   *model_cutoff = MODEL_CUTOFF; /* cutoff distance in angstroms */
   //printf("%f\n",*model_cutoff);

   // Dipole
   ntype = 3; nptype = 6;
   //natom = (int*) KIM_API_get_data(pkim, "numberOfParticles", &ier);
   natom = (int*) KIM_API_get_data(pkim, "numberOfParticles", &ier);

   dp_alpha = new double[ntype];
   dp_b     = new double[nptype];
   dp_c     = new double[nptype];
   E_statx  = new double[*natom];
   E_staty  = new double[*natom];
   E_statz  = new double[*natom];
   E_indx   = new double[*natom];
   E_indy   = new double[*natom];
   E_indz   = new double[*natom];
   E_oldx   = new double[*natom];
   E_oldy   = new double[*natom];
   E_oldz   = new double[*natom];
   E_totx   = new double[*natom];
   E_toty   = new double[*natom];
   E_totz   = new double[*natom];
   p_srx    = new double[*natom];
   p_sry    = new double[*natom];
   p_srz    = new double[*natom];
   p_indx   = new double[*natom];
   p_indy   = new double[*natom];
   p_indz   = new double[*natom];
   charge   = new double[ntype];
   buck_a   = new double[nptype];
   buck_s   = new double[nptype];
   buck_c   = new double[nptype];
   buffer = (struct model_buffer*) malloc(sizeof(struct model_buffer));

   buffer->nparticle = *natom;
   printf("---Initialize--- %d\n",buffer->nparticle);

   //Parameter for atom == 0:Zr, 1:O, 2:Y
   //Parameter for pair == 0:Zr-O, 1:O-O 2:Y-O (OLD)
   //Parameter for pair == 0:Zr-Zr, 1:Zr-O 2:Zr-Y, 3:O-O, 4:O-Y, 5:Y-Y (NEW)
   dp_kappa = 0.3;
   charge[0] = 1.909694;
   charge[1] =-0.954847;
   charge[2] = 1.432271;
   dp_alpha[0] = 0.100000;
   dp_alpha[1] = 0.045925;
   dp_alpha[2] = -0.000039;
   dp_b[0] = 179.571341;
   dp_b[1] = 7.621117;
   dp_b[2] = 292.814067;
   dp_b[3] = 126.606590;
   dp_b[4] = 147.375709;
   dp_b[5] = 290.957563;
   dp_c[0] = -173.277831;
   dp_c[1] = -197.989171;
   dp_c[2] = -67.427023;
   dp_c[3] = -195.776039;
   dp_c[4] = -223.136948;
   dp_c[5] = -341.661260;
   buck_a[0] = 4605.262793;
   buck_a[1] = 5025.220932;
   buck_a[2] = 26487.268870;
   buck_a[3] = 1125.669874;
   buck_a[4] = 12306.103355;
   buck_a[5] = 8291.238373;
   buck_s[0] = 0.372011;
   buck_s[1] = 0.247792;
   buck_s[2] = 0.330892;
   buck_s[3] = 0.279490;
   buck_s[4] = 0.225804;
   buck_s[5] = 0.415872;
   buck_c[0] = 213777.333397;
   buck_c[1] = 299999.999988;
   buck_c[2] = 1361734.175424;
   buck_c[3] = -53727.892471;
   buck_c[4] = 418385.411891;
   buck_c[5] = 890858.245624;

   buffer->dp_alpha = dp_alpha;
   buffer->dp_b = dp_b;
   buffer->dp_c = dp_c;
   buffer->E_statx = E_statx;
   buffer->E_staty = E_staty;
   buffer->E_statz = E_statz;
   buffer->E_indx = E_indx;
   buffer->E_indy = E_indy;
   buffer->E_indz = E_indz;
   buffer->E_oldx = E_oldx;
   buffer->E_oldy = E_oldy;
   buffer->E_oldz = E_oldz;
   buffer->E_totx = E_totx;
   buffer->E_toty = E_toty;
   buffer->E_totz = E_totz;
   buffer->p_srx = p_srx;
   buffer->p_sry = p_sry;
   buffer->p_srz = p_srz;
   buffer->p_indx = p_indx;
   buffer->p_indy = p_indy;
   buffer->p_indz = p_indz;
   buffer->ratio = ratio;
   buffer->charge = charge;
   buffer->last_charge = last_charge;
   buffer->dp_kappa = dp_kappa;
   buffer->buck_a = buck_a;
   buffer->buck_s = buck_s;
   buffer->buck_c = buck_c;
   buffer->sw_kappa = sw_kappa;
   buffer->initialize = initialize;
   buffer->ntype = ntype;
   buffer->nptype = nptype;

   KIM_API_set_model_buffer(pkim, (void*) buffer, &ier);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_set_model_buffer", ier);
   }

  return;
}
