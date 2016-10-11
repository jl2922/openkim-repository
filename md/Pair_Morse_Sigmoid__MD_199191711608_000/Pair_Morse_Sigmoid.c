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
 * Copyright (c) 2013--2015, Regents of the University of Minnesota.
 * All rights reserved.
 *
 * Contributors:
 *    Ryan S. Elliott
 *    Ellad B. Tadmor
 *    Valeriu Smirichinski
 *    Stephen M. Whalen
 *    Hao Xu
 */

/*******************************************************************************
 *
 *  Pair_Morse_Sigmoid
 *
 *  Morse pair potential KIM Model Driver
 *
 *  cutoff function is:
 *  phi(r) = epsilon * ( - exp[ -2*C*(r - Rzero) ] + 2*exp[ -C*(r - Rzero) ] )
 *                   * S(r, logisticWidth, logisticR0)
 *  where S(r, logisticWidth, logisticR0)
 *          =  1/(1 + exp( -logisticWidth*[r - logisticR0]))
 *
 *  logisticR0 = cutoff + 7/logisticWidth
 *
 *  Language: C
 *
 ******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "KIM_API_C.h"
#include "KIM_API_status.h"

/******************************************************************************
 * Below are the definitions and values of all Model parameters
 ******************************************************************************/
#define DIM 3       /* dimensionality of space */
#define SPECCODE 1  /* internal species code */
#define GORDER 3    /* number of coeffs for g(r) */

/* Define prototypes for Model Driver init */
/**/
int model_driver_init(void* km, char* paramfile_names, int* nmstrlen,
                      int* numparamfiles);

/* Define prototypes for Model (Driver) reinit, compute, and destroy */
/* defined as static to avoid namespace clashes with other Models    */
/**/
static int reinit(void* km);
static int destroy(void* km);
static int compute(void* km);
/**/
static void calc_phi(double const epsilon,
                     double const C,
                     double const Rzero,
                     double const logisticWidth,
                     double const cutoff,
                     double const r,
                     double* const phi);
static void calc_phi_dphi(double const epsilon,
                          double const C,
                          double const Rzero,
                          double const logisticWidth,
                          double const cutoff,
                          double const r,
                          double* const phi,
                          double* const dphi);
static void calc_phi_d2phi(double const epsilon,
                           double const C,
                           double const Rzero,
                           double const logisticWidth,
                           double const cutoff,
                           double const r,
                           double* const phi,
                           double* const dphi,
                           double* const d2phi);

/* Define model_buffer structure */
struct model_buffer {
  int NBC;
  int HalfOrFull;
  int IterOrLoca;
  int energy_ind;
  int forces_ind;
  int particleEnergy_ind;
  int process_dEdr_ind;
  int process_d2Edr2_ind;
  int model_index_shift;
  int numberOfParticles_ind;
  int particleSpecies_ind;
  int coordinates_ind;
  int numberContributingParticles_ind;
  int boxSideLengths_ind;
  int get_neigh_ind;
  int cutoff_ind;

  double Pcutoff;
  double transition;
  double cutsq;
  double epsilon;
  double C;
  double Rzero;
  double logisticWidth;
};


/* Calculate pair potential phi(r) */
static void calc_phi(double const epsilon,
                     double const C,
                     double const Rzero,
                     double const logisticWidth,
                     double const cutoff,
                     double const r,
                     double* const phi)
{
  /* local variables */
  double const logisticR0 = cutoff + 7.0/logisticWidth;
  double const S = 1.0/(1.0+exp(-logisticWidth*(r-logisticR0)));
  double const ep = exp(-C*(r-Rzero));
  double const ep2 = ep*ep;

  if (r < cutoff)
  {
    *phi = epsilon*( -ep2 + 2.0*ep )*S;
  }
  else
  {
    /* Argument exceeds cutoff radius */
    *phi = 0.0;
  }

  return;
}

/* Calculate pair potential phi(r) and its derivative dphi(r) */
static void calc_phi_dphi(double const epsilon,
                          double const C,
                          double const Rzero,
                          double const logisticWidth,
                          double const cutoff,
                          double const r,
                          double* const phi,
                          double* const dphi)
{
  /* local variables */
  double const logisticR0 = cutoff + 7.0/logisticWidth;
  double const S = 1.0/(1.0+exp(-logisticWidth*(r-logisticR0)));
  double const ep = exp(-C*(r-Rzero));
  double const ep2 = ep*ep;

  if (r < cutoff)
  {
    *phi  = epsilon*( -ep2 + 2.0*ep )*S;
    *dphi = 2.0*epsilon*C*( -ep + ep2 )*S
        + epsilon*( -ep2 + 2.0*ep )
        *S*S*logisticWidth*exp(-logisticWidth*(r-logisticR0));
  }
  else
  {
    /* Argument exceeds cutoff radius */
    *phi  = 0.0;
    *dphi = 0.0;
  }

  return;
}

/*
  Calculate pair potential phi(r) and its 1st & 2nd derivatives dphi(r),
  d2phi(r)
*/
static void calc_phi_d2phi(double const epsilon,
                           double const C,
                           double const Rzero,
                           double const logisticWidth,
                           double const cutoff,
                           double const r,
                           double* const phi,
                           double* const dphi,
                           double* const d2phi)
{
  /* local variables */
  double const logisticR0 = cutoff + 7.0/logisticWidth;
  double const S = 1.0/(1.0+exp(-logisticWidth*(r-logisticR0)));
  double const ep = exp(-C*(r-Rzero));
  double const ep2 = ep*ep;

  if (r < cutoff)
  {
    *phi   = epsilon*( -ep2 + 2.0*ep )*S;
    *dphi  = 2.0*epsilon*C*( -ep + ep2 )*S
        + epsilon*( -ep2 + 2.0*ep )
        *S*S*logisticWidth*exp(-logisticWidth*(r-logisticR0));
    *d2phi = 2.0*epsilon*C*C*(ep - 2.0*ep2)*S
        + 4.0*epsilon*C*( -ep + ep2 )
        *S*S*logisticWidth*exp(-logisticWidth*(r-logisticR0))
        + epsilon*( -ep2 + 2.0*ep )
        *(logisticWidth*logisticWidth
          *2*S*S*S*exp(-2*logisticWidth*(r-logisticR0))
          -logisticWidth*logisticWidth*S*S*exp(-logisticWidth*(r-logisticR0))
          );
  }
  else
  {
    /* Argument exceeds cutoff radius */
    *phi   = 0.0;
    *dphi  = 0.0;
    *d2phi = 0.0;
  }

  return;
}

/* compute function */
static int compute(void* km)
{
  /* local variables */
  intptr_t* pkim = *((intptr_t**) km);
  double R;
  double R_pairs[2];
  double* pR_pairs = &(R_pairs[0]);
  double Rsqij;
  double phi;
  double dphi;
  double d2phi;
  double dEidr;
  double d2Eidr;
  double Rij[DIM];
  double* pRij = &(Rij[0]);
  double Rij_pairs[2][3];
  double* pRij_pairs = &(Rij_pairs[0][0]);
  int ier;
  int i;
  int i_pairs[2];
  int* pi_pairs = &(i_pairs[0]);
  int j;
  int j_pairs[2];
  int* pj_pairs = &(j_pairs[0]);
  int jj;
  int k;
  int currentAtom;
  int* neighListOfCurrentAtom;
  int comp_energy;
  int comp_force;
  int comp_particleEnergy;
  int comp_process_dEdr;
  int comp_process_d2Edr2;
  int const zero = 0;
  int const one = 1;
  int request;

  double const* cutoff;
  int const* nAtoms;
  int const* particleSpecies;
  double const* Rij_list;
  double const* coords;
  double const* boxSideLengths;
  double* energy;
  double* force;
  double* particleEnergy;
  int const* numContrib;
  int numberContrib;
  int numOfAtomNeigh;
  typedef int (*get_neigh_ptr)(
      void *,int const*,int const*,int *, int *, int **, double const **);
  get_neigh_ptr get_neigh = NULL;

  /* get buf from KIM object */
  struct model_buffer const* const buf
      = (struct model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /*
    check to see if we have been asked to compute the forces, particleEnergy,
    and d1Edr
  */
  KIM_API_getm_compute_by_index(
      pkim, &ier, 5*3,
      buf->energy_ind,         &comp_energy,         1,
      buf->forces_ind,         &comp_force,          1,
      buf->particleEnergy_ind, &comp_particleEnergy, 1,
      buf->process_dEdr_ind,   &comp_process_dEdr,   1,
      buf->process_d2Edr2_ind, &comp_process_d2Edr2, 1);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute_by_index",
                         ier);
    return ier;
  }

  KIM_API_getm_data_by_index(
      pkim, &ier, 9*3,
      buf->cutoff_ind,                      &cutoff,         1,
      buf->numberOfParticles_ind,           &nAtoms,         1,
      buf->particleSpecies_ind,             &particleSpecies,1,
      buf->coordinates_ind,                 &coords,         1,
      buf->numberContributingParticles_ind, &numContrib,     (buf->HalfOrFull==1),
      buf->boxSideLengths_ind,              &boxSideLengths, (buf->NBC==2),
      buf->energy_ind,                      &energy,         comp_energy,
      buf->forces_ind,                      &force,          comp_force,
      buf->particleEnergy_ind,              &particleEnergy, comp_particleEnergy);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data_by_index",
                         ier);
    return ier;
  }
  if (buf->NBC!=3)
  {
    get_neigh = (get_neigh_ptr)
        KIM_API_get_method_by_index(pkim, buf->get_neigh_ind, &ier);
    if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__,
                           "KIM_API_get_method_by_index", ier);
      return ier;
    }
  }

  if (buf->HalfOrFull == 1)
  {
    if (3 != buf->NBC) /* non-CLUSTER cases */
    {
      numberContrib = *numContrib;
    }
    else /* CLUSTER cases */
    {
      numberContrib = *nAtoms;
    }
  }
  else
  {  /* provide initialization even if not used */
    numberContrib = *nAtoms;
  }

  /* Check to be sure that the atom types are correct */
  /**/
  ier = KIM_STATUS_FAIL; /* assume an error */
  for (i = 0; i < *nAtoms; ++i)
  {
    if ( SPECCODE != particleSpecies[i])
    {
      KIM_API_report_error(__LINE__, __FILE__,
                           "Unexpected species type detected", ier);
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

  /* Initialize neighbor handling for CLUSTER NBC */
  if (3 == buf->NBC) /* CLUSTER */
  {
    neighListOfCurrentAtom = (int *) malloc((*nAtoms)*sizeof(int));
  }

  /* Initialize neighbor handling for Iterator mode */

  if (1 == buf->IterOrLoca)
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

  /* Compute energy and forces */

  /* loop over particles and compute enregy and forces */
  i = -1;
  while( 1 )
  {

    /* Set up neighbor list for next atom for all NBC methods */
    if (1 == buf->IterOrLoca) /* ITERATOR mode */
    {
      ier = (*get_neigh)(&pkim, &zero, &one, &currentAtom, &numOfAtomNeigh,
                         &neighListOfCurrentAtom, &Rij_list);
      /* the end of the list, terminate loop */
      if (KIM_STATUS_NEIGH_ITER_PAST_END == ier)
      {
        break;
      }
      if (KIM_STATUS_OK > ier) /* some sort of problem, exit */
      {
        KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", ier);
        return ier;
      }

      i = currentAtom + buf->model_index_shift;
    }
    else
    {
      i++;
      if (*nAtoms <= i) /* incremented past end of list, terminate loop */
      {
        break;
      }

      if (3 == buf->NBC) /* CLUSTER NBC method */
      {
        numOfAtomNeigh = *nAtoms - (i + 1);
        for (k = 0; k < numOfAtomNeigh; ++k)
        {
          neighListOfCurrentAtom[k] = i + k + 1 - buf->model_index_shift;
        }
        ier = KIM_STATUS_OK;
      }
      else
      {
        request = i - buf->model_index_shift;
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

    /* loop over the neighbors of atom i */
    for (jj = 0; jj < numOfAtomNeigh; ++ jj)
    {
      /* get neighbor ID */
      j = neighListOfCurrentAtom[jj] + buf->model_index_shift;

      /* compute relative position vector and squared distance */
      Rsqij = 0.0;
      for (k = 0; k < DIM; ++k)
      {
        if (0 != buf->NBC) /* all methods except NEIGH_RVEC */
        {
          Rij[k] = coords[j*DIM + k] - coords[i*DIM + k];
        }
        else          /* NEIGH_RVEC_F method */
        {
          Rij[k] = Rij_list[jj*DIM + k];
        }

        /* apply periodic boundary conditions if required */
        if (2 == buf->NBC)
        {
          if (abs(Rij[k]) > 0.5*boxSideLengths[k])
          {
            Rij[k] -= (Rij[k]/fabs(Rij[k]))*boxSideLengths[k];
          }
        }

        /* compute squared distance */
        Rsqij += Rij[k]*Rij[k];
      }

      /* compute energy and force */
      if (Rsqij < buf->cutsq) /* particles are interacting ? */
      {
        R = sqrt(Rsqij);
        if (comp_process_d2Edr2)
        {
          /* compute pair potential and its derivatives */
          calc_phi_d2phi(buf->epsilon,
                         buf->C,
                         buf->Rzero,
                         buf->logisticWidth,
                         *cutoff,
                         R, &phi, &dphi, &d2phi);

          /* compute dEidr */
          if ((1 == buf->HalfOrFull) && (j < numberContrib))
          {
            /* Half mode -- double contribution */
            dEidr = dphi;
            d2Eidr = d2phi;
          }
          else
          {
            /* Full mode -- regular contribution */
            dEidr = 0.5*dphi;
            d2Eidr = 0.5*d2phi;
          }
        }
        else if (comp_force || comp_process_dEdr)
        {
          /* compute pair potential and its derivative */
          calc_phi_dphi(buf->epsilon,
                        buf->C,
                        buf->Rzero,
                        buf->logisticWidth,
                        *cutoff,
                        R, &phi, &dphi);

          /* compute dEidr */
          if ((1 == buf->HalfOrFull) && (j < numberContrib))
          {
            /* Half mode -- double contribution */
            dEidr = dphi;
          }
          else
          {
            /* Full mode -- regular contribution */
            dEidr = 0.5*dphi;
          }
        }
        else
        {
          /* compute just pair potential */
          calc_phi(buf->epsilon,
                   buf->C,
                   buf->Rzero,
                   buf->logisticWidth,
                   *cutoff,
                   R, &phi);
        }

        /* contribution to energy */
        if (comp_particleEnergy)
        {
          particleEnergy[i] += 0.5*phi;
          /* if half list add energy for the other atom in the pair */
          if ((1 == buf->HalfOrFull) && (j < numberContrib))
          {
            particleEnergy[j] += 0.5*phi;
          }
        }
        if (comp_energy)
        {
          if ((1 == buf->HalfOrFull) && (j < numberContrib))
          {
            /* Half mode -- add v to total energy */
            *energy += phi;
          }
          else
          {
            /* Full mode -- add half v to total energy */
            *energy += 0.5*phi;
          }
        }

        /* contribution to process_dEdr */
        if (comp_process_dEdr)
        {
          ier = KIM_API_process_dEdr(km, &dEidr, &R, &pRij, &i, &j);
        }

        /* contribution to process_d2Edr2 */
        if (comp_process_d2Edr2)
        {
          R_pairs[0] = R_pairs[1] = R;
          Rij_pairs[0][0] = Rij_pairs[1][0] = Rij[0];
          Rij_pairs[0][1] = Rij_pairs[1][1] = Rij[1];
          Rij_pairs[0][2] = Rij_pairs[1][2] = Rij[2];
          i_pairs[0] = i_pairs[1] = i;
          j_pairs[0] = j_pairs[1] = j;

          ier = KIM_API_process_d2Edr2(km, &d2Eidr, &pR_pairs, &pRij_pairs,
                                       &pi_pairs, &pj_pairs);
        }

        /* contribution to forces */
        if (comp_force)
        {
          for (k = 0; k < DIM; ++k)
          {  /* accumulate force on atom i */
            force[i*DIM + k] += dEidr*Rij[k]/R;
            /* accumulate force on atom j */
            force[j*DIM + k] -= dEidr*Rij[k]/R;
          }
        }
      }
    } /* loop on jj */
  }    /* infinite while loop (terminated by break statements above */

  /* Free temporary storage */
  if (3 == buf->NBC)
  {
    free(neighListOfCurrentAtom);
  }

  /* everything is great */
  ier = KIM_STATUS_OK;

  return ier;
}

/* Initialization function */
int model_driver_init(void *km, char* paramfile_names, int* nmstrlen,
                      int* numparamfiles)
{
  /* KIM variables */
  intptr_t* pkim = *((intptr_t**) km);
  char* paramfile1name;

  /* Local variables */
  FILE* fid;
  double cutoff;
  double logisticWidth;
  double epsilon;
  double C;
  double Rzero;
  double* model_cutoff;
  int ier;
  struct model_buffer* buf;
  const char* NBCstr;

  /* set paramfile1name */
  if (*numparamfiles != 1)
  {
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__,
                         "Incorrect number of parameter files.", ier);
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
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_setm_data", ier);
    return ier;
  }

  /* Read in model parameters from parameter file */
  fid = fopen(paramfile1name, "r");
  if (fid == NULL)
  {
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__,
                         "Unable to open parameter file for Morse parameters",
                         ier);
    return ier;
  }

  ier = fscanf(fid, "%lf \n%lf \n%lf \n%lf \n%lf",
               &cutoff,         /* cutoff distance in angstroms */
               &logisticWidth,  /* width of sigmoid */
               &epsilon,        /* Morse epsilon in eV */
               &C,              /* Morse C in 1/Angstroms */
               &Rzero);         /* Morse Rzero in Angstroms */
  fclose(fid);

  /* check that we read the right number of parameters */
  if (5 != ier)
  {
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__,
                         "Unable to read all Morse parameters", ier);
    return ier;
  }

  /* convert to appropriate units */
  cutoff *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "ps",
                                        1.0, 0.0,  0.0, 0.0, 0.0, &ier);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit",
                         ier);
    return ier;
  }

  epsilon *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "ps",
                                         0.0, 1.0,  0.0, 0.0, 0.0, &ier);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit",
                         ier);
    return ier;
  }

  C *= KIM_API_convert_to_act_unit(pkim, "A",  "eV", "e", "K", "ps",
                                   -1.0, 0.0,  0.0, 0.0, 0.0, &ier);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit",
                         ier);
    return ier;
  }

  Rzero *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "ps",
                                       1.0, 0.0,  0.0, 0.0, 0.0, &ier);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit",
                         ier);
    return ier;
  }

  /* store model cutoff in KIM object */
  model_cutoff = (double*) KIM_API_get_data(pkim, "cutoff", &ier);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", ier);
    return ier;
  }
  *model_cutoff = cutoff;

  /* allocate buf */
  buf = (struct model_buffer*) malloc(sizeof(struct model_buffer));
  if (NULL == buf)
  {
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
    return ier;
  }

  /* setup buf */
  /* set value of parameters */
  buf->Pcutoff = *model_cutoff;
  buf->logisticWidth = logisticWidth;
  buf->cutsq = (*model_cutoff)*(*model_cutoff);
  buf->epsilon = epsilon;
  buf->C = C;
  buf->Rzero = Rzero;

  /* Determine neighbor list boundary condition (NBC) */
  ier = KIM_API_get_NBC_method(pkim, &NBCstr);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_NBC_method", ier);
    return ier;
  }
  if ((!strcmp("NEIGH_RVEC_H",NBCstr)) || (!strcmp("NEIGH_RVEC_F",NBCstr)))
  {
    buf->NBC = 0;
  }
  else if ((!strcmp("NEIGH_PURE_H",NBCstr)) || (!strcmp("NEIGH_PURE_F",NBCstr)))
  {
    buf->NBC = 1;
  }
  else if ((!strcmp("MI_OPBC_H",NBCstr)) || (!strcmp("MI_OPBC_F",NBCstr)))
  {
    buf->NBC = 2;
  }
  else if (!strcmp("CLUSTER",NBCstr))
  {
    buf->NBC = 3;
  }
  else
  {
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__, "Unknown NBC method", ier);
    return ier;
  }

  /* Determine if Half or Full neighbor lists are being used */
  /*****************************
   * HalfOrFull = 1 -- Half
   *            = 2 -- Full
   *****************************/
  if (KIM_API_is_half_neighbors(pkim, &ier))
  {
    buf->HalfOrFull = 1;
  }
  else
  {
    buf->HalfOrFull = 2;
  }

  /* determine neighbor list handling mode */
  if (buf->NBC != 3)
  {
    /*****************************
     * IterOrLoca = 1 -- Iterator
     *            = 2 -- Locator
     *****************************/
    buf->IterOrLoca = KIM_API_get_neigh_mode(pkim, &ier);
    if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh_mode",
                           ier);
      return ier;
    }
    if ((buf->IterOrLoca != 1) && (buf->IterOrLoca != 2))
    {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__,
                           "Unsupported IterOrLoca mode", ier);
      return ier;
    }
  }
  else
  {
    buf->IterOrLoca = 2;   /* for CLUSTER NBC */
  }

  buf->model_index_shift = KIM_API_get_model_index_shift(pkim);

  KIM_API_getm_index(
      pkim, &ier, 12*3,
      "cutoff",                      &(buf->cutoff_ind),                      1,
      "numberOfParticles",           &(buf->numberOfParticles_ind),           1,
      "particleSpecies",             &(buf->particleSpecies_ind),             1,
      "numberContributingParticles", &(buf->numberContributingParticles_ind), 1,
      "coordinates",                 &(buf->coordinates_ind),                 1,
      "get_neigh",                   &(buf->get_neigh_ind),                   1,
      "boxSideLengths",              &(buf->boxSideLengths_ind),              1,
      "energy",                      &(buf->energy_ind),                      1,
      "forces",                      &(buf->forces_ind),                      1,
      "particleEnergy",              &(buf->particleEnergy_ind),              1,
      "process_dEdr",                &(buf->process_dEdr_ind),                1,
      "process_d2Edr2",              &(buf->process_d2Edr2_ind),              1
                     );
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_index", ier);
    return ier;
  }
  /* end setup buf */

  /* store in model buf */
  KIM_API_set_model_buffer(pkim, (void*) buf, &ier);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_set_model_buffer", ier);
    return ier;
  }

  /* set pointers to parameters in KIM object */
  KIM_API_setm_data(
      pkim, &ier, 6*4,
      "PARAM_FREE_cutoff",        1, &(buf->Pcutoff),       1,
      "PARAM_FREE_logisticWidth", 1, &(buf->logisticWidth), 1,
      "PARAM_FREE_epsilon",       1, &(buf->epsilon),       1,
      "PARAM_FREE_C",             1, &(buf->C),             1,
      "PARAM_FREE_Rzero",         1, &(buf->Rzero),         1,
      "PARAM_FIXED_cutsq",        1, &(buf->cutsq),         1
                    );
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_setm_data", ier);
    return ier;
  }

  return KIM_STATUS_OK;
}

/* Reinitialization function */
static int reinit(void *km)
{
  /* Local variables */
  intptr_t* pkim = *((intptr_t**) km);
  int ier;
  double *cutoff;
  struct model_buffer* buf;

  /* get buf from KIM object */
  buf = (struct model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* set new values in KIM object     */
  /*                                  */
  /* store model cutoff in KIM object */
  cutoff = KIM_API_get_data(pkim, "cutoff", &ier);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", ier);
    return ier;
  }
  *cutoff = buf->Pcutoff;

  /* set value of parameter cutsq */
  buf->cutsq = (*cutoff)*(*cutoff);

  ier = KIM_STATUS_OK;
  return ier;
}

/* destroy function */
static int destroy(void *km)
{
  /* Local variables */
  intptr_t* pkim = *((intptr_t**) km);
  struct model_buffer* buf;
  int ier;

  /* get model buf from KIM object */
  buf = (struct model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* destroy the buf */
  free(buf);

  ier = KIM_STATUS_OK;
  return ier;
}
