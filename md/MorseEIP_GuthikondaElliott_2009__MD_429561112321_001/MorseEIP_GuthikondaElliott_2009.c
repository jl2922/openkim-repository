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
 * Copyright (c) 2013--2016, Regents of the University of Minnesota.
 * All rights reserved.
 *
 * Contributors:
 *    Ryan S. Elliott
 *    Andrew Akerson
 *    Ellad B. Tadmor
 *    Valeriu Smirichinski
 *    Stephen M. Whalen
 *
 */


/*******************************************************************************
 *
 *  MorseEIP_GuthikondaElliott_2009
 *
 *  Temperature-dependent Morse pair potential KIM Model Driver
 *  Shifted to have zero energy at the cutoff radius
 *
 *  Language: C
 *
 ******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include "KIM_API_C.h"
#include "KIM_API_status.h"

/******************************************************************************
 * Below are the definitions and values of all Model parameters
 ******************************************************************************/
#define DIM 3  /* dimensionality of space */
#define MAXLINE 1024  /* max characters in line */

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
static void calc_phi(double const* const A,
                     double const* const B,
                     double const* const rHat,
                     double const* const shift,
                     double const* const cutoff,
                     double const r,
                     double* const phi);

static void calc_phi_dphi(double const* const A,
                          double const* const B,
                          double const* const rHat,
                          double const* const shift,
                          double const* const cutoff,
                          double const r,
                          double* const phi,
                          double* const dphi);

static void calc_phi_d2phi(double const* const A,
                           double const* const B,
                           double const* const rHat,
                           double const* const shift,
                           double const* const cutoff,
                           double const r,
                           double* const phi,
                           double* const dphi,
                           double* const d2phi);

double calc_A(double const* const A1,
              double const* const A2,
              double const* const A3,
              double const theta);

double calc_B(double const* const B1,
              double const* const B2,
              double const* const B3,
              double const theta);

double calc_rHat(double const* const r1,
                 double const* const r2,
                 double const* const r3,
                 double const theta);

double** AllocateAndInitialize2DArray(int const extentZero,
                                      int const extentOne);

double* AllocateAndInitialize1DArray(int const numberModelSpecies);
void Deallocate2DArrays(int const n, ...);
void Deallocate1DArrays(int const n, ...);
void getNextDataLine(FILE* const filePtr, char* nextLinePtr,
                     int const maxSize, int* const endOfFileFlag);

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

  double temperature;

  double* cutoffs;
  double* A1s;
  double* A2s;
  double* A3s;
  double* B1s;
  double* B2s;
  double* B3s;
  double* r1s;
  double* r2s;
  double* r3s;

  double** cutsq2D;
  double** As2D;
  double** Bs2D;
  double** rHats2D;
  double** shifts2D;
};

/* Calculate Guthikonda Elliott Paramaters */
double calc_A(double const* const A1,
              double const* const A2,
              double const* const A3,
              double const theta)
{
  double A;

  A = *A1 + (*A2)*(pow(theta, *A3) - 1.0);
  return A;
}

double calc_B(double const* const B1,
              double const* const B2,
              double const* const B3,
              double const theta)
{
  double B;

  B = *B1 + (*B2)*(pow(theta, *B3) - 1.0);
  return B;
}

double calc_rHat(double const* const r1,
                 double const* const r2,
                 double const* const r3,
                 double const theta)
{
  double rHat;

  rHat = *r1 + (*r2)*(exp((*r3)*(theta - 1.0)) - 1.0);
  return rHat;
}

/* Calculate pair potential phi(r) */
static void calc_phi(double const* const A,
                     double const* const B,
                     double const* const rHat,
                     double const* const shift,
                     double const* const cutoff,
                     double const r,
                     double* const phi)
{
  /* local variables */
  double ep;
  double ep2;
  double epsilon;
  double C;

  epsilon = -(*A);
  C = (*B)/(*rHat);

  ep  = exp(-(C)*(r-*rHat));
  ep2 = ep*ep;

  if (r > *cutoff)
  {
    /* Argument exceeds cutoff radius */
    *phi = 0.0;
  }
  else
  {
    *phi   = (epsilon)*( -ep2 + 2.0*ep ) + *shift;
  }

  return;
}

/* Calculate pair potential phi(r) and its derivative dphi(r) */
static void calc_phi_dphi(double const* const A,
                          double const* const B,
                          double const* const rHat,
                          double const* const shift,
                          double const* const cutoff,
                          double const r,
                          double* const phi,
                          double* const dphi)
{
  /* local variables */
  double ep;
  double ep2;
  double epsilon;
  double C;

  epsilon = -(*A);
  C = (*B)/(*rHat);

  ep  = exp(-(C)*(r-*rHat));
  ep2 = ep*ep;

  if (r > *cutoff)
  {
    /* Argument exceeds cutoff radius */
    *phi  = 0.0;
    *dphi = 0.0;
  }
  else
  {
    *phi  = (epsilon)*( -ep2 + 2.0*ep ) + *shift;
    *dphi = 2.0*(epsilon)*(C)*( -ep + ep2 );
  }

  return;
}

/*
  Calculate pair potential phi(r) and its 1st & 2nd derivatives dphi(r),
  d2phi(r)
*/
static void calc_phi_d2phi(double const* const A,
                           double const* const B,
                           double const* const rHat,
                           double const* const shift,
                           double const* const cutoff,
                           double const r,
                           double* const phi,
                           double* const dphi,
                           double* const d2phi)
{
  /* local variables */
  double ep;
  double ep2;
  double epsilon;
  double C;

  epsilon = -(*A);
  C = (*B)/(*rHat);
  ep  = exp(-(C)*(r-*rHat));
  ep2 = ep*ep;

  if (r > *cutoff)
  {
    /* Argument exceeds cutoff radius */
    *phi   = 0.0;
    *dphi  = 0.0;
    *d2phi = 0.0;
  }
  else
  {
    *phi   = (epsilon)*( -ep2 + 2.0*ep ) + *shift;
    *dphi  = 2.0*(epsilon)*(C)*( -ep + ep2 );
    *d2phi = 2.0*(epsilon)*(C)*(C)*(ep - 2.0*ep2);
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
  double *pR_pairs = &(R_pairs[0]);
  double Rsqij;
  double phi;
  double dphi;
  double d2phi;
  double dEidr;
  double d2Eidr;
  double Rij[DIM];
  double *pRij = &(Rij[0]);
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
  int iSpecies, jSpecies;
  int currentAtom;
  int* neighListOfCurrentAtom;
  struct model_buffer* buffer;
  int comp_energy;
  int comp_force;
  int comp_particleEnergy;
  int comp_process_dEdr;
  int comp_process_d2Edr2;
  int NBC;
  int HalfOrFull;
  int IterOrLoca;
  int model_index_shift;
  int zero = 0;
  int one = 1;
  int request;

  int* nAtoms;
  int* particleSpecies;
  double* cutoff;
  double ijcutoff;
  double** cutsq2D;
  double** As2D;
  double** Bs2D;
  double** rHats2D;
  double** shifts2D;
  double* Rij_list;
  double* coords;
  double* energy;
  double* force;
  double* particleEnergy;
  double* boxSideLengths;
  int* numContrib;
  int numberContrib;
  int numOfAtomNeigh;
  typedef int (*get_neigh_ptr)(void *,int *,int *,int *, int *, int **,
                               double **);
  get_neigh_ptr get_neigh = NULL;


  /* get buffer from KIM object */
  buffer = (struct model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* unpack info from the buffer */
  NBC = buffer->NBC;
  HalfOrFull = buffer->HalfOrFull;
  IterOrLoca = buffer->IterOrLoca;
  model_index_shift = buffer->model_index_shift;
  cutsq2D = (buffer->cutsq2D);
  As2D = (buffer->As2D);
  Bs2D = (buffer->Bs2D);
  rHats2D = (buffer->rHats2D);
  shifts2D = (buffer->shifts2D);

  /*
    check to see if we have been asked to compute the forces, particleEnergy,
    and d1Edr
  */
  KIM_API_getm_compute_by_index(
      pkim, &ier, 5*3,
      buffer->energy_ind,         &comp_energy,         1,
      buffer->forces_ind,         &comp_force,          1,
      buffer->particleEnergy_ind, &comp_particleEnergy, 1,
      buffer->process_dEdr_ind,   &comp_process_dEdr,   1,
      buffer->process_d2Edr2_ind, &comp_process_d2Edr2, 1);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute_by_index",
                         ier);
    return ier;
  }

  KIM_API_getm_data_by_index(
      pkim, &ier, 9*3,
      buffer->cutoff_ind,                      &cutoff,          1,
      buffer->numberOfParticles_ind,           &nAtoms,          1,
      buffer->particleSpecies_ind,             &particleSpecies, 1,
      buffer->coordinates_ind,                 &coords,          1,
      buffer->numberContributingParticles_ind, &numContrib,     (HalfOrFull==1),
      buffer->boxSideLengths_ind,              &boxSideLengths, (NBC==2),
      buffer->energy_ind,                      &energy,         comp_energy,
      buffer->forces_ind,                      &force,          comp_force,
      buffer->particleEnergy_ind,         &particleEnergy, comp_particleEnergy);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data_by_index",
                         ier);
    return ier;
  }
  if (NBC!=3)
  {
    get_neigh = (get_neigh_ptr)
        KIM_API_get_method_by_index(pkim, buffer->get_neigh_ind, &ier);
    if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__,
                           "KIM_API_get_method_by_index", ier);
      return ier;
    }
  }

  if (HalfOrFull == 1)
  {
    if (3 != NBC)  /* non-CLUSTER cases */
    {
      numberContrib = *numContrib;
    }
    else  /* CLUSTER cases */
    {
      numberContrib = *nAtoms;
    }
  }
  else
  { /* provide initialization even if not used */
    numberContrib = *nAtoms;
  }

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

  /* Compute energy and forces */

  /* loop over particles and compute enregy and forces */
  i = -1;
  while( 1 )
  {

    /* Set up neighbor list for next atom for all NBC methods */
    if (1 == IterOrLoca)  /* ITERATOR mode */
    {
      ier = (*get_neigh)(&pkim, &zero, &one, &currentAtom, &numOfAtomNeigh,
                         &neighListOfCurrentAtom, &Rij_list);
      /* the end of the list, terminate loop */
      if (KIM_STATUS_NEIGH_ITER_PAST_END == ier)
      {
        break;
      }
      if (KIM_STATUS_OK > ier)  /* some sort of problem, exit */
      {
        KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", ier);
        return ier;
      }

      i = currentAtom + model_index_shift;
    }
    else
    {
      i++;
      if (*nAtoms <= i)  /* incremented past end of list, terminate loop */
      {
        break;
      }

      if (3 == NBC)  /* CLUSTER NBC method */
      {
        numOfAtomNeigh = *nAtoms - (i + 1);
        for (k = 0; k < numOfAtomNeigh; ++k)
        {
          neighListOfCurrentAtom[k] = i + k + 1 - model_index_shift;
        }
        ier = KIM_STATUS_OK;
      }
      else
      {
        request = i - model_index_shift;
        ier = (*get_neigh)(&pkim, &one, &request,
                           &currentAtom, &numOfAtomNeigh,
                           &neighListOfCurrentAtom, &Rij_list);
        if (KIM_STATUS_OK != ier)  /* some sort of problem, exit */
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
      /* get neighbor ID */
      j = neighListOfCurrentAtom[jj] + model_index_shift;
      jSpecies = particleSpecies[j];
      /* compute relative position vector and squared distance */
      Rsqij = 0.0;
      for (k = 0; k < DIM; ++k)
      {
        if (0 != NBC)  /* all methods except NEIGH_RVEC */
        {
          Rij[k] = coords[j*DIM + k] - coords[i*DIM + k];
        }
        else  /* NEIGH_RVEC_F method */
        {
          Rij[k] = Rij_list[jj*DIM + k];
        }

        /* apply periodic boundary conditions if required */
        if (2 == NBC)
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
      if (Rsqij < cutsq2D[iSpecies][jSpecies])  /* particles are interacting? */
      {
        ijcutoff = sqrt(cutsq2D[iSpecies][jSpecies]);
        R = sqrt(Rsqij);
        if (comp_process_d2Edr2)
        {
          /* compute pair potential and its derivatives */
          calc_phi_d2phi(&As2D[iSpecies][jSpecies],
                         &Bs2D[iSpecies][jSpecies],
                         &rHats2D[iSpecies][jSpecies],
                         &shifts2D[iSpecies][jSpecies],
                         &ijcutoff, R, &phi, &dphi, &d2phi);

          /* compute dEidr */
          if ((1 == HalfOrFull) && (j < numberContrib))
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
          calc_phi_dphi(&(As2D[iSpecies][jSpecies]),
                        &(Bs2D[iSpecies][jSpecies]),
                        &(rHats2D[iSpecies][jSpecies]),
                        &(shifts2D[iSpecies][jSpecies]),
                        &ijcutoff, R, &phi, &dphi);

          /* compute dEidr */
          if ((1 == HalfOrFull) && (j < numberContrib))
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
          calc_phi(&(As2D[iSpecies][jSpecies]),
                   &(Bs2D[iSpecies][jSpecies]),
                   &(rHats2D[iSpecies][jSpecies]),
                   &(shifts2D[iSpecies][jSpecies]),
                   &ijcutoff, R, &phi);
        }

        /* contribution to energy */
        if (comp_particleEnergy)
        {
          particleEnergy[i] += 0.5*phi;
          /* if half list add energy for the other atom in the pair */
          if ((1 == HalfOrFull) && (j < numberContrib))
          {
            particleEnergy[j] += 0.5*phi;
          }
        }
        if (comp_energy)
        {
          if ((1 == HalfOrFull) && (j < numberContrib))
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
          { /* accumulate force on atom i */
            force[i*DIM + k] += dEidr*Rij[k]/R;
            /* accumulate force on atom j */
            force[j*DIM + k] -= dEidr*Rij[k]/R;
          }
        }
      }
    }  /* loop on jj */
  }  /* infinite while loop (terminated by break statements above */

  /* Free temporary storage */
  if (3 == NBC)
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
  double* cutoffs;
  double *A1s, *A2s, *A3s;
  double *B1s, *B2s, *B3s;
  double *r1s, *r2s, *r3s;

  int ier;
  int i, j;
  int dummyint;
  struct model_buffer* buffer;
  const char* NBCstr;
  const char* kimModelParticleSpecies;
  const char** particleNames;

  int numberModelSpecies;
  int N;
  int endOfFileFlag;
  int fileLineCount;
  int iIndex, jIndex, indx;
  int numberUniqueSpeciesPairs;
  char spec1[MAXLINE], spec2[MAXLINE], nextLine[MAXLINE];
  char *nextLinePtr;
  double initialTemp;
  double nextCutoff;
  double nextA1, nextA2, nextA3;
  double nextB1, nextB2, nextB3;
  double nextr1, nextr2, nextr3;

  nextLinePtr = nextLine;
  KIM_API_get_num_model_species(pkim, &numberModelSpecies, &dummyint);

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

  /* Open Parameter File */

  /* set paramfile1name */
  if (*numparamfiles != 1)
  {
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__,
                         "Incorrect number of parameter files.", ier);
    return ier;
  }
  paramfile1name = paramfile_names;
  fid = fopen(paramfile1name, "r");
  if (fid == NULL)
  {
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__,
                         "Unable to open parameter file for Morse parameters",
                         ier);
    return ier;
  }

  /* Read line 0 of Parameter File */
  getNextDataLine(fid, nextLinePtr, MAXLINE, &endOfFileFlag);
  ier = sscanf(nextLine, "%d %lg", &N, &initialTemp);
  if (ier != 2)
  {
    sprintf(nextLine, "unable to read first line of the parameter file");
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__, nextLine, ier);
    fclose(fid);
    return ier;
  }
  if (N != numberModelSpecies)
  {
    sprintf(nextLine, "The value for N from the parameter file is inconsistent "
            "with numberModelSpecies");
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__, nextLine, ier);
    fclose(fid);
    return ier;
  }

  /* get and correctly order the particle names */
  particleNames = malloc(numberModelSpecies*sizeof(kimModelParticleSpecies));
  for (i = 0; i < numberModelSpecies; ++i)
  {
    ier = KIM_API_get_model_species(pkim, i, &kimModelParticleSpecies);
    if (ier < KIM_STATUS_OK)
    {
      KIM_API_report_error(__LINE__, __FILE__, "get_model_species", ier);
      free(particleNames);
      fclose(fid);
      return ier;
    }
    indx = KIM_API_get_species_code(pkim, kimModelParticleSpecies, &ier);
    if (indx >= numberModelSpecies)
    {
      KIM_API_report_error(__LINE__, __FILE__, "get_species_code",
                           KIM_STATUS_FAIL);
      free(particleNames);
      fclose(fid);
      return KIM_STATUS_FAIL;
    }
    particleNames[indx] = kimModelParticleSpecies;
  }

  cutoffs = AllocateAndInitialize1DArray(numberModelSpecies);
  A1s = AllocateAndInitialize1DArray(numberModelSpecies);
  A2s = AllocateAndInitialize1DArray(numberModelSpecies);
  A3s = AllocateAndInitialize1DArray(numberModelSpecies);
  B1s = AllocateAndInitialize1DArray(numberModelSpecies);
  B2s = AllocateAndInitialize1DArray(numberModelSpecies);
  B3s = AllocateAndInitialize1DArray(numberModelSpecies);
  r1s = AllocateAndInitialize1DArray(numberModelSpecies);
  r2s = AllocateAndInitialize1DArray(numberModelSpecies);
  r3s = AllocateAndInitialize1DArray(numberModelSpecies);

  numberUniqueSpeciesPairs = ((numberModelSpecies+1)*numberModelSpecies)/2;
  indx = 0;
  /* set all values of cutoffs1D to -1 for check later */
  for (i = 0; i<numberUniqueSpeciesPairs; i++)
  {
    cutoffs[i] = -1.0;
  }

  fileLineCount = 0;
  /* Read and process data lines. */
  getNextDataLine(fid, nextLinePtr, MAXLINE, &endOfFileFlag);
  while (endOfFileFlag == 0)
  {
    fileLineCount++;
    if (1 == fileLineCount % 3)
    {
      ier = sscanf(nextLine, "%s  %s  %lf  %lf  %lf  %lf",
                   spec1, spec2, &nextCutoff, &nextA1, &nextA2, &nextA3);
      if (ier != 6)
      {
        sprintf(nextLine, "error reading lines of the parameter file");
        KIM_API_report_error(__LINE__, __FILE__, nextLine, KIM_STATUS_FAIL);
        fclose(fid);
        free(particleNames);
        Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                           r1s, r2s, r3s);
        return KIM_STATUS_FAIL;
      }
      iIndex = jIndex = -1;
      for (i = 0; i <  N; i++)
      {
        if (strcmp(spec1, particleNames[i]) == 0)
        {
          iIndex = i;
        }
        if (strcmp(spec2, particleNames[i]) == 0)
        {
          jIndex = i;
        }
      }
      if ((iIndex == -1) || (jIndex == -1))
      {
        sprintf(nextLine, "Unsupported Species name found in parameter file");
        KIM_API_report_error(__LINE__, __FILE__, nextLine, KIM_STATUS_FAIL);
        fclose(fid);
        free(particleNames);
        Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                           r1s, r2s, r3s);
        return KIM_STATUS_FAIL;
      }
      if (iIndex >= jIndex)
      {
        indx = jIndex*N + iIndex - (jIndex*jIndex + jIndex)/2;
      }
      else
      {
        indx = iIndex*N + jIndex - (iIndex*iIndex + iIndex)/2;
      }
      cutoffs[indx] = nextCutoff;
      A1s[indx] = nextA1;
      A2s[indx] = nextA2;
      A3s[indx] = nextA3;
    }
    else if(2 == fileLineCount % 3)
    {
      ier = sscanf(nextLine, "%lf %lf %lf",
                   &nextB1, &nextB2, &nextB3);
      if (ier != 3)
      {
        sprintf(nextLine, "error reading lines of the parameter file");
        KIM_API_report_error(__LINE__, __FILE__, nextLine, KIM_STATUS_FAIL);
        fclose(fid);
        free(particleNames);
        Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                           r1s, r2s, r3s);
        return KIM_STATUS_FAIL;
      }
      B1s[indx] = nextB1;
      B2s[indx] = nextB2;
      B3s[indx] = nextB3;
    }
    else
    {
      ier = sscanf(nextLine, "%lf %lf %lf",
                   &nextr1, &nextr2, &nextr3);
      if (ier != 3)
      {
        sprintf(nextLine, "error reading lines of the parameter file");
        KIM_API_report_error(__LINE__, __FILE__, nextLine, KIM_STATUS_FAIL);
        fclose(fid);
        free(particleNames);
        Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                           r1s, r2s, r3s);
        return KIM_STATUS_FAIL;
      }
      r1s[indx] = nextr1;
      r2s[indx] = nextr2;
      r3s[indx] = nextr3;
    }

    getNextDataLine(fid, nextLinePtr, MAXLINE, &endOfFileFlag);
  }
  /* close parameter file */
  fclose(fid);

  /* check that we got all pairs */
  ier = 0;
  sprintf(nextLine, "There are not values for the pairs: \n");
  for (i=0; i<numberModelSpecies; i++)
  {
    for (j = i; j<numberModelSpecies; j++)
    {
      indx = i*N + j - (i*i + i)/2;
      if (cutoffs[indx] == -1)
      {
        strcat(nextLine, particleNames[i]);
        strcat(nextLine, "     ");
        strcat(nextLine, particleNames[j]);
        strcat(nextLine, "\n");
        ier = 1;
      }
    }
  }
  free(particleNames);
  if (ier == 1)
  {
    KIM_API_report_error(__LINE__, __FILE__, nextLine, KIM_STATUS_FAIL);
    Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                       r1s, r2s, r3s);
    return KIM_STATUS_FAIL;
  }

  /* convert to appropriate units */
  for (i = 0; i < numberUniqueSpeciesPairs; i++)
  {
    cutoffs[i] *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "ps",
                                              1.0, 0.0,  0.0, 0.0, 0.0, &ier);
    if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit",
                           ier);
      Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                         r1s, r2s, r3s);
      return ier;
    }

    A1s[i] *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "ps",
                                          0.0, 1.0,  0.0, 0.0, 0.0, &ier);
    A2s[i] *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "ps",
                                          0.0, 1.0,  0.0, 0.0, 0.0, &ier);
    if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit",
                           ier);
      Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                         r1s, r2s, r3s);
      return ier;
    }

    r1s[i] *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "ps",
                                          1.0, 0.0,  0.0, 0.0, 0.0, &ier);
    r2s[i] *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "ps",
                                          1.0, 0.0,  0.0, 0.0, 0.0, &ier);
    if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit",
                           ier);
      Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                         r1s, r2s, r3s);
      return ier;
    }
  }

  /* allocate buffer */
  buffer = (struct model_buffer*) malloc(sizeof(struct model_buffer));
  if (NULL == buffer)
  {
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
    Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                       r1s, r2s, r3s);
    return ier;
  }

  /* setup buffer */
  buffer->cutsq2D
      = AllocateAndInitialize2DArray(numberModelSpecies, numberModelSpecies);
  buffer->As2D
      = AllocateAndInitialize2DArray(numberModelSpecies, numberModelSpecies);
  buffer->Bs2D
      = AllocateAndInitialize2DArray(numberModelSpecies, numberModelSpecies);
  buffer->rHats2D
      = AllocateAndInitialize2DArray(numberModelSpecies, numberModelSpecies);
  buffer->shifts2D
      = AllocateAndInitialize2DArray(numberModelSpecies, numberModelSpecies);

  /* set value of parameters */
  buffer->temperature = initialTemp;
  buffer->cutoffs = cutoffs;
  buffer->A1s = A1s;
  buffer->A2s = A2s;
  buffer->A3s = A3s;
  buffer->B1s = B1s;
  buffer->B2s = B2s;
  buffer->B3s = B3s;
  buffer->r1s = r1s;
  buffer->r2s = r2s;
  buffer->r3s = r3s;

  /* Determine neighbor list boundary condition (NBC) */
  ier = KIM_API_get_NBC_method(pkim, &NBCstr);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_NBC_method", ier);
    Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                       r1s, r2s, r3s);
    Deallocate2DArrays(5, buffer->cutsq2D, buffer->As2D, buffer->Bs2D,
                       buffer->rHats2D, buffer->shifts2D);
    free(buffer);
    return ier;
  }
  if ((!strcmp("NEIGH_RVEC_H",NBCstr)) || (!strcmp("NEIGH_RVEC_F",NBCstr)))
  {
    buffer->NBC = 0;
  }
  else if ((!strcmp("NEIGH_PURE_H",NBCstr)) || (!strcmp("NEIGH_PURE_F",NBCstr)))
  {
    buffer->NBC = 1;
  }
  else if ((!strcmp("MI_OPBC_H",NBCstr)) || (!strcmp("MI_OPBC_F",NBCstr)))
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
    Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                       r1s, r2s, r3s);
    Deallocate2DArrays(5, buffer->cutsq2D, buffer->As2D, buffer->Bs2D,
                       buffer->rHats2D, buffer->shifts2D);
    free(buffer);
    return ier;
  }

  /* Determine if Half or Full neighbor lists are being used */
  /*****************************
   * HalfOrFull = 1 -- Half
   *            = 2 -- Full
   *****************************/
  if (KIM_API_is_half_neighbors(pkim, &ier))
  {
    buffer->HalfOrFull = 1;
  }
  else
  {
    buffer->HalfOrFull = 2;
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
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh_mode",
                           ier);
      Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                         r1s, r2s, r3s);
      Deallocate2DArrays(5, buffer->cutsq2D, buffer->As2D, buffer->Bs2D,
                         buffer->rHats2D, buffer->shifts2D);
      free(buffer);
      return ier;
    }
    if ((buffer->IterOrLoca != 1) && (buffer->IterOrLoca != 2))
    {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__,
                           "Unsupported IterOrLoca mode", ier);
      Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                         r1s, r2s, r3s);
      Deallocate2DArrays(5, buffer->cutsq2D, buffer->As2D, buffer->Bs2D,
                         buffer->rHats2D, buffer->shifts2D);
      free(buffer);
      return ier;
    }
  }
  else
  {
    buffer->IterOrLoca = 2;   /* for CLUSTER NBC */
  }

  buffer->model_index_shift = KIM_API_get_model_index_shift(pkim);

  KIM_API_getm_index(
      pkim, &ier, 12*3,
      "cutoff",                      &(buffer->cutoff_ind),            1,
      "numberOfParticles",           &(buffer->numberOfParticles_ind), 1,
      "particleSpecies",             &(buffer->particleSpecies_ind),   1,
      "numberContributingParticles", &(buffer->numberContributingParticles_ind), 1,
      "coordinates",                 &(buffer->coordinates_ind),       1,
      "get_neigh",                   &(buffer->get_neigh_ind),         1,
      "boxSideLengths",              &(buffer->boxSideLengths_ind),    1,
      "energy",                      &(buffer->energy_ind),            1,
      "forces",                      &(buffer->forces_ind),            1,
      "particleEnergy",              &(buffer->particleEnergy_ind),    1,
      "process_dEdr",                &(buffer->process_dEdr_ind),      1,
      "process_d2Edr2",              &(buffer->process_d2Edr2_ind),    1
                     );
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_index", ier);
    Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                       r1s, r2s, r3s);
    Deallocate2DArrays(5, buffer->cutsq2D, buffer->As2D, buffer->Bs2D,
                       buffer->rHats2D, buffer->shifts2D);
    free(buffer);
    return ier;
  }
  /* end setup buffer */

  /* store in model buffer */
  KIM_API_set_model_buffer(pkim, (void*) buffer, &ier);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_set_model_buffer", ier);
    Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                       r1s, r2s, r3s);
    Deallocate2DArrays(5, buffer->cutsq2D, buffer->As2D, buffer->Bs2D,
                       buffer->rHats2D, buffer->shifts2D);
    free(buffer);
    return ier;
  }

  /* set pointers to parameters in KIM object */
  KIM_API_setm_data(
      pkim, &ier, 11*4,
      "PARAM_FREE_temperature",                     1,&(buffer->temperature), 1,
      "PARAM_FREE_cutoffs",  numberUniqueSpeciesPairs, (buffer->cutoffs),     1,
      "PARAM_FREE_A1s",      numberUniqueSpeciesPairs, (buffer->A1s),         1,
      "PARAM_FREE_A2s",      numberUniqueSpeciesPairs, (buffer->A2s),         1,
      "PARAM_FREE_A3s",      numberUniqueSpeciesPairs, (buffer->A3s),         1,
      "PARAM_FREE_B1s",      numberUniqueSpeciesPairs, (buffer->B1s),         1,
      "PARAM_FREE_B2s",      numberUniqueSpeciesPairs, (buffer->B2s),         1,
      "PARAM_FREE_B3s",      numberUniqueSpeciesPairs, (buffer->B3s),         1,
      "PARAM_FREE_r1s",      numberUniqueSpeciesPairs, (buffer->r1s),         1,
      "PARAM_FREE_r2s",      numberUniqueSpeciesPairs, (buffer->r2s),         1,
      "PARAM_FREE_r3s",      numberUniqueSpeciesPairs, (buffer->r3s),         1
                    );
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_setm_data", ier);
    Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                       r1s, r2s, r3s);
    Deallocate2DArrays(5, buffer->cutsq2D, buffer->As2D, buffer->Bs2D,
                       buffer->rHats2D, buffer->shifts2D);
    free(buffer);
    return ier;
  }
  /* Call the reinit Function */
  ier = reinit(km);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "reinit error", ier);
    Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                       r1s, r2s, r3s);
    Deallocate2DArrays(5, buffer->cutsq2D, buffer->As2D, buffer->Bs2D,
                       buffer->rHats2D, buffer->shifts2D);
    free(buffer);
    return ier;
  }

  return KIM_STATUS_OK;
}

/******************************************************************************/
/* Reinitialization function */
static int reinit(void *km)
{
  /* Local variables */
  intptr_t* pkim = *((intptr_t**) km);
  int ier;
  double cutoff;
  double *model_cutoff;
  double *cutoffs;
  double *A1s, *A2s, *A3s;
  double *B1s, *B2s, *B3s;
  double *r1s, *r2s, *r3s;
  double temperature;
  double theta;

  double nextShift;
  double dummy;
  int dummyint;
  int numberSpecies, maxStringLength;
  const char* simSpecies;
  int numberModelSpecies;
  int iIndex, jIndex, indx;
  int i, j;
  struct model_buffer* buffer;

  /* get buffer from KIM object */
  buffer = (struct model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* set value for 2D parameters */
  KIM_API_get_num_model_species(pkim, &numberModelSpecies, &dummyint);

  cutoffs = buffer->cutoffs;
  temperature = buffer->temperature;
  A1s = buffer->A1s;
  A2s = buffer->A2s;
  A3s = buffer->A3s;
  B1s = buffer->B1s;
  B2s = buffer->B2s;
  B3s = buffer->B3s;
  r1s = buffer->r1s;
  r2s = buffer->r2s;
  r3s = buffer->r3s;

  theta = temperature/333.15;
  for (i = 0; i < numberModelSpecies; ++i)
  {
    for (j = i; j < numberModelSpecies ; ++j)
    {
      indx = i*numberModelSpecies + j - (i*i + i)/2;
      buffer->cutsq2D[i][j] = buffer->cutsq2D[j][i]
          = (cutoffs[indx]*cutoffs[indx]);

      buffer->As2D[i][j] = buffer->As2D[j][i]
          = calc_A(&A1s[indx],
                   &A2s[indx],
                   &A3s[indx],
                   theta);

      buffer->Bs2D[i][j] = buffer->Bs2D[j][i]
          = calc_B(&B1s[indx],
                   &B2s[indx],
                   &B3s[indx],
                   theta);

      buffer->rHats2D[i][j] = buffer->rHats2D[j][i]
          = calc_rHat(&r1s[indx],
                      &r2s[indx],
                      &r3s[indx],
                      theta);
    }
  }

  /* store model cutoff in KIM object */
  cutoff = 0.0;
  ier = KIM_API_get_num_sim_species(pkim, &numberSpecies, &maxStringLength);
  if (ier < KIM_STATUS_OK)
  {
    KIM_API_report_error(__LINE__, __FILE__, "get_num_sim_species", ier);
    return ier;
  }
  for (i = 0; i < numberSpecies; i++)
  {
    for(j = i; j < numberSpecies; j++)
    {
      ier = KIM_API_get_sim_species( pkim, i, &simSpecies);
      if (ier < KIM_STATUS_OK)
      {
        KIM_API_report_error(__LINE__, __FILE__, "get_num_sim_species", ier);
        return ier;
      }
      iIndex = KIM_API_get_species_code(pkim, simSpecies, &ier);
      if (iIndex >= numberModelSpecies || ier<KIM_STATUS_OK )
      {
        KIM_API_report_error(__LINE__, __FILE__, "get_species_code",
                             KIM_STATUS_FAIL);
        return KIM_STATUS_FAIL;
      }
      ier = KIM_API_get_sim_species( pkim, j, &simSpecies);
      if (ier < KIM_STATUS_OK)
      {
        KIM_API_report_error(__LINE__, __FILE__, "get_num_sim_species", ier);
        return ier;
      }
      jIndex = KIM_API_get_species_code(pkim, simSpecies, &ier);
      if (jIndex >= numberModelSpecies || ier<KIM_STATUS_OK )
      {
        KIM_API_report_error(__LINE__, __FILE__, "get_species_code",
                             KIM_STATUS_FAIL);
        return KIM_STATUS_FAIL;
      }

      if (iIndex >= jIndex)
      {
        indx = jIndex*numberModelSpecies + iIndex - (jIndex*jIndex + jIndex)/2;
      }
      else
      {
        indx = iIndex*numberModelSpecies + jIndex - (iIndex*iIndex + iIndex)/2;
      }
      if (cutoff < cutoffs[indx])
      {
        cutoff = cutoffs[indx];
      }
    }
  }
  model_cutoff = (double*) KIM_API_get_data(pkim, "cutoff", &ier);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", ier);
    return ier;
  }
  *model_cutoff = cutoff;

  /* Set Values for Shifts */
  dummy = 0.0;
  cutoff += 1.0;  /* add a bit to avoid rounding problem */
  for (i = 0 ; i < numberModelSpecies; i++)
  {
    for (j = i; j < numberModelSpecies; j++)
    {
      indx = i*numberModelSpecies + j - (i*i + i)/2;
      /* call calc_phi with r=cutoff and shift=0.0 */
      calc_phi(&(buffer->As2D[i][j]),
               &(buffer->Bs2D[i][j]),
               &(buffer->rHats2D[i][j]),
               &dummy,
               &cutoff, cutoffs[indx], &nextShift);
      /* set shift to -shift */
      buffer->shifts2D[i][j] = buffer->shifts2D[j][i] = -nextShift;
    }
  }
  ier = KIM_STATUS_OK;
  return ier;
}

/******************************************************************************/
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

  /* destroy the buffer */
  Deallocate2DArrays(5,
                     buffer->cutsq2D,
                     buffer->As2D,
                     buffer->Bs2D,
                     buffer->rHats2D,
                     buffer->shifts2D);
  Deallocate1DArrays(10,
                     buffer->cutoffs,
                     buffer->A1s,
                     buffer->A2s,
                     buffer->A3s,
                     buffer->B1s,
                     buffer->B2s,
                     buffer->B3s,
                     buffer->r1s,
                     buffer->r2s,
                     buffer->r3s);
  free(buffer);

  ier = KIM_STATUS_OK;
  return ier;
}

/******************************************************************************/
void getNextDataLine( FILE* const filePtr, char* nextLinePtr,
                      int const maxSize, int* const endOfFileFlag)
{
  *endOfFileFlag = 0;
  do
  {
    if(fgets(nextLinePtr, maxSize, filePtr) == NULL)
    {
      *endOfFileFlag = 1;
      break;
    }
    while ((nextLinePtr[0] == ' ' || nextLinePtr[0] == '\t') ||
           (nextLinePtr[0] == '\n' || nextLinePtr[0] == '\r' ))
    {
      nextLinePtr = (nextLinePtr + 1);
    }
  }
  while ((strncmp("#", nextLinePtr, 1) == 0) || (strlen(nextLinePtr) == 0));
}

/******************************************************************************/
double**  AllocateAndInitialize2DArray(int const extentZero,
                                       int const extentOne)
{
  double** arrayPtr;
  int i, j;
  /* allocate memory and set pointers */
  arrayPtr = malloc(extentZero*sizeof(double*));
  arrayPtr[0] = malloc((extentZero * extentOne)*sizeof(double));
  for (i = 1; i < extentZero; ++i)
  {
    arrayPtr[i] = arrayPtr[i-1] + extentOne;
  }

  /* initialize */
  for (i = 0; i < extentZero; ++i)
  {
    for (j = 0; j < extentOne; ++j)
    {
      arrayPtr[i][j] = 0.0;
    }
  }
  return arrayPtr;
}

/******************************************************************************/
double* AllocateAndInitialize1DArray(int const numberModelSpecies)
{ double* arrayPtr;
  int numberUniqueSpeciesPairs;

  /* allocate memory and set pointers */
  numberUniqueSpeciesPairs = ((numberModelSpecies+1)*numberModelSpecies)/2;
  arrayPtr = calloc(numberUniqueSpeciesPairs, sizeof(double));
  return arrayPtr;
}

/******************************************************************************/
void Deallocate1DArrays(int const n, ...)
{
  int i;
  va_list pointerArgs;

  va_start(pointerArgs, n);
  for (i = 0; i < n; i++)
  {
    free(va_arg(pointerArgs, double*));
  }
  va_end(pointerArgs);
}

/******************************************************************************/
void Deallocate2DArrays(int const n, ...)
{
  double** nextArray;
  int i;
  va_list doublePointerArgs;

  va_start(doublePointerArgs, n);
  for (i = 0; i < n; i++)
  {
    nextArray = va_arg(doublePointerArgs, double**);
    free(nextArray[0]);
    free(nextArray);
  }
  va_end(doublePointerArgs);
}
