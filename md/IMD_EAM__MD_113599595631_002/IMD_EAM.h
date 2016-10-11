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
* Copyright (c) 2012,   Institute for Theoretical and Applied Physics
*      			University of Stuttgart, D-70550 Stuttgart, Germany .
* 			All rights reserved.
*
* Contributors:
*    Daniel Schopf
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "KIM_API_C.h"
#include "KIM_API_status.h"

/* dimensions */
#define DIM 3

/* macro definitions */
#define PSTEP 50
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define POS_TRUNC(x) ((int)(x))
#define PTR_2D(var,i,j,dim_i,dim_j) (((var) + ((i)*(dim_j)) + (j)))

/* macro "functions" used by the original IMD code */

/****************************************************************
 *
 *  Evaluate potential table with quadratic interpolation.
 *  Returns the potential value and twice the derivative.
 *  Note: we need (1/r)(dV/dr) = 2 * dV/dr^2 --> use with equidistant r^2
 *  col is p_typ * ntypes + q_typ
 *
 ****************************************************************/

#define PAIR_INT(pot, grad, pt, col, inc, r2, is_short)                      \
{                                                                            \
  double r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;                         \
  int  k;                                                                    \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2a * istep;                                                       \
  k     = POS_TRUNC(r2a);                                                    \
  chi   = r2a - k;                                                           \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k, (col), (pt).maxsteps, (inc));                  \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr;                                                                \
  dv  = p1 - p0;                                                             \
  d2v = p2 - 2 * p1 + p0;                                                    \
                                                                             \
  /* potential and twice the derivative */                                   \
  pot  = p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;                        \
  grad = 2 * istep * (dv + (chi - 0.5) * d2v);                               \
}

/****************************************************************
 *
 *  Evaluate potential table with quadratic interpolation.
 *  Returns the potential value and twice the derivative.
 *  Note: we need (1/r)(dV/dr) = 2 * dV/dr^2 --> use with equidistant r^2
 *  col is p_typ * ntypes + q_typ
 *
 ****************************************************************/

#define PAIR_INT2(pot, grad, pt, col, inc, r2)                      	     \
{                                                                            \
  double r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;                         \
  int  k;                                                                    \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2a * istep;                                                       \
  k     = POS_TRUNC(r2a);                                                    \
  chi   = r2a - k;                                                           \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k, (col), (pt).maxsteps, (inc));                  \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr;                                                                \
  dv  = p1 - p0;                                                             \
  d2v = p2 - 2 * p1 + p0;                                                    \
                                                                             \
  /* potential and twice the derivative */                                   \
  pot  = p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;                        \
  grad = 2 * istep * (dv + (chi - 0.5) * d2v);                               \
}

/****************************************************************
 *
 *  Evaluate tabulated function with quadratic interpolation.
 *
 ****************************************************************/

#define VAL_FUNC(val, pt, col, inc, r2, is_short)                            \
{                                                                            \
  double r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;                         \
  int  k;                                                                    \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2a * istep;                                                       \
  k     = POS_TRUNC(r2a);                                                    \
  chi   = r2a - k;                                                           \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k, (col), (pt).maxsteps, (inc));                  \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr;                                                                \
  dv  = p1 - p0;                                                             \
  d2v = p2 - 2 * p1 + p0;                                                    \
                                                                             \
  /* the function value */                                                   \
  val = p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;                         \
}

/****************************************************************
 *
 *  Evaluate the derivative of a function with quadratic interpolation.
 *  Returns *twice* the derivative.
 *
 ****************************************************************/

#define DERIV_FUNC(grad, pt, col, inc, r2, is_short)                         \
{                                                                            \
  double r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;                         \
  int  k;                                                                    \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a   = 0;                                                               \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2a * istep;                                                       \
  k     = POS_TRUNC(r2a);                                                    \
  /* k     = MIN( POS_TRUNC(r2a), (pt).len[col]-3 ); */                      \
  chi   = r2a - k;                                                           \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k, (col), (pt).maxsteps, (inc));                  \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr;                                                                \
  dv  = p1 - p0;                                                             \
  d2v = p2 - 2 * p1 + p0;                                                    \
                                                                             \
  /* twice the derivative */                                                 \
  grad = 2 * istep * (dv + (chi - 0.5) * d2v);                               \
}

/* a potential table containing all relevant information (from IMD code) */
typedef struct {
  double *begin;		/* first value in the table */
  double *end;			/* last value in the table (followed by extra zeros) */
  double *step;			/* table increment */
  double *invstep;		/* inverse of increment */
  int  *len;			/* length of the individual columns */
  int   ncols;			/* number of columns in the table */
  int   maxsteps;		/* physical length of the table */
  double *table;		/* the actual data */
} pot_table_t;

/* model buffer for storing additional data within the KIM API */
typedef struct {
  /* indices for data in the KIM API object */
  int   coordinates_ind;
  int   energy_ind;
  int   forces_ind;
  int   get_neigh_ind;
  int   numberContributingParticles_ind;
  int   numberOfParticles_ind;
  int   numberOfSpecies_ind;
  int   particleEnergy_ind;
  int   particleSpecies_ind;
  int   process_dEdr_ind;

  /* additional data which is needed frequently */
  int   model_index_shift;
  int   ntypes;

  /* two arrays and their length, storing the density and the derivative of the embedding function */
  int   table_len;
  double *rho_val;
  double *dF_val;

  /* function pointer to the correct compute routine */
  int  (*compute) (void *);

  /* the IMD potential tables */
  pot_table_t pair_pot;
  pot_table_t transfer_pot;
  pot_table_t embed_pot;
} model_buffer;


/* function prototypes for initialization and computes */

/* model driver initialization routine */
int   model_driver_init(void *, char *, int *, int *);

/* general compute function */
int  compute(void *);

/* compute function for different kind of neighbor lists */
int  compute_cluster(void *);
/* pure lists full */
int  compute_pure_full_loca(void *);
int  compute_pure_full_iter(void *);
/* pure lists half */
int  compute_pure_half_loca(void *);
int  compute_pure_half_iter(void *);
/* lists with rvec full */
int  compute_rvec_full_loca(void *);
int  compute_rvec_full_iter(void *);
/* lists with rvec half */
int  compute_rvec_half_loca(void *);
int  compute_rvec_half_iter(void *);

/* destroy function */
int  destroy(void *);

/* reinit function */
int  reinit(void *);

/* original IMD functions for reading the potential files */
void  read_pot_table(pot_table_t *, char *, int, int, int);
void  read_pot_table1(pot_table_t *, int, int, char *, FILE *, int);
void  read_pot_table2(pot_table_t *, int, int, char *, FILE *, int);
void  init_threepoint(pot_table_t *, int);

/* warning for short distances */
void  short_dist_warning(int, int, int, int *, double *, double);
