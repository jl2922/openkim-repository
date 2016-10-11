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
* Copyright (c) 2013,   Institute for Theoretical and Applied Physics
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

#include "IMD_EAM.h"

/* functions implementing the IMD EAM force routines */

/****************************************************************
 *
 * General compute entry function for the API to call,
 * does not compute anything but selects the correct
 * function according to the neighbor list and calls it
 *
 ****************************************************************/

int compute(void *km)
{
  intptr_t *pkim = *((intptr_t **) km);
  int ier;

  /* pointer to model buffer */
  model_buffer *buffer;

  /* get buffer from KIM object */
  buffer = (model_buffer *) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* call compute routine from buffer */
  return buffer->compute(km);
}

/****************************************************************
 *
 * compute function for neighbor list = CLUSTER
 *
 ****************************************************************/

int compute_cluster(void *km)
{
  intptr_t *pkim = *((intptr_t **) km);

  /* general variables */
  int ier;
  int   i, j, l;

  /* pointer to model buffer */
  model_buffer *buffer;

  /* flags for available computes */
  int   comp_energy;
  int   comp_force;
  int   comp_particleEnergy;
  int   comp_process_dEdr;

  /* pointers to objects in the KIM API */
  int  *nAtoms;
  int  *particleSpecies;
  double *coords;
  double *energy;
  double *force;
  double *particleEnergy;

  /* variables / pointers for objects in the buffer */
  int   ntypes;
  double *dF_val;
  double *rho_val;

  /* variables for dealing with the KIM API */
  double Rij[DIM];
  double *pRij = &(Rij[0]);

  /* variables used for calculations */
  int   col1, col2;
  int   inc;
  int   is_short = 0;
  int   it, jt;
  double r2;
  double rho = 0.;
  double rho_i_prime, rho_j_prime;
  double R;
  double phi;
  double dphi;

  /* get buffer from KIM object */
  buffer = (model_buffer *) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* check to see if we have been asked to compute the energy, forces, and particleEnergies */
  /* *INDENT-OFF* */
  KIM_API_getm_compute_by_index(pkim, &ier, 4 * 3,
		buffer->energy_ind, 		&comp_energy, 		1,
		buffer->forces_ind, 		&comp_force, 		1,
		buffer->particleEnergy_ind, 	&comp_particleEnergy, 	1,
		buffer->process_dEdr_ind, 	&comp_process_dEdr,   	1);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute_by_index", ier);
    return ier;
  }

  /* pull out information by index */
  /* *INDENT-OFF* */
  KIM_API_getm_data_by_index(pkim, &ier, 6 * 3,
		buffer->numberOfParticles_ind, 		&nAtoms,         	1,
		buffer->particleSpecies_ind, 		&particleSpecies,  	1,
		buffer->coordinates_ind, 		&coords,         	1,
		buffer->energy_ind,                     &energy,         	comp_energy,
		buffer->forces_ind,                     &force,          	comp_force,
		buffer->particleEnergy_ind,             &particleEnergy, 	comp_particleEnergy);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data_by_index", ier);
    return ier;
  }

  ntypes = buffer->ntypes;
  inc = ntypes * ntypes;

  /* Check to be sure that the atom types are correct */
  for (i = 0; i < *nAtoms; ++i)
    if (particleSpecies[i] > ntypes) {
      KIM_API_report_error(__LINE__, __FILE__, "Unexpected species type detected", KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }

  /* check if the arrays for rho and dF are big enough */
  if (buffer->table_len < *nAtoms) {
    buffer->rho_val = (double *)realloc(buffer->rho_val, (*nAtoms) * sizeof(double));
    buffer->dF_val = (double *)realloc(buffer->dF_val, (*nAtoms) * sizeof(double));
    if (NULL == buffer->rho_val || NULL == buffer->dF_val) {
      KIM_API_report_error(__LINE__, __FILE__, "realloc", KIM_STATUS_FAIL);
      return EXIT_FAILURE;
    }
    buffer->table_len = *nAtoms;
  }
  rho_val = buffer->rho_val;
  dF_val = buffer->dF_val;

  /* reset arrays */
  for (i = 0; i < *nAtoms; i++) {
    rho_val[i] = 0.0;
    dF_val[i] = 0.0;
  }
  if (comp_energy)
    *energy = 0.0;
  if (comp_force)
    for (i = 0; i < DIM * (*nAtoms); i++)
      force[i] = 0.0;
  if (comp_particleEnergy)
    for (i = 0; i < *nAtoms; ++i)
      particleEnergy[i] = 0.0;

  /* first loop over atoms, calculates pair energies and forces as well as eam energies */
  for (i = 0; i < *nAtoms; i++) {
    it = particleSpecies[i];

    /* loop over all remaining atoms, like in a half neighbor list */
    for (j = i + 1; j < *nAtoms; j++) {

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = it * ntypes + jt;
      col2 = jt * ntypes + it;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = coords[j * DIM + l] - coords[i * DIM + l];
	r2 += Rij[l] * Rij[l];
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if (r2 <= buffer->pair_pot.end[col1]) {

	/* calculate pair potential and the derivative at r2 */
	PAIR_INT(phi, dphi, buffer->pair_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(0, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_energy)
	  *energy += phi;

	if (comp_particleEnergy) {
	  particleEnergy[i] += 0.5 * phi;
	  particleEnergy[j] += 0.5 * phi;
	}

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &j);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }

      /* calculate contribution to density */
      if (r2 < buffer->transfer_pot.end[col1]) {
	VAL_FUNC(rho, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(1, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	rho_val[i] += rho;
      }
      /* assign also to neighbor (half neighbor list) */
      if (it == jt) {
	if (r2 < buffer->transfer_pot.end[col1])
	  rho_val[j] += rho;
      } else {
	if (r2 < buffer->transfer_pot.end[col2]) {
	  VAL_FUNC(rho, buffer->transfer_pot, col2, inc, r2, is_short);
	  if (is_short) {
	    short_dist_warning(1, i, j, particleSpecies, Rij, R);
	    is_short = 0;
	  }
	  rho_val[j] += rho;
	}
      }
    }				/* loop over neighbors */

    /* calculate the embedding energies for all atoms */
    PAIR_INT2(phi, dF_val[i], buffer->embed_pot, particleSpecies[i], ntypes, rho_val[i]);

    if (comp_energy)
      *energy += phi;

    if (comp_particleEnergy)
      particleEnergy[i] += phi;
  }				/* loop over atoms */

  /* second loop over atoms, calculates eam forces */
  for (i = 0; i < *nAtoms; i++) {
    it = particleSpecies[i];

    /* loop over all remaining atoms, like in a half neighbor list */
    for (j = i + 1; j < *nAtoms; j++) {

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = jt * ntypes + it;
      col2 = it * ntypes + jt;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = coords[j * DIM + l] - coords[i * DIM + l];
	r2 += Rij[l] * Rij[l];
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if ((r2 < buffer->transfer_pot.end[col1]) || (r2 < buffer->transfer_pot.end[col2])) {

	/* caculate derivative of the density function for both atoms */
	DERIV_FUNC(rho_i_prime, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(2, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	if (col1 == col2) {
	  rho_j_prime = rho_i_prime;
	} else {
	  DERIV_FUNC(rho_j_prime, buffer->transfer_pot, col2, inc, r2, is_short);
	  if (is_short) {
	    short_dist_warning(2, i, j, particleSpecies, Rij, R);
	    is_short = 0;
	  }
	}

	/* combine all contributions, dF_val[] is too big by a factor of 2 */
	dphi = 0.5 * (dF_val[i] * rho_j_prime + dF_val[j] * rho_i_prime);

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &j);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }
    }				/* loop over neighbors */
  }				/* loop over atoms */

  return KIM_STATUS_OK;
}

/****************************************************************
 *
 * compute function for neighbor list = NEIGH_PURE_F
 * and locator mode
 *
 ****************************************************************/

int compute_pure_full_loca(void *km)
{
  intptr_t *pkim = *((intptr_t **) km);

  /* general variables */
  int ier;
  int   i, j, jj, l;

  /* pointer to model buffer */
  model_buffer *buffer;

  /* flags for available computes */
  int   comp_energy;
  int   comp_force;
  int   comp_particleEnergy;
  int   comp_process_dEdr;

  /* pointers to objects in the KIM API */
  typedef int (*get_neigh_ptr) (void *, int *, int *, int *, int *, int **, double **);
  get_neigh_ptr get_neigh;
  int  *nAtoms;
  int  *particleSpecies;
  double *coords;
  double *energy;
  double *force;
  double *particleEnergy;

  /* variables / pointers for objects in the buffer */
  int   model_index_shift;
  int   ntypes;
  double *dF_val;
  double *rho_val;

  /* variables for dealing with the KIM API */
  int   currentAtom;
  int  *neighListOfCurrentAtom;
  int   numOfAtomNeigh;
  int   one = 1;
  int   request;
  double Rij[DIM];
  double *pRij = &(Rij[0]);
  double *Rij_list;

  /* variables used for calculations */
  int   col1, col2;
  int   inc;
  int   is_short = 0;
  int   it, jt;
  double r2;
  double rho = 0.;
  double rho_i_prime, rho_j_prime;
  double R;
  double phi;
  double dphi;

  /* get buffer from KIM object */
  buffer = (model_buffer *) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* unpack info from the buffer */
  model_index_shift = buffer->model_index_shift;

  /* check to see if we have been asked to compute the energy, forces, and particleEnergies */
  /* *INDENT-OFF* */
  KIM_API_getm_compute_by_index(pkim, &ier, 4 * 3,
		buffer->energy_ind,    		&comp_energy, 		1,
		buffer->forces_ind, 		&comp_force, 		1,
		buffer->particleEnergy_ind, 	&comp_particleEnergy, 	1,
		buffer->process_dEdr_ind, 	&comp_process_dEdr, 	1);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute_by_index", ier);
    return ier;
  }

  /* pull out information by index */
  /* *INDENT-OFF* */
  KIM_API_getm_data_by_index(pkim, &ier, 6 * 3,
		buffer->coordinates_ind, 	&coords, 		1,
		buffer->numberOfParticles_ind,  &nAtoms, 		1,
		buffer->particleSpecies_ind, 	&particleSpecies, 	1,
		buffer->energy_ind, 		&energy, 		comp_energy,
		buffer->forces_ind, 		&force, 		comp_force,
		buffer->particleEnergy_ind, 	&particleEnergy, 	comp_particleEnergy);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data_by_index", ier);
    return ier;
  }
  get_neigh = (get_neigh_ptr) KIM_API_get_method_by_index(pkim, buffer->get_neigh_ind, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_method_by_index", ier);
    return ier;
  }

  ntypes = buffer->ntypes;
  inc = ntypes * ntypes;

  /* Check to be sure that the atom types are correct */
  for (i = 0; i < *nAtoms; ++i)
    if (particleSpecies[i] > ntypes) {
      KIM_API_report_error(__LINE__, __FILE__, "Unexpected species type detected", KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }

  /* check if the arrays for rho and dF are big enough */
  if (buffer->table_len < *nAtoms) {
    buffer->rho_val = (double *)realloc(buffer->rho_val, (*nAtoms) * sizeof(double));
    buffer->dF_val = (double *)realloc(buffer->dF_val, (*nAtoms) * sizeof(double));
    if (NULL == buffer->rho_val || NULL == buffer->dF_val) {
      KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }
    buffer->table_len = *nAtoms;
  }
  rho_val = buffer->rho_val;
  dF_val = buffer->dF_val;

  /* reset arrays */
  for (i = 0; i < *nAtoms; i++) {
    rho_val[i] = 0.0;
    dF_val[i] = 0.0;
  }
  if (comp_energy)
    *energy = 0.0;
  if (comp_force)
    for (i = 0; i < DIM * (*nAtoms); i++)
      force[i] = 0.0;
  if (comp_particleEnergy)
    for (i = 0; i < *nAtoms; ++i)
      particleEnergy[i] = 0.0;

  /* first loop over atoms, calculates pair energies and forces as well as eam energies */
  for (i = 0; i < *nAtoms; i++) {
    request = i - model_index_shift;
    ier = (*get_neigh) (km, &one, &request, &currentAtom,
      &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
    if (KIM_STATUS_OK > ier) {
      KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
      return ier;
    }
    it = particleSpecies[i];

    /* loop over all atoms in the neighbor list */
    for (jj = 0; jj < numOfAtomNeigh; jj++) {
      j = neighListOfCurrentAtom[jj] + model_index_shift;

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = it * ntypes + jt;
      col2 = jt * ntypes + it;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = coords[j * DIM + l] - coords[i * DIM + l];
	r2 += Rij[l] * Rij[l];
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if (r2 <= buffer->pair_pot.end[col1]) {

	/* calculate pair potential and the derivative at r2 */
	PAIR_INT(phi, dphi, buffer->pair_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(0, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	/* we have a full neighbor list, so we only contribute 0.5 */
	phi *= 0.5;
	dphi *= 0.5;

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_energy)
	  *energy += phi;

	if (comp_particleEnergy)
	  particleEnergy[i] += phi;

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &j);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }

      /* calculate contribution to density */
      if (r2 < buffer->transfer_pot.end[col1]) {
	VAL_FUNC(rho, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(1, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	rho_val[i] += rho;
      }
      /* do NOT assign to neighbor (full neighbor list) */
    }				/* loop over neighbors */

    /* calculate the embedding energies for all contributing atoms */
    if (numOfAtomNeigh > 0) {
      PAIR_INT2(phi, dF_val[i], buffer->embed_pot, particleSpecies[i], ntypes, rho_val[i]);

      if (comp_energy)
	*energy += phi;

      if (comp_particleEnergy)
	particleEnergy[i] += phi;
    }
  }				/* loop over atoms */

  /* second loop over atoms, calculates eam forces */
  for (i = 0; i < *nAtoms; i++) {
    request = i - model_index_shift;
    ier = (*get_neigh) (km, &one, &request, &currentAtom,
      &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
    if (KIM_STATUS_OK > ier) {
      KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
      return ier;
    }
    it = particleSpecies[i];

    /* loop over all atoms in the neighbor list */
    for (jj = 0; jj < numOfAtomNeigh; jj++) {
      j = neighListOfCurrentAtom[jj] + model_index_shift;

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = jt * ntypes + it;
      col2 = it * ntypes + jt;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = coords[j * DIM + l] - coords[i * DIM + l];
	r2 += (Rij[l] * Rij[l]);
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if ((r2 < buffer->transfer_pot.end[col1]) || (r2 < buffer->transfer_pot.end[col2])) {

	/* caculate derivative of the density function for both atoms */
	DERIV_FUNC(rho_i_prime, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(2, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	if (col1 == col2) {
	  rho_j_prime = rho_i_prime;
	} else {
	  DERIV_FUNC(rho_j_prime, buffer->transfer_pot, col2, inc, r2, is_short);
	  if (is_short) {
	    short_dist_warning(2, i, j, particleSpecies, Rij, R);
	    is_short = 0;
	  }
	}

	/* combine all contributions, dF_val[] is too big by a factor of 2 */
	dphi = 0.5 * dF_val[i] * rho_j_prime;

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &j);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }
    }				/* loop over neighbors */
  }				/* loop over atoms */

  return KIM_STATUS_OK;
}

/****************************************************************
 *
 * compute function for neighbor list = NEIGH_PURE_F
 * and iterator mode
 *
 ****************************************************************/

int compute_pure_full_iter(void *km)
{
  intptr_t *pkim = *((intptr_t **) km);

  /* general variables */
  int ier;
  int   i, j, jj, l;

  /* pointer to model buffer */
  model_buffer *buffer;

  /* flags for available computes */
  int   comp_energy;
  int   comp_force;
  int   comp_particleEnergy;
  int   comp_process_dEdr;

  /* pointers to objects in the KIM API */
  typedef int   (*get_neigh_ptr) (void *, int *, int *, int *, int *, int **, double **);
  get_neigh_ptr get_neigh;
  int  *nAtoms;
  int  *particleSpecies;
  double *coords;
  double *energy;
  double *force;
  double *particleEnergy;

  /* variables / pointers for objects in the buffer */
  int   model_index_shift;
  int   ntypes;
  double *dF_val;
  double *rho_val;

  /* variables for dealing with the KIM API */
  int   currentAtom;
  int  *neighListOfCurrentAtom;
  int   numOfAtomNeigh;
  int   one = 1;
  int   zero = 0;
  double Rij[DIM];
  double *pRij = &(Rij[0]);
  double *Rij_list;

  /* variables used for calculations */
  int   col1, col2;
  int   inc;
  int   is_short = 0;
  int   it, jt;
  double r2;
  double rho = 0.0;
  double rho_i_prime, rho_j_prime;
  double R;
  double phi;
  double dphi;

  /* get buffer from KIM object */
  buffer = (model_buffer *) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* unpack info from the buffer */
  model_index_shift = buffer->model_index_shift;

  /* check to see if we have been asked to compute the energy, forces, and particleEnergies */
  /* *INDENT-OFF* */
  KIM_API_getm_compute_by_index(pkim, &ier, 4 * 3,
		buffer->energy_ind, 		&comp_energy, 		1,
		buffer->forces_ind, 		&comp_force, 		1,
		buffer->particleEnergy_ind, 	&comp_particleEnergy, 	1,
		buffer->process_dEdr_ind, 	&comp_process_dEdr, 	1);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute_by_index", ier);
    return ier;
  }

  /* pull out information by index */
  /* *INDENT-OFF* */
  KIM_API_getm_data_by_index(pkim, &ier, 6 * 3,
		buffer->coordinates_ind, 	&coords,         	1,
		buffer->numberOfParticles_ind, 	&nAtoms,         	1,
		buffer->particleSpecies_ind, 	&particleSpecies,  	1,
		buffer->energy_ind, 		&energy,         	comp_energy,
		buffer->forces_ind, 		&force,          	comp_force,
		buffer->particleEnergy_ind, 	&particleEnergy, 	comp_particleEnergy);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data_by_index", ier);
    return ier;
  }
  get_neigh = (get_neigh_ptr) KIM_API_get_method_by_index(pkim,buffer->get_neigh_ind, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_method_by_index", ier);
    return ier;
  }

  ntypes = buffer->ntypes;
  inc = ntypes * ntypes;

  /* Check to be sure that the atom types are correct */
  for (i = 0; i < *nAtoms; ++i)
    if (particleSpecies[i] > ntypes) {
      KIM_API_report_error(__LINE__, __FILE__, "Unexpected species type detected", KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }

  /* check if the arrays for rho and dF are big enough */
  if (buffer->table_len < *nAtoms) {
    buffer->rho_val = (double *)realloc(buffer->rho_val, (*nAtoms) * sizeof(double));
    buffer->dF_val = (double *)realloc(buffer->dF_val, (*nAtoms) * sizeof(double));
    if (NULL == buffer->rho_val || NULL == buffer->dF_val) {
      KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }
    buffer->table_len = *nAtoms;
  }
  rho_val = buffer->rho_val;
  dF_val = buffer->dF_val;

  /* reset arrays */
  for (i = 0; i < *nAtoms; i++) {
    rho_val[i] = 0.0;
    dF_val[i] = 0.0;
  }
  if (comp_energy)
    *energy = 0.0;
  if (comp_force)
    for (i = 0; i < DIM * (*nAtoms); i++)
      force[i] = 0.0;
  if (comp_particleEnergy)
    for (i = 0; i < *nAtoms; ++i)
      particleEnergy[i] = 0.0;

  /* initialize the iterator */
  ier = (*get_neigh) (km, &zero, &zero, &currentAtom, &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
  if (KIM_STATUS_NEIGH_ITER_INIT_OK != ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", ier);
    ier = KIM_STATUS_FAIL;
    return ier;
  }

  i = -1;
  /* first infinite iterator loop, calculates pair energies and forces as well as eam energies */
  while (1) {

    /* get the next neighbor */
    ier = (*get_neigh) (km, &zero, &one, &currentAtom, &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
    /* abort the infinite loop if we are past the end of the list */
    if (KIM_STATUS_NEIGH_ITER_PAST_END == ier)
      break;
    if (KIM_STATUS_OK != ier) {
      KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
      ier = KIM_STATUS_FAIL;
      return ier;
    }
    i = currentAtom + model_index_shift;
    it = particleSpecies[i];

    /* loop over all atoms in the neighbor list */
    for (jj = 0; jj < numOfAtomNeigh; jj++) {
      j = neighListOfCurrentAtom[jj] + model_index_shift;

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = it * ntypes + jt;
      col2 = jt * ntypes + it;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = coords[j * DIM + l] - coords[i * DIM + l];
	r2 += Rij[l] * Rij[l];
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if (r2 <= buffer->pair_pot.end[col1]) {

	/* calculate pair potential and the derivative at r2 */
	PAIR_INT(phi, dphi, buffer->pair_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(0, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	/* we have a full neighbor list, so we only contribute 0.5 */
	phi *= 0.5;
	dphi *= 0.5;

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_energy)
	  *energy += phi;

	if (comp_particleEnergy)
	  particleEnergy[i] += phi;

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &j);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }

      /* calculate contribution to density */
      if (r2 < buffer->transfer_pot.end[col1]) {
	VAL_FUNC(rho, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(1, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	rho_val[i] += rho;
      }
      /* do NOT assign to neighbor (full neighbor list) */
    }				/* loop over neighbors */

    /* calculate the embedding energies for all contributing atoms */
    if (numOfAtomNeigh > 0) {
      PAIR_INT2(phi, dF_val[i], buffer->embed_pot, particleSpecies[i], ntypes, rho_val[i]);

      if (comp_energy)
	*energy += phi;

      if (comp_particleEnergy)
	particleEnergy[i] += phi;
    }
  }				/* infinite loop */

  /* initialize the iterator */
  ier = (*get_neigh) (km, &zero, &zero, &currentAtom, &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
  if (KIM_STATUS_NEIGH_ITER_INIT_OK != ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", ier);
    return KIM_STATUS_FAIL;
  }

  i = -1;
  /* second infinite iterator loop, calculates eam forces */
  while (1) {

    /* get the next neighbor */
    ier = (*get_neigh) (km, &zero, &one, &currentAtom, &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
    /* abort the infinite loop if we are past the end of the list */
    if (KIM_STATUS_NEIGH_ITER_PAST_END == ier)
      break;
    if (KIM_STATUS_OK > ier) {
      KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
      return KIM_STATUS_FAIL;
    }
    i = currentAtom + model_index_shift;
    it = particleSpecies[i];

    /* loop over all atoms in the neighbor list */
    for (jj = 0; jj < numOfAtomNeigh; jj++) {
      j = neighListOfCurrentAtom[jj] + model_index_shift;

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = jt * ntypes + it;
      col2 = it * ntypes + jt;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = coords[j * DIM + l] - coords[i * DIM + l];
	r2 += Rij[l] * Rij[l];
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if ((r2 < buffer->transfer_pot.end[col1]) || (r2 < buffer->transfer_pot.end[col2])) {

	/* caculate derivative of the density function for both atoms */
	DERIV_FUNC(rho_i_prime, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(2, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	if (col1 == col2) {
	  rho_j_prime = rho_i_prime;
	} else {
	  DERIV_FUNC(rho_j_prime, buffer->transfer_pot, col2, inc, r2, is_short);
	  if (is_short) {
	    short_dist_warning(2, i, j, particleSpecies, Rij, R);
	    is_short = 0;
	  }
	}

	/* combine all contributions, dF_val[] is too big by a factor of 2 */
	dphi = 0.5 * dF_val[i] * rho_j_prime;

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &j);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }
    }				/* loop over neighbors */
  }				/* infinite loop */

  return KIM_STATUS_OK;
}

/****************************************************************
 *
 * compute function for neighbor list = NEIGH_PURE_H
 * and locator mode
 *
 ****************************************************************/

int compute_pure_half_loca(void *km)
{
  intptr_t *pkim = *((intptr_t **) km);

  /* general variables */
  int ier;
  int   i, j, jj, l;

  /* pointer to model buffer */
  model_buffer *buffer;

  /* flags for available computes */
  int   comp_energy;
  int   comp_force;
  int   comp_particleEnergy;
  int   comp_process_dEdr;

  /* pointers to objects in the KIM API */
  typedef int   (*get_neigh_ptr) (void *, int *, int *, int *, int *, int **, double **);
  get_neigh_ptr get_neigh;
  int  *nAtoms;
  int  *numContrib;
  int  *particleSpecies;
  double *coords;
  double *energy;
  double *force;
  double *particleEnergy;

  /* variables / pointers for objects in the buffer */
  int   model_index_shift;
  int   ntypes;
  double *dF_val;
  double *rho_val;

  /* variables for dealing with the KIM API */
  int   currentAtom;
  int  *neighListOfCurrentAtom;
  int   numOfAtomNeigh;
  int   one = 1;
  int   request;
  double Rij[DIM];
  double *pRij = &(Rij[0]);
  double *Rij_list;

  /* variables used for calculations */
  int   col1, col2;
  int   inc;
  int   is_short = 0;
  int   it, jt;
  int   numberContrib;
  double r2;
  double rho = 0.0;
  double rho_i_prime, rho_j_prime;
  double R;
  double phi;
  double dphi;

  /* get buffer from KIM object */
  buffer = (model_buffer *) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* unpack info from the buffer */
  model_index_shift = buffer->model_index_shift;

  /* check to see if we have been asked to compute the energy, forces, and particleEnergies */
  /* *INDENT-OFF* */
  KIM_API_getm_compute_by_index(pkim, &ier, 4 * 3,
		buffer->energy_ind, 		&comp_energy, 		1,
		buffer->forces_ind, 		&comp_force, 		1,
		buffer->particleEnergy_ind, 	&comp_particleEnergy, 	1,
		buffer->process_dEdr_ind, 	&comp_process_dEdr, 	1);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute_by_index", ier);
    return ier;
  }

  /* pull out information by index */
  /* *INDENT-OFF* */
  KIM_API_getm_data_by_index(pkim, &ier, 7 * 3,
		buffer->coordinates_ind, 			&coords, 		1,
		buffer->numberOfParticles_ind, 			&nAtoms, 		1,
		buffer->numberContributingParticles_ind, 	&numContrib, 		1,
		buffer->particleSpecies_ind, 			&particleSpecies, 	1,
		buffer->energy_ind, 				&energy,    		comp_energy,
		buffer->forces_ind, 				&force, 		comp_force,
    		buffer->particleEnergy_ind, 			&particleEnergy, 	comp_particleEnergy);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data_by_index", ier);
    return ier;
  }
  get_neigh = (get_neigh_ptr) KIM_API_get_method_by_index(pkim, buffer->get_neigh_ind, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_method_by_index", ier);
    return ier;
  }

  numberContrib = *numContrib;
  ntypes = buffer->ntypes;
  inc = ntypes * ntypes;

  /* Check to be sure that the atom types are correct */
  for (i = 0; i < *nAtoms; ++i)
    if (particleSpecies[i] > ntypes) {
      KIM_API_report_error(__LINE__, __FILE__, "Unexpected species type detected", KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }

  /* check if the arrays for rho and dF are big enough */
  if (buffer->table_len < *nAtoms) {
    buffer->rho_val = (double *)realloc(buffer->rho_val, (*nAtoms) * sizeof(double));
    buffer->dF_val = (double *)realloc(buffer->dF_val, (*nAtoms) * sizeof(double));
    if (NULL == buffer->rho_val || NULL == buffer->dF_val) {
      KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }
    buffer->table_len = *nAtoms;
  }
  rho_val = buffer->rho_val;
  dF_val = buffer->dF_val;

  /* reset arrays */
  for (i = 0; i < *nAtoms; i++) {
    rho_val[i] = 0.0;
    dF_val[i] = 0.0;
  }
  if (comp_energy)
    *energy = 0.0;
  if (comp_force)
    for (i = 0; i < DIM * (*nAtoms); i++)
      force[i] = 0.0;
  if (comp_particleEnergy)
    for (i = 0; i < *nAtoms; ++i)
      particleEnergy[i] = 0.0;

  /* first loop over contributing atoms, calculates pair energies and forces as well as eam energies */
  for (i = 0; i < numberContrib; i++) {
    request = i - model_index_shift;
    ier = (*get_neigh) (km, &one, &request, &currentAtom,
      &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
    if (KIM_STATUS_OK != ier) {
      KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
      return KIM_STATUS_FAIL;
    }
    it = particleSpecies[i];

    /* loop over all atoms in the neighbor list */
    for (jj = 0; jj < numOfAtomNeigh; jj++) {
      j = neighListOfCurrentAtom[jj] + model_index_shift;

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = it * ntypes + jt;
      col2 = jt * ntypes + it;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = coords[j * DIM + l] - coords[i * DIM + l];
	r2 += Rij[l] * Rij[l];
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if (r2 <= buffer->pair_pot.end[col1]) {

	/* calculate pair potential and the derivative at r2 */
	PAIR_INT(phi, dphi, buffer->pair_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(0, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	phi *= 0.5;
	if (j >= numberContrib)
	  dphi *= 0.5;

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_energy) {
	  *energy += phi;
	  if (j < numberContrib)
	    *energy += phi;
	}

	if (comp_particleEnergy) {
	  particleEnergy[i] += phi;
	  particleEnergy[j] += phi;
	}

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &jj);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }

      /* calculate contribution to density */
      if (r2 < buffer->transfer_pot.end[col1]) {
	VAL_FUNC(rho, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(1, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	rho_val[i] += rho;
      }
      /* assign also to neighbor (half neighbor list) */
      if (it == jt) {
	if (r2 < buffer->transfer_pot.end[col2])
	  rho_val[j] += rho;
      } else {
	if (r2 < buffer->transfer_pot.end[col2]) {
	  VAL_FUNC(rho, buffer->transfer_pot, col2, inc, r2, is_short);
	  if (is_short) {
	    short_dist_warning(1, i, j, particleSpecies, Rij, R);
	    is_short = 0;
	  }
	  rho_val[j] += rho;
	}
      }
    }				/* loop over neighbors */

    /* calculate the embedding energies for all contributing atoms */
    if (i < numberContrib) {
      PAIR_INT2(phi, dF_val[i], buffer->embed_pot, particleSpecies[i], ntypes, rho_val[i]);

      if (comp_energy)
	*energy += phi;

      if (comp_particleEnergy)
	particleEnergy[i] += phi;
    }
  }				/* loop over contributing atoms */

  /* second loop over contributing atoms, calculates eam forces */
  for (i = 0; i < numberContrib; i++) {
    request = i - model_index_shift;
    ier = (*get_neigh) (km, &one, &request, &currentAtom,
      &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
    if (KIM_STATUS_OK != ier) {
      KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
      return KIM_STATUS_FAIL;
    }
    it = particleSpecies[i];

    /* loop over all atoms in the neighbor list */
    for (jj = 0; jj < numOfAtomNeigh; jj++) {
      j = neighListOfCurrentAtom[jj] + model_index_shift;

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = jt * ntypes + it;
      col2 = it * ntypes + jt;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = coords[j * DIM + l] - coords[i * DIM + l];
	r2 += Rij[l] * Rij[l];
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if ((r2 < buffer->transfer_pot.end[col1]) || (r2 < buffer->transfer_pot.end[col2])) {

	/* caculate derivative of the density function for both atoms */
	DERIV_FUNC(rho_i_prime, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(2, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	if (col1 == col2) {
	  rho_j_prime = rho_i_prime;
	} else {
	  DERIV_FUNC(rho_j_prime, buffer->transfer_pot, col2, inc, r2, is_short);
	  if (is_short) {
	    short_dist_warning(2, i, j, particleSpecies, Rij, R);
	    is_short = 0;
	  }
	}

	/* combine all contributions, dF_val[] is too big by a factor of 2 */
	dphi = 0.5 * (dF_val[i] * rho_j_prime + dF_val[j] * rho_i_prime);

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &j);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }
    }				/* loop over neighbors */
  }				/* loop over contributing atoms */

  return KIM_STATUS_OK;
}

/****************************************************************
 *
 * compute function for neighbor list = NEIGH_PURE_H
 * and iterator mode
 *
 ****************************************************************/

int compute_pure_half_iter(void *km)
{
  intptr_t *pkim = *((intptr_t **) km);

  /* general variables */
  int ier;
  int   i, j, jj, l;

  /* pointer to model buffer */
  model_buffer *buffer;

  /* flags for available computes */
  int   comp_energy;
  int   comp_force;
  int   comp_particleEnergy;
  int   comp_process_dEdr;

  /* pointers to objects in the KIM API */
  typedef int   (*get_neigh_ptr) (void *, int *, int *, int *, int *, int **, double **);
  get_neigh_ptr get_neigh;
  int  *nAtoms;
  int  *numContrib;
  int  *particleSpecies;
  double *coords;
  double *energy;
  double *force;
  double *particleEnergy;

  /* variables / pointers for objects in the buffer */
  int   model_index_shift;
  int   ntypes;
  double *dF_val;
  double *rho_val;

  /* variables for dealing with the KIM API */
  int   currentAtom;
  int  *neighListOfCurrentAtom;
  int   numOfAtomNeigh;
  int   one = 1;
  int   zero = 0;
  double Rij[DIM];
  double *pRij = &(Rij[0]);
  double *Rij_list;

  /* variables used for calculations */
  int   col1, col2;
  int   inc;
  int   is_short = 0;
  int   it, jt;
  int   numberContrib;
  double r2;
  double rho = 0.0;
  double rho_i_prime, rho_j_prime;
  double R;
  double phi;
  double dphi;

  /* get buffer from KIM object */
  buffer = (model_buffer *) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* unpack info from the buffer */
  model_index_shift = buffer->model_index_shift;

  /* check to see if we have been asked to compute the energy, forces, and particleEnergies */
  /* *INDENT-OFF* */
  KIM_API_getm_compute_by_index(pkim, &ier, 4 * 3,
	  	buffer->energy_ind, 		&comp_energy, 		1,
	  	buffer->forces_ind, 		&comp_force, 		1,
		buffer->particleEnergy_ind, 	&comp_particleEnergy, 	1,
		buffer->process_dEdr_ind, 	&comp_process_dEdr, 	1);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute_by_index", ier);
    return ier;
  }

  /* pull out information by index */
  /* *INDENT-OFF* */
  KIM_API_getm_data_by_index(pkim, &ier, 7 * 3,
		buffer->coordinates_ind,                 	&coords,         	1,
               	buffer->numberOfParticles_ind,           	&nAtoms,         	1,
		buffer->numberContributingParticles_ind, 	&numContrib, 		1,
                buffer->particleSpecies_ind,               	&particleSpecies,  	1,
                buffer->energy_ind,                      	&energy,         	comp_energy,
                buffer->forces_ind,                      	&force,          	comp_force,
                buffer->particleEnergy_ind,              	&particleEnergy, 	comp_particleEnergy);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data_by_index", ier);
    return ier;
  }
  get_neigh = (get_neigh_ptr) KIM_API_get_method_by_index(pkim, buffer->get_neigh_ind, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_method_by_index", ier);
    return ier;
  }

  numberContrib = *numContrib;
  ntypes = buffer->ntypes;
  inc = ntypes * ntypes;

  /* Check to be sure that the atom types are correct */
  for (i = 0; i < *nAtoms; ++i)
    if (particleSpecies[i] > ntypes) {
      KIM_API_report_error(__LINE__, __FILE__, "Unexpected species type detected", KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }

  /* check if the arrays for rho and dF are big enough */
  if (buffer->table_len < *nAtoms) {
    buffer->rho_val = (double *)realloc(buffer->rho_val, (*nAtoms) * sizeof(double));
    buffer->dF_val = (double *)realloc(buffer->dF_val, (*nAtoms) * sizeof(double));
    if (NULL == buffer->rho_val || NULL == buffer->dF_val) {
      KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }
    buffer->table_len = *nAtoms;
  }
  rho_val = buffer->rho_val;
  dF_val = buffer->dF_val;

  /* reset arrays */
  for (i = 0; i < *nAtoms; i++) {
    rho_val[i] = 0.0;
    dF_val[i] = 0.0;
  }
  if (comp_energy)
    *energy = 0.0;
  if (comp_force)
    for (i = 0; i < DIM * (*nAtoms); i++)
      force[i] = 0.0;
  if (comp_particleEnergy)
    for (i = 0; i < *nAtoms; ++i)
      particleEnergy[i] = 0.0;

  /* initialize the iterator */
  ier = (*get_neigh) (km, &zero, &zero, &currentAtom, &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
  if (KIM_STATUS_NEIGH_ITER_INIT_OK != ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", ier);
    return KIM_STATUS_FAIL;
  }

  i = -1;
  /* first infinite iterator loop, calculates pair energies and forces as well as eam energies */
  while (1) {

    /* get the next neighbor */
    ier = (*get_neigh) (km, &zero, &one, &currentAtom, &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
    /* abort the infinite loop if we are past the end of the list */
    if (KIM_STATUS_NEIGH_ITER_PAST_END == ier)
      break;
    if (KIM_STATUS_OK != ier) {
      KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
      return KIM_STATUS_FAIL;
    }
    i = currentAtom + model_index_shift;
    it = particleSpecies[i];

    /* loop over all atoms in the neighbor list */
    for (jj = 0; jj < numOfAtomNeigh; jj++) {
      j = neighListOfCurrentAtom[jj] + model_index_shift;

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = it * ntypes + jt;
      col2 = jt * ntypes + it;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = coords[j * DIM + l] - coords[i * DIM + l];
	r2 += Rij[l] * Rij[l];
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if (r2 <= buffer->pair_pot.end[col1]) {

	/* calculate pair potential and the derivative at r2 */
	PAIR_INT(phi, dphi, buffer->pair_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(0, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	phi *= 0.5;
	if (j >= numberContrib)
	  dphi *= 0.5;

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_energy) {
	  *energy += phi;
	  if (j < numberContrib)
	    *energy += phi;
	}

	if (comp_particleEnergy) {
	  particleEnergy[i] += phi;
	  particleEnergy[j] += phi;
	}

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &j);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }

      /* calculate contribution to density */
      if (r2 < buffer->transfer_pot.end[col1]) {
	VAL_FUNC(rho, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(1, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	rho_val[i] += rho;
      }
      /* assign also to neighbor (half neighbor list) */
      if (it == jt) {
	if (r2 < buffer->transfer_pot.end[col2])
	  rho_val[j] += rho;
      } else {
	if (r2 < buffer->transfer_pot.end[col2]) {
	  VAL_FUNC(rho, buffer->transfer_pot, col2, inc, r2, is_short);
	  if (is_short) {
	    short_dist_warning(1, i, j, particleSpecies, Rij, R);
	    is_short = 0;
	  }
	  rho_val[j] += rho;
	}
      }
    }				/* loop over neighbors */

    /* calculate the embedding energies for all contributing atoms */
    if (i < numberContrib) {
      PAIR_INT2(phi, dF_val[i], buffer->embed_pot, particleSpecies[i], ntypes, rho_val[i]);

      if (comp_energy)
	*energy += phi;

      if (comp_particleEnergy)
	particleEnergy[i] += phi;
    }
  }				/* infinite loop */

  /* initialize the iterator */
  ier = (*get_neigh) (km, &zero, &zero, &currentAtom, &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
  if (KIM_STATUS_NEIGH_ITER_INIT_OK != ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", ier);
    return KIM_STATUS_FAIL;
  }

  i = -1;
  /* second infinite iterator loop, calculates eam forces */
  while (1) {

    /* get the next neighbor */
    ier = (*get_neigh) (km, &zero, &one, &currentAtom, &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
    /* abort the infinite loop if we are past the end of the list */
    if (KIM_STATUS_NEIGH_ITER_PAST_END == ier)
      break;
    if (KIM_STATUS_OK != ier) {
      KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
      return KIM_STATUS_FAIL;
    }
    i = currentAtom + model_index_shift;
    it = particleSpecies[i];

    /* loop over all atoms in the neighbor list */
    for (jj = 0; jj < numOfAtomNeigh; jj++) {
      j = neighListOfCurrentAtom[jj] + model_index_shift;

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = jt * ntypes + it;
      col2 = it * ntypes + jt;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = coords[j * DIM + l] - coords[i * DIM + l];
	r2 += Rij[l] * Rij[l];
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if ((r2 < buffer->transfer_pot.end[col1]) || (r2 < buffer->transfer_pot.end[col2])) {

	/* caculate derivative of the density function for both atoms */
	DERIV_FUNC(rho_i_prime, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(2, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	if (col1 == col2) {
	  rho_j_prime = rho_i_prime;
	} else {
	  DERIV_FUNC(rho_j_prime, buffer->transfer_pot, col2, inc, r2, is_short);
	  if (is_short) {
	    short_dist_warning(2, i, j, particleSpecies, Rij, R);
	    is_short = 0;
	  }
	}

	/* combine all contributions, dF_val[] is too big by a factor of 2 */
	dphi = 0.5 * (dF_val[i] * rho_j_prime + dF_val[j] * rho_i_prime);

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &j);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }
    }				/* loop over neighbors */
  }				/* infinite loop */

  return KIM_STATUS_OK;
}

/****************************************************************
 *
 * compute function for neighbor list = NEIGH_RVEC_F
 * and locator mode
 *
 ****************************************************************/

int compute_rvec_full_loca(void *km)
{
  intptr_t *pkim = *((intptr_t **) km);

  /* general variables */
  int ier;
  int   i, j, jj, l;

  /* pointer to model buffer */
  model_buffer *buffer;

  /* flags for available computes */
  int   comp_energy;
  int   comp_force;
  int   comp_particleEnergy;
  int   comp_process_dEdr;

  /* pointers to objects in the KIM API */
  typedef int   (*get_neigh_ptr) (void *, int *, int *, int *, int *, int **, double **);
  get_neigh_ptr get_neigh;
  int  *nAtoms;
  int  *particleSpecies;
  double *energy;
  double *force;
  double *particleEnergy;

  /* variables / pointers for objects in the buffer */
  int   model_index_shift;
  int   ntypes;
  double *dF_val;
  double *rho_val;

  /* variables for dealing with the KIM API */
  int   currentAtom;
  int  *neighListOfCurrentAtom;
  int   numOfAtomNeigh;
  int   one = 1;
  int   request;
  double Rij[DIM];
  double *pRij = &(Rij[0]);
  double *Rij_list;

  /* variables used for calculations */
  int   col1, col2;
  int   inc;
  int   is_short = 0;
  int   it, jt;
  double r2;
  double rho = 0.0;
  double rho_i_prime, rho_j_prime;
  double R;
  double phi;
  double dphi;

  /* get buffer from KIM object */
  buffer = (model_buffer *) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* unpack info from the buffer */
  model_index_shift = buffer->model_index_shift;

  /* check to see if we have been asked to compute the energy, forces, and particleEnergies */
  /* *INDENT-OFF* */
  KIM_API_getm_compute_by_index(pkim, &ier, 4 * 3,
		  buffer->energy_ind,    	&comp_energy, 		1,
		  buffer->forces_ind, 		&comp_force, 		1,
		  buffer->particleEnergy_ind, 	&comp_particleEnergy, 	1,
		  buffer->process_dEdr_ind, 	&comp_process_dEdr, 	1);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute_by_index", ier);
    return ier;
  }

  /* pull out information by index */
  /* *INDENT-OFF* */
  KIM_API_getm_data_by_index(pkim, &ier, 5 * 3,
		  buffer->numberOfParticles_ind,    		&nAtoms, 		1,
		  buffer->particleSpecies_ind, 			&particleSpecies, 	1,
		  buffer->energy_ind, 				&energy, 		comp_energy,
		  buffer->forces_ind, 				&force, 		comp_force,
		  buffer->particleEnergy_ind, 			&particleEnergy, 	comp_particleEnergy);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data_by_index", ier);
    return ier;
  }
  get_neigh = (get_neigh_ptr) KIM_API_get_method_by_index(pkim, buffer->get_neigh_ind, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_method_by_index", ier);
    return ier;
  }

  ntypes = buffer->ntypes;
  inc = ntypes * ntypes;

  /* Check to be sure that the atom types are correct */
  for (i = 0; i < *nAtoms; ++i)
    if (particleSpecies[i] > ntypes) {
      KIM_API_report_error(__LINE__, __FILE__, "Unexpected species type detected", KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }

  /* check if the arrays for rho and dF are big enough */
  if (buffer->table_len < *nAtoms) {
    buffer->rho_val = (double *)realloc(buffer->rho_val, (*nAtoms) * sizeof(double));
    buffer->dF_val = (double *)realloc(buffer->dF_val, (*nAtoms) * sizeof(double));
    if (NULL == buffer->rho_val || NULL == buffer->dF_val) {
      KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }
    buffer->table_len = *nAtoms;
  }
  rho_val = buffer->rho_val;
  dF_val = buffer->dF_val;

  /* reset arrays */
  for (i = 0; i < *nAtoms; i++) {
    rho_val[i] = 0.0;
    dF_val[i] = 0.0;
  }
  if (comp_energy)
    *energy = 0.0;
  if (comp_force)
    for (i = 0; i < DIM * (*nAtoms); i++)
      force[i] = 0.0;
  if (comp_particleEnergy)
    for (i = 0; i < *nAtoms; ++i)
      particleEnergy[i] = 0.0;

  /* first loop over atoms, calculates pair energies and forces as well as eam energies */
  for (i = 0; i < *nAtoms; i++) {
    request = i - model_index_shift;
    ier = (*get_neigh) (km, &one, &request, &currentAtom,
      &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
    if (KIM_STATUS_OK != ier) {
      KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
      return KIM_STATUS_FAIL;
    }
    it = particleSpecies[i];

    /* loop over all atoms in the neighbor list */
    for (jj = 0; jj < numOfAtomNeigh; jj++) {
      j = neighListOfCurrentAtom[jj] + model_index_shift;

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = it * ntypes + jt;
      col2 = jt * ntypes + it;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = Rij_list[jj * DIM + l];
	r2 += Rij[l] * Rij[l];
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if (r2 <= buffer->pair_pot.end[col1]) {

	/* calculate pair potential and the derivative at r2 */
	PAIR_INT(phi, dphi, buffer->pair_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(0, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	/* we have a full neighbor list, so we only contribute 0.5 */
	phi *= 0.5;
	dphi *= 0.5;

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_energy)
	  *energy += phi;

	if (comp_particleEnergy)
	  particleEnergy[i] += phi;

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &j);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }

      /* calculate contribution to density */
      if (r2 < buffer->transfer_pot.end[col1]) {
	VAL_FUNC(rho, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(1, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	rho_val[i] += rho;
      }
      /* do NOT assign to neighbor (full neighbor list) */
    }				/* loop over neighbors */

    /* calculate the embedding energies for all contributing atoms */
    if (numOfAtomNeigh > 0) {
      PAIR_INT2(phi, dF_val[i], buffer->embed_pot, particleSpecies[i], ntypes, rho_val[i]);

      if (comp_energy)
	*energy += phi;

      if (comp_particleEnergy)
	particleEnergy[i] += phi;
    }
  }				/* loop over atoms */

  /* second loop over atoms, calculates eam forces */
  for (i = 0; i < *nAtoms; i++) {
    request = i - model_index_shift;
    ier = (*get_neigh) (km, &one, &request, &currentAtom,
      &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
    if (KIM_STATUS_OK != ier) {
      KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
      ier = KIM_STATUS_FAIL;
      return ier;
    }
    it = particleSpecies[i];

    /* loop over all atoms in the neighbor list */
    for (jj = 0; jj < numOfAtomNeigh; jj++) {
      j = neighListOfCurrentAtom[jj] + model_index_shift;

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = jt * ntypes + it;
      col2 = it * ntypes + jt;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = Rij_list[jj * DIM + l];
	r2 += (Rij[l] * Rij[l]);
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if ((r2 < buffer->transfer_pot.end[col1]) || (r2 < buffer->transfer_pot.end[col2])) {

	/* caculate derivative of the density function for both atoms */
	DERIV_FUNC(rho_i_prime, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(2, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	if (col1 == col2) {
	  rho_j_prime = rho_i_prime;
	} else {
	  DERIV_FUNC(rho_j_prime, buffer->transfer_pot, col2, inc, r2, is_short);
	  if (is_short) {
	    short_dist_warning(2, i, j, particleSpecies, Rij, R);
	    is_short = 0;
	  }
	}

	/* combine all contributions, dF_val[] is too big by a factor of 2 */
	dphi = 0.5 * dF_val[i] * rho_j_prime;

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &j);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }
    }				/* loop over neighbors */
  }				/* loop over atoms */

  return KIM_STATUS_OK;
}

/****************************************************************
 *
 * compute function for neighbor list = NEIGH_RVEC_F
 * and iterator mode
 *
 ****************************************************************/

int compute_rvec_full_iter(void *km)
{
  intptr_t *pkim = *((intptr_t **) km);

  /* general variables */
  int ier;
  int   i, j, jj, l;

  /* pointer to model buffer */
  model_buffer *buffer;

  /* flags for available computes */
  int   comp_energy;
  int   comp_force;
  int   comp_particleEnergy;
  int   comp_process_dEdr;

  /* pointers to objects in the KIM API */
  typedef int   (*get_neigh_ptr) (void *, int *, int *, int *, int *, int **, double **);
  get_neigh_ptr get_neigh;
  int  *nAtoms;
  int  *particleSpecies;
  double *energy;
  double *force;
  double *particleEnergy;

  /* variables / pointers for objects in the buffer */
  int   model_index_shift;
  int   ntypes;
  double *dF_val;
  double *rho_val;

  /* variables for dealing with the KIM API */
  int   currentAtom;
  int  *neighListOfCurrentAtom;
  int   numOfAtomNeigh;
  int   one = 1;
  int   zero = 0;
  double Rij[DIM];
  double *pRij = &(Rij[0]);
  double *Rij_list;

  /* variables used for calculations */
  int   col1, col2;
  int   inc;
  int   is_short = 0;
  int   it, jt;
  double r2;
  double rho = 0.0;
  double rho_i_prime, rho_j_prime;
  double R;
  double phi;
  double dphi;

  /* get buffer from KIM object */
  buffer = (model_buffer *) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* unpack info from the buffer */
  model_index_shift = buffer->model_index_shift;

  /* check to see if we have been asked to compute the energy, forces, and particleEnergies */
  /* *INDENT-OFF* */
  KIM_API_getm_compute_by_index(pkim, &ier, 4 * 3,
		  	buffer->energy_ind, 		&comp_energy, 		1,
		  	buffer->forces_ind, 		&comp_force, 		1,
			buffer->particleEnergy_ind, 	&comp_particleEnergy, 	1,
			buffer->process_dEdr_ind, 	&comp_process_dEdr, 	1);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute_by_index", ier);
    return ier;
  }

  /* pull out information by index */
  /* *INDENT-OFF* */
  KIM_API_getm_data_by_index(pkim, &ier, 5 * 3,
                              buffer->numberOfParticles_ind, 	&nAtoms,         	1,
                              buffer->particleSpecies_ind,        &particleSpecies,  	1,
                              buffer->energy_ind,               &energy,         	comp_energy,
                              buffer->forces_ind,               &force,          	comp_force,
                              buffer->particleEnergy_ind,       &particleEnergy, 	comp_particleEnergy);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data_by_index", ier);
    return ier;
  }
  get_neigh = (get_neigh_ptr) KIM_API_get_method_by_index(pkim, buffer->get_neigh_ind, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_method_by_index", ier);
    return ier;
  }

  ntypes = buffer->ntypes;
  inc = ntypes * ntypes;

  /* Check to be sure that the atom types are correct */
  for (i = 0; i < *nAtoms; ++i)
    if (particleSpecies[i] > ntypes) {
      KIM_API_report_error(__LINE__, __FILE__, "Unexpected species type detected", KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }

  /* check if the arrays for rho and dF are big enough */
  if (buffer->table_len < *nAtoms) {
    buffer->rho_val = (double *)realloc(buffer->rho_val, (*nAtoms) * sizeof(double));
    buffer->dF_val = (double *)realloc(buffer->dF_val, (*nAtoms) * sizeof(double));
    if (NULL == buffer->rho_val || NULL == buffer->dF_val) {
      KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }
    buffer->table_len = *nAtoms;
  }
  rho_val = buffer->rho_val;
  dF_val = buffer->dF_val;

  /* reset arrays */
  for (i = 0; i < *nAtoms; i++) {
    rho_val[i] = 0.0;
    dF_val[i] = 0.0;
  }
  if (comp_energy)
    *energy = 0.0;
  if (comp_force)
    for (i = 0; i < DIM * (*nAtoms); i++)
      force[i] = 0.0;
  if (comp_particleEnergy)
    for (i = 0; i < *nAtoms; ++i)
      particleEnergy[i] = 0.0;

  /* initialize the iterator */
  ier = (*get_neigh) (km, &zero, &zero, &currentAtom, &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
  if (KIM_STATUS_NEIGH_ITER_INIT_OK != ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", ier);
    return KIM_STATUS_FAIL;
  }

  i = -1;
  /* first infinite iterator loop, calculates pair energies and forces as well as eam energies */
  while (1) {

    /* get the next neighbor */
    ier = (*get_neigh) (km, &zero, &one, &currentAtom, &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
    /* abort the infinite loop if we are past the end of the list */
    if (KIM_STATUS_NEIGH_ITER_PAST_END == ier)
      break;
    if (KIM_STATUS_OK != ier) {
      KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
      return KIM_STATUS_FAIL;
    }
    i = currentAtom + model_index_shift;
    it = particleSpecies[i];

    /* loop over all atoms in the neighbor list */
    for (jj = 0; jj < numOfAtomNeigh; jj++) {
      j = neighListOfCurrentAtom[jj] + model_index_shift;

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = it * ntypes + jt;
      col2 = jt * ntypes + it;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = Rij_list[jj * DIM + l];
	r2 += Rij[l] * Rij[l];
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if (r2 <= buffer->pair_pot.end[col1]) {

	/* calculate pair potential and the derivative at r2 */
	PAIR_INT(phi, dphi, buffer->pair_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(0, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	/* we have a full neighbor list, so we only contribute 0.5 */
	phi *= 0.5;
	dphi *= 0.5;

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_energy)
	  *energy += phi;

	if (comp_particleEnergy)
	  particleEnergy[i] += phi;

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &j);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }

      /* calculate contribution to density */
      if (r2 < buffer->transfer_pot.end[col1]) {
	VAL_FUNC(rho, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(1, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	rho_val[i] += rho;
      }
      /* do NOT assign to neighbor (full neighbor list) */
    }				/* loop over neighbors */

    /* calculate the embedding energies for all contributing atoms */
    if (numOfAtomNeigh > 0) {
      PAIR_INT2(phi, dF_val[i], buffer->embed_pot, particleSpecies[i], ntypes, rho_val[i]);

      if (comp_energy)
	*energy += phi;

      if (comp_particleEnergy)
	particleEnergy[i] += phi;
    }
  }				/* infinite loop */

  /* initialize the iterator */
  ier = (*get_neigh) (km, &zero, &zero, &currentAtom, &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
  if (KIM_STATUS_NEIGH_ITER_INIT_OK != ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", ier);
    return KIM_STATUS_FAIL;
  }

  i = -1;
  /* second infinite iterator loop, calculates eam forces */
  while (1) {

    /* get the next neighbor */
    ier = (*get_neigh) (km, &zero, &one, &currentAtom, &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
    /* abort the infinite loop if we are past the end of the list */
    if (KIM_STATUS_NEIGH_ITER_PAST_END == ier)
      break;
    if (KIM_STATUS_OK != ier) {
      KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
      return KIM_STATUS_FAIL;
    }
    i = currentAtom + model_index_shift;
    it = particleSpecies[i];

    /* loop over all atoms in the neighbor list */
    for (jj = 0; jj < numOfAtomNeigh; jj++) {
      j = neighListOfCurrentAtom[jj] + model_index_shift;

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = jt * ntypes + it;
      col2 = it * ntypes + jt;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = Rij_list[jj * DIM + l];
	r2 += Rij[l] * Rij[l];
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if ((r2 < buffer->transfer_pot.end[col1]) || (r2 < buffer->transfer_pot.end[col2])) {

	/* caculate derivative of the density function for both atoms */
	DERIV_FUNC(rho_i_prime, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(2, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	if (col1 == col2) {
	  rho_j_prime = rho_i_prime;
	} else {
	  DERIV_FUNC(rho_j_prime, buffer->transfer_pot, col2, inc, r2, is_short);
	  if (is_short) {
	    short_dist_warning(2, i, j, particleSpecies, Rij, R);
	    is_short = 0;
	  }
	}

	/* combine all contributions, dF_val[] is too big by a factor of 2 */
	dphi = 0.5 * dF_val[i] * rho_j_prime;

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &j);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }
    }				/* loop over neighbors */
  }				/* infinite loop */

  return KIM_STATUS_OK;
}

/****************************************************************
 *
 * compute function for neighbor list = NEIGH_RVEC_H
 * and locator mode
 *
 ****************************************************************/

int compute_rvec_half_loca(void *km)
{
  intptr_t *pkim = *((intptr_t **) km);

  /* general variables */
  int ier;
  int   i, j, jj, l;

  /* pointer to model buffer */
  model_buffer *buffer;

  /* flags for available computes */
  int   comp_energy;
  int   comp_force;
  int   comp_particleEnergy;
  int   comp_process_dEdr;

  /* pointers to objects in the KIM API */
  typedef int   (*get_neigh_ptr) (void *, int *, int *, int *, int *, int **, double **);
  get_neigh_ptr get_neigh;
  int  *nAtoms;
  int  *numContrib;
  int  *particleSpecies;
  double *energy;
  double *force;
  double *particleEnergy;

  /* variables / pointers for objects in the buffer */
  int   model_index_shift;
  int   ntypes;
  double *dF_val;
  double *rho_val;

  /* variables for dealing with the KIM API */
  int   currentAtom;
  int  *neighListOfCurrentAtom;
  int   numOfAtomNeigh;
  int   one = 1;
  int   request;
  double Rij[DIM];
  double *pRij = &(Rij[0]);
  double *Rij_list;

  /* variables used for calculations */
  int   col1, col2;
  int   inc;
  int   is_short = 0;
  int   it, jt;
  int 	numberContrib;
  double r2;
  double rho = 0.0;
  double rho_i_prime, rho_j_prime;
  double R;
  double phi;
  double dphi;

  /* get buffer from KIM object */
  buffer = (model_buffer *) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* unpack info from the buffer */
  model_index_shift = buffer->model_index_shift;

  /* check to see if we have been asked to compute the energy, forces, and particleEnergies */
  /* *INDENT-OFF* */
  KIM_API_getm_compute_by_index(pkim, &ier, 4 * 3,
		  buffer->energy_ind,    	&comp_energy, 		1,
		  buffer->forces_ind, 		&comp_force, 		1,
		  buffer->particleEnergy_ind, 	&comp_particleEnergy, 	1,
		  buffer->process_dEdr_ind, 	&comp_process_dEdr, 	1);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute_by_index", ier);
    return ier;
  }

  /* pull out information by index */
  /* *INDENT-OFF* */
  KIM_API_getm_data_by_index(pkim, &ier, 6 * 3,
		  buffer->numberOfParticles_ind,    		&nAtoms, 		1,
		  buffer->numberContributingParticles_ind, 	&numContrib, 		1,
		  buffer->particleSpecies_ind, 			&particleSpecies, 	1,
		  buffer->energy_ind, 				&energy, 		comp_energy,
		  buffer->forces_ind, 				&force, 		comp_force,
		  buffer->particleEnergy_ind, 			&particleEnergy, 	comp_particleEnergy);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data_by_index", ier);
    return ier;
  }
  get_neigh = (get_neigh_ptr) KIM_API_get_method_by_index(pkim, buffer->get_neigh_ind, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_method_by_index", ier);
    return ier;
  }

  numberContrib = *numContrib;
  ntypes = buffer->ntypes;
  inc = ntypes * ntypes;

  /* Check to be sure that the atom types are correct */
  for (i = 0; i < *nAtoms; ++i)
    if (particleSpecies[i] > ntypes) {
      KIM_API_report_error(__LINE__, __FILE__, "Unexpected species type detected", KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }

  /* check if the arrays for rho and dF are big enough */
  if (buffer->table_len < *nAtoms) {
    buffer->rho_val = (double *)realloc(buffer->rho_val, (*nAtoms) * sizeof(double));
    buffer->dF_val = (double *)realloc(buffer->dF_val, (*nAtoms) * sizeof(double));
    if (NULL == buffer->rho_val || NULL == buffer->dF_val) {
      KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }
    buffer->table_len = *nAtoms;
  }
  rho_val = buffer->rho_val;
  dF_val = buffer->dF_val;

  /* reset arrays */
  for (i = 0; i < *nAtoms; i++) {
    rho_val[i] = 0.0;
    dF_val[i] = 0.0;
  }
  if (comp_energy)
    *energy = 0.0;
  if (comp_force)
    for (i = 0; i < DIM * (*nAtoms); i++)
      force[i] = 0.0;
  if (comp_particleEnergy)
    for (i = 0; i < *nAtoms; ++i)
      particleEnergy[i] = 0.0;

  /* first loop over atoms, calculates pair energies and forces as well as eam energies */
  for (i = 0; i < numberContrib; i++) {
    request = i - model_index_shift;
    ier = (*get_neigh) (km, &one, &request, &currentAtom,
      &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
    if (KIM_STATUS_OK != ier) {
      KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
      return KIM_STATUS_FAIL;
    }
    it = particleSpecies[i];

    /* loop over all atoms in the neighbor list */
    for (jj = 0; jj < numOfAtomNeigh; jj++) {
      j = neighListOfCurrentAtom[jj] + model_index_shift;

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = it * ntypes + jt;
      col2 = jt * ntypes + it;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = Rij_list[jj * DIM + l];
	r2 += Rij[l] * Rij[l];
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if (r2 <= buffer->pair_pot.end[col1]) {

	/* calculate pair potential and the derivative at r2 */
	PAIR_INT(phi, dphi, buffer->pair_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(0, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	phi *= 0.5;
	if (j >= numberContrib)
	  dphi *= 0.5;

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_energy) {
	  *energy += phi;
	  if (j < numberContrib)
	    *energy += phi;
	}

	if (comp_particleEnergy) {
	  particleEnergy[i] += phi;
	  particleEnergy[j] += phi;
	}

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &j);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }

      /* calculate contribution to density */
      if (r2 < buffer->transfer_pot.end[col1]) {
	VAL_FUNC(rho, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(1, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	rho_val[i] += rho;
      }
      /* assign also to neighbor (half neighbor list) */
      if (it == jt) {
	if (r2 < buffer->transfer_pot.end[col2])
	  rho_val[j] += rho;
      } else {
	if (r2 < buffer->transfer_pot.end[col2]) {
	  VAL_FUNC(rho, buffer->transfer_pot, col2, inc, r2, is_short);
	  if (is_short) {
	    short_dist_warning(1, i, j, particleSpecies, Rij, R);
	    is_short = 0;
	  }
	  rho_val[j] += rho;
	}
      }
    }				/* loop over neighbors */

    /* calculate the embedding energies for all contributing atoms */
    if (i < numberContrib) {
      PAIR_INT2(phi, dF_val[i], buffer->embed_pot, particleSpecies[i], ntypes, rho_val[i]);

      if (comp_energy)
	*energy += phi;

      if (comp_particleEnergy)
	particleEnergy[i] += phi;
    }
  }				/* loop over atoms */

  /* second loop over atoms, calculates eam forces */
  for (i = 0; i < numberContrib; i++) {
    request = i - model_index_shift;
    ier = (*get_neigh) (km, &one, &request, &currentAtom,
      &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
    if (KIM_STATUS_OK != ier) {
      KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
      return KIM_STATUS_FAIL;
    }
    it = particleSpecies[i];

    /* loop over all atoms in the neighbor list */
    for (jj = 0; jj < numOfAtomNeigh; jj++) {
      j = neighListOfCurrentAtom[jj] + model_index_shift;

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = jt * ntypes + it;
      col2 = it * ntypes + jt;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = Rij_list[jj * DIM + l];
	r2 += (Rij[l] * Rij[l]);
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if ((r2 < buffer->transfer_pot.end[col1]) || (r2 < buffer->transfer_pot.end[col2])) {

	/* caculate derivative of the density function for both atoms */
	DERIV_FUNC(rho_i_prime, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(2, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	if (col1 == col2) {
	  rho_j_prime = rho_i_prime;
	} else {
	  DERIV_FUNC(rho_j_prime, buffer->transfer_pot, col2, inc, r2, is_short);
	  if (is_short) {
	    short_dist_warning(2, i, j, particleSpecies, Rij, R);
	    is_short = 0;
	  }
	}

	/* combine all contributions, dF_val[] is too big by a factor of 2 */
	dphi = 0.5 * (dF_val[i] * rho_j_prime + dF_val[j] * rho_i_prime);

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &j);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }
    }				/* loop over neighbors */
  }				/* loop over atoms */

  return KIM_STATUS_OK;
}

/****************************************************************
 *
 * compute function for neighbor list = NEIGH_RVEC_H
 * and iterator mode
 *
 ****************************************************************/

int compute_rvec_half_iter(void *km)
{
  intptr_t *pkim = *((intptr_t **) km);

  /* general variables */
  int ier;
  int   i, j, jj, l;

  /* pointer to model buffer */
  model_buffer *buffer;

  /* flags for available computes */
  int   comp_energy;
  int   comp_force;
  int   comp_particleEnergy;
  int   comp_process_dEdr;

  /* pointers to objects in the KIM API */
  typedef int   (*get_neigh_ptr) (void *, int *, int *, int *, int *, int **, double **);
  get_neigh_ptr get_neigh;
  int  *nAtoms;
  int  *numContrib;
  int  *particleSpecies;
  double *energy;
  double *force;
  double *particleEnergy;

  /* variables / pointers for objects in the buffer */
  int   model_index_shift;
  int   ntypes;
  double *dF_val;
  double *rho_val;

  /* variables for dealing with the KIM API */
  int   currentAtom;
  int  *neighListOfCurrentAtom;
  int   numOfAtomNeigh;
  int   one = 1;
  int   zero = 0;
  double Rij[DIM];
  double *pRij = &(Rij[0]);
  double *Rij_list;

  /* variables used for calculations */
  int   col1, col2;
  int   inc;
  int   is_short = 0;
  int   it, jt;
  int 	numberContrib;
  double r2;
  double rho = 0.0;
  double rho_i_prime, rho_j_prime;
  double R;
  double phi;
  double dphi;

  /* get buffer from KIM object */
  buffer = (model_buffer *) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* unpack info from the buffer */
  model_index_shift = buffer->model_index_shift;

  /* check to see if we have been asked to compute the energy, forces, and particleEnergies */
  /* *INDENT-OFF* */
  KIM_API_getm_compute_by_index(pkim, &ier, 4 * 3,
		  	buffer->energy_ind, 		&comp_energy, 		1,
		  	buffer->forces_ind, 		&comp_force, 		1,
			buffer->particleEnergy_ind, 	&comp_particleEnergy, 	1,
			buffer->process_dEdr_ind, 	&comp_process_dEdr, 	1);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute_by_index", ier);
    return ier;
  }

  /* pull out information by index */
  /* *INDENT-OFF* */
  KIM_API_getm_data_by_index(pkim, &ier, 6 * 3,
                              buffer->numberOfParticles_ind,           &nAtoms,         	1,
			      buffer->numberContributingParticles_ind, &numContrib, 		1,
                              buffer->particleSpecies_ind,               &particleSpecies,  	1,
                              buffer->energy_ind,                      &energy,         	comp_energy,
                              buffer->forces_ind,                      &force,          	comp_force,
                              buffer->particleEnergy_ind,              &particleEnergy, 	comp_particleEnergy);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data_by_index", ier);
    return ier;
  }
  get_neigh = (get_neigh_ptr) KIM_API_get_method_by_index(pkim, buffer->get_neigh_ind, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_method_by_index", ier);
    return ier;
  }

  numberContrib = *numContrib;
  ntypes = buffer->ntypes;
  inc = ntypes * ntypes;

  /* Check to be sure that the atom types are correct */
  for (i = 0; i < *nAtoms; ++i)
    if (particleSpecies[i] > ntypes) {
      KIM_API_report_error(__LINE__, __FILE__, "Unexpected species type detected", KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }

  /* check if the arrays for rho and dF are big enough */
  if (buffer->table_len < *nAtoms) {
    buffer->rho_val = (double *)realloc(buffer->rho_val, (*nAtoms) * sizeof(double));
    buffer->dF_val = (double *)realloc(buffer->dF_val, (*nAtoms) * sizeof(double));
    if (NULL == buffer->rho_val || NULL == buffer->dF_val) {
      KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }
    buffer->table_len = *nAtoms;
  }
  rho_val = buffer->rho_val;
  dF_val = buffer->dF_val;

  /* reset arrays */
  for (i = 0; i < *nAtoms; i++) {
    rho_val[i] = 0.0;
    dF_val[i] = 0.0;
  }
  if (comp_energy)
    *energy = 0.0;
  if (comp_force)
    for (i = 0; i < DIM * (*nAtoms); i++)
      force[i] = 0.0;
  if (comp_particleEnergy)
    for (i = 0; i < *nAtoms; ++i)
      particleEnergy[i] = 0.0;

  /* initialize the iterator */
  ier = (*get_neigh) (km, &zero, &zero, &currentAtom, &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
  if (KIM_STATUS_NEIGH_ITER_INIT_OK != ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", ier);
    return KIM_STATUS_FAIL;
  }

  i = -1;
  /* first infinite iterator loop, calculates pair energies and forces as well as eam energies */
  while (1) {

    /* get the next neighbor */
    ier = (*get_neigh) (km, &zero, &one, &currentAtom, &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
    /* abort the infinite loop if we are past the end of the list */
    if (KIM_STATUS_NEIGH_ITER_PAST_END == ier)
      break;
    if (KIM_STATUS_OK != ier) {
      KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
      return KIM_STATUS_FAIL;
    }
    i = currentAtom + model_index_shift;
    it = particleSpecies[i];

    /* loop over all atoms in the neighbor list */
    for (jj = 0; jj < numOfAtomNeigh; jj++) {
      j = neighListOfCurrentAtom[jj] + model_index_shift;

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = it * ntypes + jt;
      col2 = jt * ntypes + it;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = Rij_list[jj * DIM + l];
	r2 += Rij[l] * Rij[l];
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if (r2 <= buffer->pair_pot.end[col1]) {

	/* calculate pair potential and the derivative at r2 */
	PAIR_INT(phi, dphi, buffer->pair_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(0, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	phi *= 0.5;
	if (j >= numberContrib)
	  dphi *= 0.5;

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_energy) {
	  *energy += phi;
	  if (j < numberContrib)
	    *energy += phi;
	}

	if (comp_particleEnergy) {
	  particleEnergy[i] += phi;
	  particleEnergy[j] += phi;
	}

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &j);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }

      /* calculate contribution to density */
      if (r2 < buffer->transfer_pot.end[col1]) {
	VAL_FUNC(rho, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(1, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	rho_val[i] += rho;
      }
      /* assign also to neighbor (half neighbor list) */
      if (it == jt) {
	if (r2 < buffer->transfer_pot.end[col2])
	  rho_val[j] += rho;
      } else {
	if (r2 < buffer->transfer_pot.end[col2]) {
	  VAL_FUNC(rho, buffer->transfer_pot, col2, inc, r2, is_short);
	  if (is_short) {
	    short_dist_warning(1, i, j, particleSpecies, Rij, R);
	    is_short = 0;
	  }
	  rho_val[j] += rho;
	}
      }
    }				/* loop over neighbors */

    /* calculate the embedding energies for all contributing atoms */
    if (i < numberContrib) {
      PAIR_INT2(phi, dF_val[i], buffer->embed_pot, particleSpecies[i], ntypes, rho_val[i]);

      if (comp_energy)
	*energy += phi;

      if (comp_particleEnergy)
	particleEnergy[i] += phi;
    }
  }				/* infinite loop */

  /* initialize the iterator */
  ier = (*get_neigh) (km, &zero, &zero, &currentAtom, &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
  if (KIM_STATUS_NEIGH_ITER_INIT_OK != ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", ier);
    return KIM_STATUS_FAIL;
  }

  i = -1;
  /* second infinite iterator loop, calculates eam forces */
  while (1) {

    /* get the next neighbor */
    ier = (*get_neigh) (km, &zero, &one, &currentAtom, &numOfAtomNeigh, &neighListOfCurrentAtom, &Rij_list);
    /* abort the infinite loop if we are past the end of the list */
    if (KIM_STATUS_NEIGH_ITER_PAST_END == ier)
      break;
    if (KIM_STATUS_OK != ier) {
      KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
      return KIM_STATUS_FAIL;
    }
    i = currentAtom + model_index_shift;
    it = particleSpecies[i];

    /* loop over all atoms in the neighbor list */
    for (jj = 0; jj < numOfAtomNeigh; jj++) {
      j = neighListOfCurrentAtom[jj] + model_index_shift;

      /* set up particle types and cols */
      jt = particleSpecies[j];
      col1 = jt * ntypes + it;
      col2 = it * ntypes + jt;

      /* calculate distance */
      r2 = 0.0;
      for (l = 0; l < DIM; l++) {
	Rij[l] = Rij_list[jj * DIM + l];
	r2 += Rij[l] * Rij[l];
      }
      R = sqrt(r2);

      /* check if we are within the cutoff radius */
      if ((r2 < buffer->transfer_pot.end[col1]) || (r2 < buffer->transfer_pot.end[col2])) {

	/* caculate derivative of the density function for both atoms */
	DERIV_FUNC(rho_i_prime, buffer->transfer_pot, col1, inc, r2, is_short);
	if (is_short) {
	  short_dist_warning(2, i, j, particleSpecies, Rij, R);
	  is_short = 0;
	}
	if (col1 == col2) {
	  rho_j_prime = rho_i_prime;
	} else {
	  DERIV_FUNC(rho_j_prime, buffer->transfer_pot, col2, inc, r2, is_short);
	  if (is_short) {
	    short_dist_warning(2, i, j, particleSpecies, Rij, R);
	    is_short = 0;
	  }
	}

	/* combine all contributions, dF_val[] is too big by a factor of 2 */
	dphi = 0.5 * (dF_val[i] * rho_j_prime + dF_val[j] * rho_i_prime);

	if (comp_force)
	  for (l = 0; l < DIM; l++) {
	    force[i * DIM + l] += Rij[l] * dphi;
	    force[j * DIM + l] -= Rij[l] * dphi;
	  }

	if (comp_process_dEdr) {
	  dphi *= R;
	  ier = KIM_API_process_dEdr(km, &dphi, &R, &pRij, &i, &j);
	  if (KIM_STATUS_OK > ier) {
	    KIM_API_report_error(__LINE__, __FILE__, "error in process_dEdr routine", ier);
	    return ier;
	  }
	}
      }
    }				/* loop over neighbors */
  }				/* infinite loop */

  return KIM_STATUS_OK;
}

/****************************************************************
 *
 * model driver initialization function
 * this is called by the KIM API and reads the parameter file
 * of the model implementing this model driver
 *
 ****************************************************************/

int model_driver_init(void *km, char* paramfile_names, int* nmstrlen, int* numparamfiles)
{
  intptr_t *pkim = *((intptr_t **) km);

  /* local variables */
  char  msg[255];
  const char *NBCstr;
  int   i, ier = KIM_STATUS_OK;
  int   ntypes;
  int   dummy;
  double *model_cutoff = NULL;
  model_buffer *buffer;

  /* neighbor lists */
  int   NBC;
  int   IterOrLoca;
  int   HalfOrFull;

  /* strings to store the filenames from the parameter file */
  char *pairpot = &(paramfile_names[0]);
  char *transfer = &(paramfile_names[(*nmstrlen)]);
  char *embed = &(paramfile_names[2*(*nmstrlen)]);

  if (*numparamfiles != 3) {
     KIM_API_report_error(__LINE__, __FILE__, "Wrong number of input files", KIM_STATUS_FAIL);
     return KIM_STATUS_FAIL;
  }

  /* get number of types */
  ier = KIM_API_get_num_model_species(pkim, &ntypes, &dummy);
  if (KIM_STATUS_OK > ier) {
     KIM_API_report_error(__LINE__, __FILE__, "Unsupported number of particle types or error getting number of types.", KIM_STATUS_FAIL);
     return KIM_STATUS_FAIL;
  }

  /* store pointer to function in KIM object */
  /* *INDENT-OFF* */
  KIM_API_setm_method(pkim, &ier, 3 * 4,
		  "compute", 	1, 	&compute, 	1,
		  "destroy", 	1,  	&destroy, 	1,
		  "reinit", 	1, 	&reinit, 	1);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_setm_method", ier);
    return KIM_STATUS_FAIL;
  }

  /* allocate buffer */
  buffer = (model_buffer *) malloc(sizeof(model_buffer));
  if (NULL == buffer) {
    KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
    return KIM_STATUS_FAIL;
  }

  /* setup buffer */
  /* Determine neighbor list boundary condition (NBC) */
  /**************************************************************
   * NBCstr =  	0 -> cluster mode
   * 		1 -> pure mode (half or full)
   * 		2 -> rvec (only full)
   * 		3 -> mi_opbc (CURRENTLY NOT SUPPORTED)
   **************************************************************/
  ier = KIM_API_get_NBC_method(pkim, &NBCstr);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_NBC_method", ier);
    return ier;
  }
  if (!strcmp("CLUSTER", NBCstr)) {
    NBC = 0;
  } else if ((!strcmp("NEIGH_PURE_H", NBCstr)) || (!strcmp("NEIGH_PURE_F", NBCstr))) {
    NBC = 1;
  } else if ((!strcmp("NEIGH_RVEC_H", NBCstr)) || (!strcmp("NEIGH_RVEC_F", NBCstr))) {
    NBC = 2;
  } else {
    KIM_API_report_error(__LINE__, __FILE__, "Unknown NBC method", KIM_STATUS_FAIL);
    return KIM_STATUS_FAIL;
  }

  /* Determine if Half or Full neighbor lists are being used */
  /**************************************************************
   * HalfOrFull = 1 -- Half
   *            = 2 -- Full
   **************************************************************/
  if (KIM_API_is_half_neighbors(pkim, &ier)) {
    HalfOrFull = 1;
  } else {
    HalfOrFull = 2;
  }

  /* determine neighbor list handling mode */
  if (0 != NBC) {
  /**************************************************************
   * IterOrLoca = 1 -- Iterator
   *            = 2 -- Locator
   **************************************************************/
    IterOrLoca = KIM_API_get_neigh_mode(pkim, &ier);
    if (KIM_STATUS_OK > ier) {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh_mode", ier);
      return ier;
    }
    if ((IterOrLoca != 1) && (IterOrLoca != 2)) {
      sprintf(msg, "Unsupported IterOrLoca mode = %i\n", IterOrLoca);
      KIM_API_report_error(__LINE__, __FILE__, msg, KIM_STATUS_FAIL);
      return KIM_STATUS_FAIL;
    }
  } else {
    IterOrLoca = 2;		/* for CLUSTER NBC */
  }

  /* assign correct compute routine */
  /* cluster mode */
  if (0 == NBC) {
    buffer->compute = compute_cluster;
    /* pure neighbor lists */
  } else if (1 == NBC) {
    if (1 == HalfOrFull && 2 == IterOrLoca)
      buffer->compute = compute_pure_half_loca;
    else if (1 == HalfOrFull && 1 == IterOrLoca)
      buffer->compute = compute_pure_half_iter;
    else if (2 == HalfOrFull && 2 == IterOrLoca)
      buffer->compute = compute_pure_full_loca;
    else if (2 == HalfOrFull && 1 == IterOrLoca)
      buffer->compute = compute_pure_full_iter;
    /* rvec neighbor lists */
  } else if (2 == NBC) {
    if (1 == HalfOrFull && 2 == IterOrLoca)
      buffer->compute = compute_rvec_half_loca;
    else if (1 == HalfOrFull && 1 == IterOrLoca)
      buffer->compute = compute_rvec_half_iter;
    else if (2 == HalfOrFull && 2 == IterOrLoca)
      buffer->compute = compute_rvec_full_loca;
    else if (2 == HalfOrFull && 1 == IterOrLoca)
      buffer->compute = compute_rvec_full_iter;
  }

  /* store the model_index_shift in the buffer */
  buffer->model_index_shift = KIM_API_get_model_index_shift(pkim);

  /* read the tabulated pair potential file */
  read_pot_table(&(buffer->pair_pot), pairpot, ntypes * ntypes, ntypes, 1);
  /* read the tabulated electron density function */
  read_pot_table(&(buffer->transfer_pot), transfer, ntypes * ntypes, ntypes, 1);
  /* read the tabulated embedding energy function */
  read_pot_table(&(buffer->embed_pot), embed, ntypes, ntypes, 0);

  /* get model_cutoff pointer from KIM object */
  model_cutoff = (double *)KIM_API_get_data(pkim, "cutoff", &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data(\"cutoff\")", ier);
    return ier;
  }

  /* calculate the cutoff */
  /* the cutoff is the maximum cutoff of all potentials
   * this is needed so the model can give us all needed neighbors
   * the compute routines check for the cutoff of the different potentials
   */
  for (i = 0; i < ntypes * ntypes; i++)
    *model_cutoff = MAX(*model_cutoff, buffer->pair_pot.end[i]);
  for (i = 0; i < ntypes; i++)
    *model_cutoff = MAX(*model_cutoff, buffer->embed_pot.end[i]);
  for (i = 0; i < ntypes * ntypes; i++)
    *model_cutoff = MAX(*model_cutoff, buffer->transfer_pot.end[i]);
  /* the distance in IMD is r^2, we need the cutoff in angstrom */
  *model_cutoff = sqrt(*model_cutoff);

  /* set the pointer in the buffer and assign the correct value */
  buffer->ntypes = ntypes;

  /* allocate the arrays for the density and embedding values */
  buffer->rho_val = (double *)malloc(1 * sizeof(double));
  buffer->dF_val = (double *)malloc(1 * sizeof(double));
  if (NULL == buffer->rho_val || NULL == buffer->dF_val) {
    KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
    return KIM_STATUS_FAIL;
  }
  buffer->table_len = 1;

  /* get the indices of several data structures within the KIM API for easy access */
  /* *INDENT-OFF* */
  KIM_API_getm_index(pkim, &ier, 10 * 3,
		"coordinates", 	   		&(buffer->coordinates_ind), 			1,
		"energy", 			&(buffer->energy_ind), 				1,
		"forces",    			&(buffer->forces_ind), 				1,
		"get_neigh",                    &(buffer->get_neigh_ind),                   	1,
		"numberContributingParticles", 	&(buffer->numberContributingParticles_ind), 	1,
		"numberOfParticles", 		&(buffer->numberOfParticles_ind), 		1,
		"numberOfSpecies", 		&(buffer->numberOfSpecies_ind), 		1,
		"particleEnergy", 		&(buffer->particleEnergy_ind),    		1,
		"process_dEdr",                &(buffer->process_dEdr_ind),      1,
		"particleSpecies", 		&(buffer->particleSpecies_ind), 			1);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_index", ier);
    return ier;
  }

  /* store in model buffer */
  KIM_API_set_model_buffer(pkim, (void *)buffer, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_set_model_buffer", ier);
    return ier;
  }

  return KIM_STATUS_OK;
}

/****************************************************************
 *
 * read potential table; choose format according to header
 *
 ****************************************************************/

void read_pot_table(pot_table_t *pt, char *filename, int ncols, int ntypes, int radial)
{
  FILE *infile = NULL;
  char  buffer[1024], msg[255];
  char *res;
  int   have_header = 0, have_format = 0, end_header = 0;
  int   size = ncols;
  int   format = 2;		/* 2 for EAM2, 1 otherwise */
  int   i;

  /* read header */
  /* open file */
  infile = fopen(filename, "r");
  if (NULL == infile) {
    sprintf(msg, "Could not open potential file:\n\t\t %s", filename);
    KIM_API_report_error(__LINE__, __FILE__, msg, KIM_STATUS_FAIL);
    exit(EXIT_FAILURE);
  }

  /* read the header */
  do {
    /* read one line */
    res = fgets(buffer, 1024, infile);
    if (NULL == res) {
      sprintf(msg, "Unexpected end of file in %s", filename);
      KIM_API_report_error(__LINE__, __FILE__, msg, KIM_STATUS_FAIL);
      exit(EXIT_FAILURE);
    }
    /* see if it is a header line */
    if (buffer[0] == '#') {
      have_header = 1;
      /* stop after last header line */
      end_header = (buffer[1] == 'E');
      /* see if it is the format line */
      if (buffer[1] == 'F') {
	/* format complete? */
	if (2 != sscanf((const char *)(buffer + 2), "%d%d", &format, &size)) {
	  sprintf(msg, "Corrupted format header line in file %s", filename);
	  KIM_API_report_error(__LINE__, __FILE__, msg, KIM_STATUS_FAIL);
	  exit(EXIT_FAILURE);
	}
	/* right number of columns? */
	if (size != ncols) {
	  sprintf(msg, "Wrong number of data columns in file %%s\nShould be %d, is %d", ncols, size);
	  KIM_API_report_error(__LINE__, __FILE__, msg, KIM_STATUS_FAIL);
	  exit(EXIT_FAILURE);
	}
	/* recognized format? */
	if ((format != 1) && (format != 2)) {
	  sprintf(msg, "Unrecognized format specified for file %s", filename);
	  KIM_API_report_error(__LINE__, __FILE__, msg, KIM_STATUS_FAIL);
	  exit(EXIT_FAILURE);
	}
	have_format = 1;
      }
    } else if (have_header) {
      /* header does not end properly */
      sprintf(msg, "Corrupted header in file %s", filename);
      KIM_API_report_error(__LINE__, __FILE__, msg, KIM_STATUS_FAIL);
      exit(EXIT_FAILURE);
    } else {
      /* we have no header, stop reading further */
      end_header = 1;
    }
  } while (!end_header);

  /* did we have a format in the header */
  if ((have_header) && (!have_format)) {
    sprintf(msg, "Format not specified in header of file %s", filename);
    KIM_API_report_error(__LINE__, __FILE__, msg, KIM_STATUS_FAIL);
    exit(EXIT_FAILURE);
  }

  /* rewind if there was no header */
  if (!have_header)
    rewind(infile);

  /* warn if we have no header */
  if (!have_header) {
    fprintf(stderr, "WARNING: File %s has no header !\n", filename);
    fflush(stderr);
  }

  /* have read header */

  /* allocate info block of function table */
  pt->maxsteps = 0;
  pt->ncols = ncols;
  pt->begin = (double *)malloc(ncols * sizeof(double));
  pt->end = (double *)malloc(ncols * sizeof(double));
  pt->step = (double *)malloc(ncols * sizeof(double));
  pt->invstep = (double *)malloc(ncols * sizeof(double));
  pt->len = (int *)malloc(ncols * sizeof(int));
  if ((pt->begin == NULL) || (pt->end == NULL) || (pt->step == NULL) || (pt->invstep == NULL)
    || (pt->len == NULL)) {
    sprintf(msg, "Cannot allocate info block for function table %s.", filename);
    KIM_API_report_error(__LINE__, __FILE__, msg, KIM_STATUS_FAIL);
    exit(EXIT_FAILURE);
  }

  /* catch the case where potential is identically zero */
  for (i = 0; i < ncols; ++i) {
    pt->end[i] = 0.0;
    pt->len[i] = 0;
  }

  /* read table */
  if (format == 1)
    read_pot_table1(pt, ncols, ntypes, filename, infile, radial);
  if (format == 2)
    read_pot_table2(pt, ncols, ntypes, filename, infile, radial);
  fclose(infile);

  init_threepoint(pt, ncols);

  return;
}

/****************************************************************
 *
 * read potential in first format: each line contains
 *
 * r**2 V00 V01 V02 ... V10 V11 V12 ... VNN
 *
 * N is the number of different atom types
 *
 * Note that it is assumed that the r**2 are aequidistant.
 *
 ****************************************************************/

void read_pot_table1(pot_table_t *pt, int ncols, int ntypes, char *filename, FILE *infile, int radial)
{
  char  msg[255];
  int   i, k;
  int   tablesize, npot = 0;
  double val, delta;
  double r2, r2_start = 0.0, r2_step;

  /* allocate the function table */
  pt->maxsteps = PSTEP;
  tablesize = ncols * pt->maxsteps;
  pt->table = (double *)malloc(tablesize * sizeof(double));
  if (NULL == pt->table) {
    sprintf(msg, "Cannot allocate memory for function table %s.", filename);
    KIM_API_report_error(__LINE__, __FILE__, msg, KIM_STATUS_FAIL);
    exit(EXIT_FAILURE);
  }

  /* input loop */
  while (!feof(infile)) {

    /* still some space left? */
    if (((npot % PSTEP) == 0) && (npot > 0)) {
      pt->maxsteps += PSTEP;
      tablesize = ncols * pt->maxsteps;
      pt->table = (double *)realloc(pt->table, tablesize * sizeof(double));
      if (NULL == pt->table) {
	sprintf(msg, "Cannot extend memory for function table %s.", filename);
	KIM_API_report_error(__LINE__, __FILE__, msg, KIM_STATUS_FAIL);
	exit(EXIT_FAILURE);
      }
    }

    /*  read in potential */
    if (1 != fscanf(infile, "%lf", &r2))
      break;
    if (npot == 0)
      r2_start = r2;		/* catch first value */
    for (i = 0; i < ncols; ++i) {
      if (1 != fscanf(infile, "%lf", &val)) {
	KIM_API_report_error(__LINE__, __FILE__, "Line incomplete in potential file.", KIM_STATUS_FAIL);
	exit(EXIT_FAILURE);
      }
      *PTR_2D(pt->table, npot, i, pt->maxsteps, ncols) = val;
      if (val != 0.) {		/* catch last non-zero value */
	pt->end[i] = r2;
	pt->len[i] = npot + 1;
      }
    }
    ++npot;
  }

  r2_step = (r2 - r2_start) / (npot - 1);

  /* fill info block, and shift potential to zero */
  for (i = 0; i < ncols; ++i) {
    pt->begin[i] = r2_start;
    pt->step[i] = r2_step;
    pt->invstep[i] = 1.0 / r2_step;
    delta = *PTR_2D(pt->table, (npot - 1), i, pt->maxsteps, ncols);
    /* if function of r2, shift potential and adjust cellsz */
    if (radial) {
      if (delta != 0.) {
	printf("Potential %1d%1d shifted by %f\n", (i / ntypes), (i % ntypes), delta);
	for (k = 0; k < npot; ++k)
	  *PTR_2D(pt->table, k, i, pt->table, ncols) -= delta;
      }
    }
  }

  /* increase table size for security */
  tablesize = ncols * (pt->maxsteps + 2);
  pt->table = (double *)realloc(pt->table, tablesize * sizeof(double));
  if (NULL == pt->table) {
    sprintf(msg, "Cannot extend memory for function table %s.", filename);
    KIM_API_report_error(__LINE__, __FILE__, msg, KIM_STATUS_FAIL);
    exit(EXIT_FAILURE);
  }

  return;
}

/****************************************************************
 *
 *  read potential in second format: at the beginning <ncols> times
 *  a line of the form
 *
 *  r_begin r_end r_step,
 *
 *  then the values of the potential (one per line), first those
 *  for atom pair  00, then an empty line (for gnuplot), then 01 and so on.
 *  Analogously, if there is only one column per atom type.
 *
 *  Note that it is assumed that the r**2 are aequidistant.
 *
 ****************************************************************/

void read_pot_table2(pot_table_t *pt, int ncols, int ntypes, char *filename, FILE *infile, int radial)
{
  char  msg[255];
  int   i, k;
  int   tablesize;
  double val, numstep, delta;

  /* read the info block of the function table */
  for (i = 0; i < ncols; i++) {
    if (3 != fscanf(infile, "%lf %lf %lf", &pt->begin[i], &pt->end[i], &pt->step[i])) {
      sprintf(msg, "Info line %d in %s corrupt.", i + 1, filename);
      KIM_API_report_error(__LINE__, __FILE__, msg, KIM_STATUS_FAIL);
      exit(EXIT_FAILURE);
    }
    pt->invstep[i] = 1.0 / pt->step[i];
    numstep = 1 + (pt->end[i] - pt->begin[i]) / pt->step[i];
    pt->len[i] = (int)(numstep + 0.49);
    pt->maxsteps = MAX(pt->maxsteps, pt->len[i]);

    /* some security against rounding errors */
    if (fabs(pt->len[i] - numstep) >= 0.1) {
      fprintf(stderr, "numstep = %f rounded to %d in file %s.\n", numstep, pt->len[i], filename);
      fflush(stderr);
    }
  }

  /* allocate the function table */
  /* allow some extra values at the end for interpolation */
  tablesize = ncols * (pt->maxsteps + 2);
  pt->table = (double *)malloc(tablesize * sizeof(double));
  if (NULL == pt->table) {
    sprintf(msg, "Cannot allocate memory for function table %s.", filename);
    KIM_API_report_error(__LINE__, __FILE__, msg, KIM_STATUS_FAIL);
    exit(EXIT_FAILURE);
  }

  /* input loop */
  for (i = 0; i < ncols; i++) {
    for (k = 0; k < pt->len[i]; k++) {
      if (1 != fscanf(infile, "%lf", &val)) {
	sprintf(msg, "wrong format in file %s.", filename);
	KIM_API_report_error(__LINE__, __FILE__, msg, KIM_STATUS_FAIL);
	exit(EXIT_FAILURE);
      }
      *PTR_2D(pt->table, k, i, pt->maxsteps, ncols) = val;
    }
  }

  /* if function of r2, shift potential if necessary */
  if (radial) {
    for (i = 0; i < ncols; i++) {
      delta = *PTR_2D(pt->table, pt->len[i] - 1, i, pt->maxsteps, ncols);
      if (delta != 0.0) {
	printf("Potential %1d%1d shifted by %f\n", (i / ntypes), (i % ntypes), delta);
	for (k = 0; k < pt->len[i]; k++)
	  *PTR_2D(pt->table, k, i, pt->table, ncols) -= delta;
      }
    }
  }

  return;
}

/****************************************************************
 *
 *  init_threepoint -- initialize for 3point interpolation
 *
 ****************************************************************/

void init_threepoint(pot_table_t *pt, int ncols)
{
  int   col, n;
  double *y;

  /* loop over columns */
  for (col = 0; col < ncols; col++) {

    y = pt->table + col;
    n = pt->len[col];

    /* for security, we continue the last interpolation polynomial */
    y[n * ncols] = 3 * y[(n - 1) * ncols] - 3 * y[(n - 2) * ncols] + y[(n - 3) * ncols];
    y[(n + 1) * ncols] = 6 * y[(n - 1) * ncols] - 8 * y[(n - 2) * ncols] + 3 * y[(n - 3) * ncols];

  }

  return;
}

/****************************************************************
 *
 * model driver destroy function
 * this is called by the KIM API after the calculation is done
 *
 ****************************************************************/

int destroy(void *km)
{
  intptr_t *pkim = *((intptr_t **) km);
  model_buffer *buffer;
  int ier;

  /* get model buffer from KIM object */
  buffer = (model_buffer *) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return KIM_STATUS_FAIL;
  }

  /* destroy all variables we allocated */
  free(buffer->dF_val);
  free(buffer->rho_val);

  /* destroy the buffer */
  free(buffer);

  return KIM_STATUS_OK;
}

/****************************************************************
 *
 * model driver reinitialization function
 * this function is called by the KIM API after the potentials have changed
 *
 * for IMD potentials this does not make sense, because we have no
 * published parameters which could be changed, so we provide a
 * dummy function which doesn't do anything
 *
 ****************************************************************/

int reinit(void *km)
{
  return KIM_STATUS_OK;
}

/****************************************************************
 *
 * warning function for short distances
 * this is called every time a short distance is encountered
 *
 ****************************************************************/

void short_dist_warning(int type, int i, int j, int *types, double *Rij, double R)
{
  int   l;

  if (0 == type)
    fprintf(stderr, "Short distance in the pair potential!\n");
  else if (1 == type)
    fprintf(stderr, "Short distance in the transfer function!\n");
  else
    fprintf(stderr, "Short distance in the embedding function!\n");
  fprintf(stderr, "Involved particles are %d (type %d) and %d (type %d).\n", i, types[i], j, types[j]);
  fprintf(stderr, "Relative position vector is (");
  for (l = 0; l < (DIM - 1); l++)
    fprintf(stderr, "%f,", Rij[l]);
  fprintf(stderr, "%f), distance is %f\n\n", Rij[DIM - 1], R);
  fflush(stderr);
}
