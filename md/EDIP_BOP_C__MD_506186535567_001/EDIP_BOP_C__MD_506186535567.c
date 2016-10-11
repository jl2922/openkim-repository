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
* Author: Daniel S. Karls (implementation adapted from Martin Z. Bazant's original code)
*
* Copyright (c) 2012, Regents of the University of Minnesota.  All rights reserved.
*
*/

/*******************************************************************************
*
*  model_driver_Si_BOP_C
*
*  EDIP (Environmentally-Dependent Interatomic Potential) for Silicon
*
*  Language: C
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
#define MAX_NBRS 1024
#define DIM 3 /* dimensionality of the space */

/* A note on how this model is implemented here:

I adapted this model from Prof. Martin Z. Bazant's (bazant@math.mit.edu) code (thanks also to Prof. Joao F. Justo for an additional FORTRAN reference code).
Because of the coordination-dependency of the 2-body and 3-body terms, we're forced to do multiple loops over all of the atoms.
The first inner loop (the "Pre-pass loop") goes over all of the other atoms and determines which ones are inside the cutoff radius of the atom being considered by the outer loop.
For those within the cutoff radius, it computes part of the 2-body interaction (the part that doesn't require the coordination to be known)
and the radial parts of the 3-body interaction.  Finally, the coordination of each of the atoms within the cutoff
is added onto the total coordination, Z, for that atom.  Having Z, we now have to loop over all other atoms once more
in order to calculate the coordination-dependent portion of the 2-body energy and get the final 2-body energy for the atom
being considered in the outermost loop, V2.  Note also that we are unable to perform a "half-summation" for the two-body terms
(with an outer loop over i and an inner loop over j>i) because of the asymmetry the coordination introduces into the 2-body energy (V2(R_ij,Z_i) != V2(R_ji,Z_j)).
We also need to use Z to perform another set of nested loops in order to calculate the coordination-dependent portion of the
3-body energy, h(l,Z).  This will give us the total 3-body energy for the atom being considered in the outermost loop, V3.
Since we are forced to loop over the same atoms multiple times (first to get the coordination and then to compute the energies),
it makes sense to reduce the computational expense of the subsequent looping by only looping over the atoms which have already been
determined to fall within the cutoff radius.  To this end, this code generates its own *internal* neighbor list even if it is not provided
by the test i.e. the test only supports the CLUSTER NBC mode.

The original EDIP publications can be found at:

     1.  M. Z. Bazant and E. Kaxiras, Phys. Rev. Lett. 77, 4370 (1996).
     2.  M. Z. Bazant, E. Kaxiras, J. F. Justo, Phys. Rev. B 56, 8542 (1997).
     3.  J. F. Justo, M. Z. Bazant, E. Kaxiras, V. V. Bulatov, and S. Yip,
	   Phys. Rev. B 58, 2539 (1998).

If you are interested in an explanation of the looping structure, etc., in this code and why it was chosen, I highly recommend looking at chapter 6 of Prof. Bazant's 1997 Ph.D. thesis, available at http://web.mit.edu/bazant/www/thesis/md.html.
Page 156 of this thesis includes a pseudocode of what has been implemented here.

Written by Daniel S. Karls (University of Minnesota). Contact: karl0100@umn.edu.
*/

/* ######## Function Prototypes ######## */

/* Define prototypes for model init */
int model_driver_init(void* km, char* paramfile_names, int* nmstrlen, int* numparamfiles);

/* Define prototypes: defined as static to avoid namespace clashes with other Models */
static int reinit(void* km);
static int destroy(void* km);
static int compute(void* km);

/* ######## Function Definitions ######## /*

/* Primary compute function*/
static int compute(void* km)
{
intptr_t* pkim = *((intptr_t**) km);

int* nAtoms;
int* particleSpecies;
double* cutoff;
double* coords;
double* energy;
double* force;
double* virial;
double* Rij_list;
int* neighListOfCurrentAtom;
int numOfAtomNeigh;
const char* NBCstr;
double* boxSideLengths;
int currentAtom;
int ier;
int i;
int NBC; /* Neighbor list flag */
int IterOrLoca; /* Neighbor list handling flag */
int comp_energy, comp_force, comp_virial;
/* Parameters...a is 2-body interaction cutoff....b is 3-body interaction cutoff (taken to be equal to 'a' in this code)...u1 - u4 are parameters for the interpolation function tau(z) parameters (Ismail & Kaxiras, 1993) */
double A;
double B;
double rh;
double a;
double sig;
double lam;
double gam;
double b;
double c;
double mu;
double Qo;
double eta;
double bet;
double alp;
double u1;
double u2;
double u3;
double u4;
double* param_A;
double* param_B;
double* param_rh;
double* param_a;
double* param_sig;
double* param_lam;
double* param_gam;
double* param_b;
double* param_c;
double* param_mu;
double* param_Qo;
double* param_eta;
double* param_bet;
double* param_alp;
double* param_u1;
double* param_u2;
double* param_u3;
double* param_u4;

   /* Determine neighbor list boundary condition (NBC) */
   ier = KIM_API_get_NBC_method(pkim, &NBCstr);
   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_NBC_method", ier);
      return ier;
   }
   if (!strcmp("CLUSTER",NBCstr))
   {
      NBC = 0;
   }
   else if (!strcmp("NEIGH_PURE_F",NBCstr))
   {
      NBC = 1;
   }
   else if (!strcmp("MI_OPBC_F",NBCstr))
   {
      NBC = 2;
   }
   else
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "Unsupported NBC method", ier);
      return ier;
   }

   /* Check to see if we have been asked to compute the energy or forces */
   KIM_API_getm_compute(pkim, &ier, 3*3,
                        "energy",         &comp_energy,         1,
                        "forces",         &comp_force,          1,
			"virial",         &comp_virial,         1);

   if (KIM_STATUS_OK > ier)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute", ier);
      return ier;
   }

   /* If the user asks for the virial or particleVirial but not the forces, then we need to allocate the space for (and compute)
   the forces anyway since they are required for the computation of the virial */
   if(comp_virial)
   {
     comp_force=1;
   }

   KIM_API_getm_data(pkim, &ier, 8*3,
                     "cutoff",                      &cutoff,         1,
                     "numberOfParticles",           &nAtoms,         1,
                     "particleSpecies",               &particleSpecies,  1,
                     "coordinates",                 &coords,         1,
                     "energy",                      &energy,         (comp_energy==1),
		       "forces",			   &force,          (comp_force==1),
		       "virial",                      &virial,         (comp_virial==1),
		       "boxSideLengths",              &boxSideLengths, (NBC==2));

   if (KIM_STATUS_OK > ier){
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", ier);
      return ier;}

   /* Get model parameters from KIM API object */
   KIM_API_getm_data(pkim, &ier, 18*3,
		"PARAM_FREE_cutoff", &param_a, 1,
		"PARAM_FREE_A", &param_A, 1,
		"PARAM_FREE_B", &param_B, 1,
		"PARAM_FREE_rh", &param_rh, 1,
		"PARAM_FREE_sig", &param_sig, 1,
		"PARAM_FREE_lam", &param_lam, 1,
		"PARAM_FREE_gam", &param_gam, 1,
		"PARAM_FREE_b", &param_b, 1,
		"PARAM_FREE_c", &param_c, 1,
		"PARAM_FREE_mu", &param_mu, 1,
		"PARAM_FREE_Qo", &param_Qo, 1,
		"PARAM_FREE_eta", &param_eta, 1,
		"PARAM_FREE_bet", &param_bet, 1,
		"PARAM_FREE_alp", &param_alp, 1,
		"PARAM_FREE_u1", &param_u1, 1,
		"PARAM_FREE_u2", &param_u2, 1,
		"PARAM_FREE_u3", &param_u3, 1,
		"PARAM_FREE_u4", &param_u4, 1
		);
   if (KIM_STATUS_OK > ier){
   KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", ier);
   return ier;}


/* Set up local variables */
/* ---------------------------- */
int q,j,k,jj,j2,j3,l,nk,nl,neighCounter;
int num2[MAX_NBRS], num3[MAX_NBRS], numz[MAX_NBRS];
a = *param_a;
A = *param_A;
B = *param_B;
rh = *param_rh;
sig = *param_sig;
lam = *param_lam;
gam = *param_gam;
b = *param_b;
c = *param_c;
mu = *param_mu;
Qo = *param_Qo;
eta = *param_eta;
bet = *param_bet;
alp = *param_alp;
u1 = *param_u1;
u2 = *param_u2;
u3 = *param_u3;
u4 = *param_u4;
double bg=a; /*  Cutoff for radial part of 3-body interactions, g(r) */

typedef struct
{double t0, t1, t2, t3;
 double dx, dy, dz;
 double r;
} store2;

typedef struct
{
    double g,dg;         /* 3-body radial function and its derivative */
    double rinv;         /* 1/r */
    double dx,dy,dz;     /* unit separation vector */
    double r;		 /* bond length (only needed for virial) */
} store3;

typedef struct
{
    double df;
    double dx,dy,dz;
    double r;
} storez;

store2 s2[MAX_NBRS]={0.0};  /* Initialize to all 0.0 */
store3 s3[MAX_NBRS]={0.0};
storez sz[MAX_NBRS]={0.0};

int n2, n3, nz, dummycount;
double V2, V3;
double dx,dy,dz,r,rsqr,asqr;
double rinv,rmainv,xinv,xinv3,den,Z,fZ;
double dV2j,dV2ijx,dV2ijy,dV2ijz,pZ,dp;
double temp0,temp1,temp2;
double Qort,muhalf,u5;
double rmbinv,winv,dwinv,tau,dtau,lcos,x,H, dHdx, dhdl;
double dV3rij,dV3rijx,dV3rijy,dV3rijz;
double dV3rik,dV3rikx,dV3riky,dV3rikz;
double dV3l,dV3ljx,dV3ljy,dV3ljz,dV3lkx,dV3lky,dV3lkz;
double dV2dZ,dxdZ,dV3dZ,dVdZ_sum;
double dEdrl,dEdrlx,dEdrly,dEdrlz;
double bmc,cmbinv;
double fjx,fjy,fjz,fkx,fky,fkz;

/* Combine Coefficients */
asqr=a*a;
Qort = sqrt(Qo);
muhalf = mu*0.5;
u5 = u2*u4;
bmc = b-c;
cmbinv = 1.0/(c-b);

 /* Determine neighbor list handling mode */
   if (NBC != 0)
   {
      /*****************************
       * IterOrLoca = 1 -- Iterator      Right now, this model only supports LOCATOR mode!
       *            = 2 -- Locator
       *            = 3 -- Both
       *****************************/
      IterOrLoca = KIM_API_get_neigh_mode(pkim, &ier);

      if (KIM_STATUS_OK > ier)
      {
         KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh_mode", ier);
         return ier;
      }
      if (IterOrLoca != 2)
      {  KIM_API_print(pkim,&ier);
         printf("* ERROR: Unsupported IterOrLoca mode = %i\n", IterOrLoca);
         return KIM_STATUS_FAIL;
      }
   }
   else
   {
      IterOrLoca = 99;   /* for CLUSTER NBC...the value of 99 here isn't actually used in the code */
   }

   /* Initialize neighbor handling for CLUSTER NBC */
   if (0 == NBC) /* CLUSTER */
   {
      neighListOfCurrentAtom = (int *) malloc((*nAtoms)*sizeof(int));
   }

/* -----------------------------Start actual computations--------------------------------------------- */

	/* Initialize Forces and Energies to zero */
	if(comp_force)
        {
	  for(i=0;i<*nAtoms;i++)
	  {
	  force[i*DIM+0]=0.0;
	  force[i*DIM+1]=0.0;
	  force[i*DIM+2]=0.0;
	  }
	}

	if(comp_virial)
	{
	  for(i=0;i<6;i++)
	  {virial[i]=0.0;
	  }
	}

	V2=0.0;
	V3=0.0;

	/* --- Level 1: Outer loop over atoms --- */
       for(i=0;i<*nAtoms;i++)
	{

      /* Set up neighbor list for next atom for all NBC methods */
         if (0 == NBC)     /* CLUSTER NBC method */
         {
	  numOfAtomNeigh=*nAtoms-1;

	      /* For the CLUSTER NBC method, we set the neighbor list equal to all of the atoms (notice that we check whether j == i in the loops below) */
	      neighCounter=0;
	      for (q=0; q < *nAtoms; q++)
	      {
	        if(q!=i)
	 	{neighListOfCurrentAtom[neighCounter]=q;
	      	neighCounter++;}
	      }
         }
         else              /* All other NBCs */
         {
            ier = KIM_API_get_neigh(pkim, 1, i, &currentAtom, &numOfAtomNeigh,
                                     &neighListOfCurrentAtom, &Rij_list);
            if (KIM_STATUS_OK != ier) /* some sort of problem, exit */
            {printf("i=%i\n",i);
	     printf("currentAtom = %i\n",currentAtom);
 	     printf("numOfAtomNeigh = %i\n",numOfAtomNeigh);
 	     printf("*nAtoms = %i\n",*nAtoms);
            KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", ier);
            ier = KIM_STATUS_FAIL;
            return ier;
            }
         }

		 /* RESET COORDINATION AND NEIGHBOR NUMBERS */

   		 Z = 0.0;
   		 n2 = 0;
   		 n3 = 0;
   		 nz = 0;
                 
		for(dummycount=0;dummycount<MAX_NBRS;dummycount++)
                {numz[dummycount]=0;}

		/* --- Level 2: Prepass loop over neighbors of atom i --- */
		for (jj = 0; jj < numOfAtomNeigh; ++ jj)
        	{
       		j = neighListOfCurrentAtom[jj]; /* get neighbor ID */

		if(j!=i)
		{
		dx=coords[j*DIM+0]-coords[i*DIM+0];
		dy=coords[j*DIM+1]-coords[i*DIM+1];
		dz=coords[j*DIM+2]-coords[i*DIM+2];

		/* apply periodic boundary conditions if required */
            	if (2 == NBC)
            	{
                 if (abs(dx) > 0.5*boxSideLengths[0])
                 dx -= (dx/fabs(dx))*boxSideLengths[0];
		   if (abs(dy) > 0.5*boxSideLengths[1])
		   dy -= (dy/fabs(dy))*boxSideLengths[1];
		   if (abs(dz) > 0.5*boxSideLengths[2])
		   dz -= (dz/fabs(dz))*boxSideLengths[2];
               }
			if(dx < a)
			{
				if(dy < a)
				{
					if(dz < a)
					{
					  rsqr=dx*dx+dy*dy+dz*dz;
					  if(rsqr < asqr)
					  {
					    num2[n2]=j;
					    r=sqrt(rsqr);
				   	    rinv=1.0/r;
					    dx *= rinv;
					    dy *= rinv;
					    dz *= rinv;
					    rmainv=1.0/(r-a);
    	                		    s2[n2].t0 = A*exp(sig*rmainv);
                        		    s2[n2].t1 = pow(B*rinv,rh);
  				           s2[n2].t2 = rh*rinv;
                        		    s2[n2].t3 = sig*rmainv*rmainv;
       				    s2[n2].dx = dx;
      				           s2[n2].dy = dy;
			                  s2[n2].dz = dz;
 					    s2[n2].r = r;
					    n2++;

					    /* Also compute the radial part of the 3-body energy term */
					    if(r<bg)
					    {
					      num3[n3]=j;
					      rmbinv = 1.0/(r-bg);
					      temp1 = gam*rmbinv;
                          		      temp0 = exp(temp1);
					      s3[n3].g = temp0;
	 				      s3[n3].dg = -rmbinv*temp1*temp0;
                          		      s3[n3].dx = dx;
                         		      s3[n3].dy = dy;
                          		      s3[n3].dz = dz;
                          		      s3[n3].rinv = rinv;
                          		      s3[n3].r = r;
                          		      n3++;

							/* Coordination and neighbor function c < r < a */
							if(r<b)
							{
							  if(r<=c)
							  Z+=1.0;
							  else
                             			          {
	   						  xinv = bmc/(r-c);
    						   	  xinv3 = xinv*xinv*xinv;
	                           			  den = 1.0/(1 - xinv3);
	 						  temp1 = alp*den;
	  						  fZ = exp(temp1);
   						         Z += fZ;
							  numz[nz]=j;
 						         sz[nz].df = fZ*temp1*den*3.0*xinv3*xinv*cmbinv;   /* df/dr */
   						         sz[nz].dx = dx;
   						         sz[nz].dy = dy;
   						         sz[nz].dz = dz;
					           	  sz[nz].r = r;
							  nz++;
							  }
							} /* End if-statement of r<b */

					    } /* End if-statement r<bg */

					   }/*  End if-statement rsqrt < asqr */

					 }
				 }
			 }

		} /* End if-statement j!=i */


		} /*  End loop over j */

		dVdZ_sum=0.0;
        	temp0 = bet*Z;
   		pZ = exp(-temp0*Z); 
    		dp = -2.0*temp0*pZ;         /* derivative of bond order */

		/* --- Level 2: Second loop over pairs to get the final 2-body energy --- */
		for(j2=0; j2<n2; j2++)
		{
      		  temp0 = s2[j2].t1 - pZ;

		  /*  Two-body energy */
		  V2 += temp0*s2[j2].t0;

      		  /* two-body forces */
	          if(comp_force)
      		  {
      		    dV2j = - (s2[j2].t0) * ((s2[j2].t1)*(s2[j2].t2) + temp0 * (s2[j2].t3));   /* dV2/dr */
      		    dV2ijx = dV2j * s2[j2].dx; /*  x component of force on atom i */
      		    dV2ijy = dV2j * s2[j2].dy; /*  y component of force on atom i */
      		    dV2ijz = dV2j * s2[j2].dz; /*  z component of force on atom i */
		    force[i*DIM+0] += dV2ijx;
      		    force[i*DIM+1] += dV2ijy;
      		    force[i*DIM+2] += dV2ijz;
		    j = num2[j2];
		    force[j*DIM+0] -= dV2ijx;
      		    force[j*DIM+1] -= dV2ijy;
      		    force[j*DIM+2] -= dV2ijz;

      		    /* accumulation of pair coordination forces */
      		    dV2dZ = - dp * s2[j2].t0;
      		    dVdZ_sum += dV2dZ;
		   }

		   /* Two-body contribution to virial */
		   if(comp_virial)
		   {
		     virial[0] += s2[j2].r * dV2ijx * s2[j2].dx;
		     virial[1] += s2[j2].r * dV2ijy * s2[j2].dy;
		     virial[2] += s2[j2].r * dV2ijz * s2[j2].dz;
		     virial[3] += s2[j2].r * dV2ijy * s2[j2].dz;
		     virial[4] += s2[j2].r * dV2ijz * s2[j2].dx;
		     virial[5] += s2[j2].r * dV2ijx * s2[j2].dy;
		   }

		 } /* End loop over j2 */


		/* Coordination-dependence of three-body interactions */
    		winv = Qort*exp(-muhalf*Z); /* inverse width of angular function */
    		dwinv = -muhalf*winv;       /* its derivative */
    		temp0 = exp(-u4*Z);
    		tau = u1+u2*temp0*(u3-temp0); /* -cosine of angular minimum */
    		dtau = u5*temp0*(2*temp0-u3); /* its derivative */

		/* --- Level 2: First loop for three-body interactions --- */
		for(j3=0; j3<(n3-1); j3++)
		{j=num3[j3];

			/* --- Level 3: Second loop for three-body interactions --- */
			for(nk=j3+1; nk<n3; nk++)
			{
			  k=num3[nk];
			  /* angular function h(l,Z) */
			  lcos = s3[j3].dx * s3[nk].dx + s3[j3].dy * s3[nk].dy + s3[j3].dz * s3[nk].dz;
			  x = (lcos + tau)*winv;
        		  temp0 = exp(-x*x);

        		  H = lam*(1 - temp0 + eta*x*x);
        		  dHdx = 2*lam*x*(temp0 + eta);
			  dhdl = dHdx*winv;

			  /* three-body energy */
			  temp1 = s3[j3].g * s3[nk].g;
	 		  V3 += temp1*H;

			  if(comp_force)
		          {
			    /* (-) radial force on atom j */
			    dV3rij = s3[j3].dg * s3[nk].g * H;
			    dV3rijx = dV3rij * s3[j3].dx;
			    dV3rijy = dV3rij * s3[j3].dy;
			    dV3rijz = dV3rij * s3[j3].dz;
			    fjx = dV3rijx;
			    fjy = dV3rijy;
			    fjz = dV3rijz;

			    /* (-) radial force on atom k */
			    dV3rik = s3[j3].g * s3[nk].dg * H;
			    dV3rikx = dV3rik * s3[nk].dx;
			    dV3riky = dV3rik * s3[nk].dy;
			    dV3rikz = dV3rik * s3[nk].dz;
			    fkx = dV3rikx;
			    fky = dV3riky;
			    fkz = dV3rikz;

			    /* (-) angular force on j */
			    dV3l = temp1*dhdl;
			    dV3ljx = dV3l * (s3[nk].dx - lcos * s3[j3].dx) * s3[j3].rinv;
			    dV3ljy = dV3l * (s3[nk].dy - lcos * s3[j3].dy) * s3[j3].rinv;
			    dV3ljz = dV3l * (s3[nk].dz - lcos * s3[j3].dz) * s3[j3].rinv;
			    fjx += dV3ljx;
			    fjy += dV3ljy;
			    fjz += dV3ljz;

			    /* (-) angular force on k */
			    dV3lkx = dV3l * (s3[j3].dx - lcos * s3[nk].dx) * s3[nk].rinv;
			    dV3lky = dV3l * (s3[j3].dy - lcos * s3[nk].dy) * s3[nk].rinv;
			    dV3lkz = dV3l * (s3[j3].dz - lcos * s3[nk].dz) * s3[nk].rinv;
			    fkx += dV3lkx;
			    fky += dV3lky;
			    fkz += dV3lkz;

			    /* apply radial + angular forces to i, j, k */
			    force[j*DIM+0] -= fjx;
		           force[j*DIM+1] -= fjy;
		           force[j*DIM+2] -= fjz;
		           force[k*DIM+0] -= fkx;
		           force[k*DIM+1] -= fky;
		           force[k*DIM+2] -= fkz;
			    force[i*DIM+0] += fjx + fkx;
			    force[i*DIM+1] += fjy + fky;
			    force[i*DIM+2] += fjz + fkz;
			  }

			  /* dV3/dR contributions to virial */
			  if(comp_virial)
			  {
			    virial[0] += s3[j3].r * (fjx*s3[j3].dx) + s3[nk].r * (fkx*s3[nk].dx);
			    virial[1] += s3[j3].r * (fjy*s3[j3].dy) + s3[nk].r * (fky*s3[nk].dy);
			    virial[2] += s3[j3].r * (fjz*s3[j3].dz) + s3[nk].r * (fkz*s3[nk].dz);
			    virial[3] += s3[j3].r * (fjy*s3[j3].dz) + s3[nk].r * (fky*s3[nk].dz);
			    virial[4] += s3[j3].r * (fjz*s3[j3].dx) + s3[nk].r * (fkz*s3[nk].dx);
			    virial[5] += s3[j3].r * (fjx*s3[j3].dy) + s3[nk].r * (fkx*s3[nk].dy);
			  }

			  dxdZ = dwinv*(lcos + tau) + winv*dtau;
			  dV3dZ = temp1*dHdx*dxdZ;

			  /* accumulation of 3-body coordination forces */
        		  dVdZ_sum += dV3dZ;
			} /* End loop over k */

		} /* End loop over j3 */


		/* --- Level 2: Loop to apply coordination forces --- */
		if(comp_force)
		{
    		 for(nl=0; nl<nz; nl++)
		 {
        		dEdrl = dVdZ_sum * sz[nl].df;
        		dEdrlx = (dEdrl * sz[nl].dx); /* The EDIP currently in LAMMPS (as of the May 21, 2012 release) has the signs on dEdrlx, dEdrly, and dEdrlz reversed.  The signs used here are correct. */
        		dEdrly = (dEdrl * sz[nl].dy);
        		dEdrlz = (dEdrl * sz[nl].dz);
        		force[i*DIM+0] += dEdrlx;
        		force[i*DIM+1] += dEdrly;
        		force[i*DIM+2] += dEdrlz;
        		l = numz[nl];
        		force[l*DIM+0] -= dEdrlx;
        		force[l*DIM+1] -= dEdrly;
        		force[l*DIM+2] -= dEdrlz;

		   /* dE/dZ*dZ/dr contributions to virial */
		   if(comp_virial) /*  Recall that, above, we set comp_force to 1 if comp_virial is 1 so this nested if-statement is ok */
		   {
		     virial[0] += sz[nl].r * (dEdrlx * sz[nl].dx);
		     virial[1] += sz[nl].r * (dEdrly * sz[nl].dy);
		     virial[2] += sz[nl].r * (dEdrlz * sz[nl].dz);
		     virial[3] += sz[nl].r * (dEdrly * sz[nl].dz);
		     virial[4] += sz[nl].r * (dEdrlz * sz[nl].dx);
		     virial[5] += sz[nl].r * (dEdrlx * sz[nl].dy);
		   }
   		 }
		}


	} /* End loop over i */

if(comp_energy)
{*energy = V2+V3;}

   /* Free temporary storage */
   if (0 == NBC)
   {free(neighListOfCurrentAtom);}

   /* everything is great */
   ier = KIM_STATUS_OK;

return ier;}







/* Initialization function */
int model_driver_init(void* km, char* paramfile_names, int* nmstrlen, int* numparamfiles)
{
   /* KIM variables */
   intptr_t* pkim = *((intptr_t**) km);
   char* paramfile1name;

   /* Local variables */
   FILE* fid;
   double* model_Pcutoff;
   double* model_cutoff;
   double* model_A;
   double* model_B;
   double* model_rh;
   double* model_a;
   double* model_sig;
   double* model_lam;
   double* model_gam;
   double* model_b;
   double* model_c;
   double* model_mu;
   double* model_Qo;
   double* model_eta;
   double* model_bet;
   double* model_alp;
   double* model_u1;
   double* model_u2;
   double* model_u3;
   double* model_u4;
   double cutoff;
   double A, B, rh, a, sig, lam, gam, b, c, mu, Qo, eta, bet, alp, u1, u2, u3, u4;
   int ier;   
	
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
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_setm_data", ier);
      return ier;
   }


   /* Read in model parameters from parameter file */
   fid=fopen(paramfile1name,"r");
   if(fid==NULL)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "Unable to open parameter file for EDIP parameters", ier);
      return ier;
   }

   ier = fscanf(fid, "%lf \n%lf \n%lf \n%lf %lf \n%lf \n%lf \n%lf \n%lf \n%lf \n%lf %lf \n%lf \n%lf \n%lf %lf \n%lf \n%lf",
		&a, /* cutoff in Angstroms */
		&A, &B, &rh, &sig, &lam, &gam, &b, &c, &mu, &Qo, &eta, &bet, &alp,
		&u1, &u2, &u3, &u4);
   fclose(fid);

   /* check that we read the right number of parameters */
   if (18 != ier)
   {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "Unable to read all EDIP parameters", ier);
      printf("number of parameters read: %i\n",ier);
      return ier;
   }

   /* convert to appropriate units */
   A *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
 	                                  0.0, 1.0,  0.0, 0.0, 0.0, &ier);
 	   if (KIM_STATUS_OK > ier){
	      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
 	      return ier;}

   B *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
 	                                  1.0, 0.0,  0.0, 0.0, 0.0, &ier);
 	   if (KIM_STATUS_OK > ier){
	      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
 	      return ier;}

   /*  rh is unitless */

   a *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
 	                                  1.0, 0.0,  0.0, 0.0, 0.0, &ier);
 	   if (KIM_STATUS_OK > ier){
	      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
 	      return ier;}

   sig *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
 	                                  1.0, 0.0,  0.0, 0.0, 0.0, &ier);
 	   if (KIM_STATUS_OK > ier){
	      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
 	      return ier;}

   lam *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
 	                                  0.0, 1.0,  0.0, 0.0, 0.0, &ier);
 	   if (KIM_STATUS_OK > ier){
	      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
 	      return ier;}

   gam *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
 	                                  1.0, 0.0,  0.0, 0.0, 0.0, &ier);
 	   if (KIM_STATUS_OK > ier){
	      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
 	      return ier;}

   b *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
 	                                  1.0, 0.0,  0.0, 0.0, 0.0, &ier);
 	   if (KIM_STATUS_OK > ier){
	      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
 	      return ier;}

   c *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
 	                                  1.0, 0.0,  0.0, 0.0, 0.0, &ier);
 	   if (KIM_STATUS_OK > ier){
	      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
 	      return ier;}

   /*  mu, Qo, eta, bet, alp, u1, u2, u3, u4 are unitless */

   /* store model cutoff in KIM object */
   model_cutoff = (double*) KIM_API_get_data(pkim, "cutoff", &ier);
   if (KIM_STATUS_OK > ier){
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", ier);
      return ier;}
   *model_cutoff = a;

   /* allocate memory for parameters */
   model_Pcutoff = (double*) malloc(1*sizeof(double));
   if(NULL == model_Pcutoff){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}
   model_A = (double*) malloc(1*sizeof(double));
   if(NULL == model_A){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}
   model_B = (double*) malloc(1*sizeof(double));
   if(NULL == model_B){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}
   model_rh = (double*) malloc(1*sizeof(double));
   if(NULL == model_rh){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}
   model_sig = (double*) malloc(1*sizeof(double));
   if(NULL == model_sig){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}
   model_lam = (double*) malloc(1*sizeof(double));
   if(NULL == model_lam){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}
   model_gam = (double*) malloc(1*sizeof(double));
   if(NULL == model_gam){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}
   model_b = (double*) malloc(1*sizeof(double));
   if(NULL == model_b){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}
   model_c = (double*) malloc(1*sizeof(double));
   if(NULL == model_c){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}
   model_mu = (double*) malloc(1*sizeof(double));
   if(NULL == model_mu){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}
   model_Qo = (double*) malloc(1*sizeof(double));
   if(NULL == model_Qo){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}
   model_eta = (double*) malloc(1*sizeof(double));
   if(NULL == model_eta){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}
   model_bet = (double*) malloc(1*sizeof(double));
   if(NULL == model_bet){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}
   model_alp = (double*) malloc(1*sizeof(double));
   if(NULL == model_alp){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}
   model_u1 = (double*) malloc(1*sizeof(double));
   if(NULL == model_u1){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}
   model_u2 = (double*) malloc(1*sizeof(double));
   if(NULL == model_u2){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}
   model_u3 = (double*) malloc(1*sizeof(double));
   if(NULL == model_u3){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}
   model_u4 = (double*) malloc(1*sizeof(double));
   if(NULL == model_u4){
   	KIM_API_report_error(__LINE__, __FILE__, "malloc", KIM_STATUS_FAIL);
	return KIM_STATUS_FAIL;}

   /* store parameters in KIM object */
   KIM_API_setm_data(pkim, &ier, 18*4,
		"PARAM_FREE_cutoff", 1, model_Pcutoff, 1,
		"PARAM_FREE_A", 1, model_A, 1,
		"PARAM_FREE_B", 1, model_B, 1,
		"PARAM_FREE_rh", 1, model_rh, 1,
		"PARAM_FREE_sig", 1, model_sig, 1,
		"PARAM_FREE_lam", 1, model_lam, 1,
		"PARAM_FREE_gam", 1, model_gam, 1,
		"PARAM_FREE_b", 1, model_b, 1,
		"PARAM_FREE_c", 1, model_c, 1,
		"PARAM_FREE_mu", 1, model_mu, 1,
		"PARAM_FREE_Qo", 1, model_Qo, 1,
		"PARAM_FREE_eta", 1, model_eta, 1,
                "PARAM_FREE_bet", 1, model_bet, 1,
		"PARAM_FREE_alp", 1, model_alp, 1,
		"PARAM_FREE_u1", 1, model_u1, 1,
		"PARAM_FREE_u2", 1, model_u2, 1,
		"PARAM_FREE_u3", 1, model_u3, 1,
		"PARAM_FREE_u4", 1, model_u4, 1
		);

   /* set value of parameters */
   *model_Pcutoff = *model_cutoff;
   *model_A = A;
   *model_B = B;
   *model_rh = rh;
   *model_sig = sig;
   *model_lam = lam;
   *model_gam = gam;
   *model_b = b;
   *model_c = c;
   *model_mu = mu;
   *model_Qo = Qo;
   *model_eta = eta;
   *model_bet = bet;
   *model_alp = alp;
   *model_u1 = u1;
   *model_u2 = u2;
   *model_u3 = u3;
   *model_u4 = u4;

   ier = KIM_STATUS_OK;
return ier;
}







/* Reinitialization function */
static int reinit(void *km)
{
   /* Local variables */
   intptr_t* pkim = *((intptr_t**) km);
   int ier;
   double* cutoff;
   double* param_cutoff;
   double* param_A;
   double* param_B;
   double* param_rh;
   double* param_a;
   double* param_sig;
   double* param_lam;
   double* param_gam;
   double* param_b;
   double* param_c;
   double* param_mu;
   double* param_Qo;
   double* param_eta;
   double* param_bet;
   double* param_alp;
   double* param_u1;
   double* param_u2;
   double* param_u3;
   double* param_u4;

   /* Get parameters from KIM object */
   KIM_API_getm_data(pkim, &ier, 18*3,
		"PARAM_FREE_cutoff", &param_cutoff, 1,
		"PARAM_FREE_A", &param_A, 1,
		"PARAM_FREE_B", &param_B, 1,
		"PARAM_FREE_rh", &param_rh, 1,
		"PARAM_FREE_sig", &param_sig, 1,
		"PARAM_FREE_lam", &param_lam, 1,
		"PARAM_FREE_gam", &param_gam, 1,
		"PARAM_FREE_b", &param_b, 1,
		"PARAM_FREE_c", &param_c, 1,
		"PARAM_FREE_mu", &param_mu, 1,
		"PARAM_FREE_Qo", &param_Qo, 1,
		"PARAM_FREE_eta", &param_eta, 1,
                "PARAM_FREE_bet", &param_bet, 1,
		"PARAM_FREE_alp", &param_alp, 1,
		"PARAM_FREE_u1", &param_u1, 1,
		"PARAM_FREE_u2", &param_u2, 1,
		"PARAM_FREE_u3", &param_u3, 1,
		"PARAM_FREE_u4", &param_u4, 1
		);
   if (KIM_STATUS_OK > ier){
   KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", ier);
   return ier;}

   /* set new values in KIM object     */
   /*                                  */
   /* store model cutoff in KIM object */
   cutoff = KIM_API_get_data(pkim, "cutoff", &ier);
   *cutoff = *param_cutoff;

   ier = KIM_STATUS_OK;
return ier;
}







static int destroy(void *km)
{
   /* Local variables */
   intptr_t* pkim = *((intptr_t**) km);
   int ier;
   double* param_cutoff;
   double* param_A;
   double* param_B;
   double* param_rh;
   double* param_a;
   double* param_sig;
   double* param_lam;
   double* param_gam;
   double* param_b;
   double* param_c;
   double* param_mu;
   double* param_Qo;
   double* param_eta;
   double* param_bet;
   double* param_alp;
   double* param_u1;
   double* param_u2;
   double* param_u3;
   double* param_u4;

   /* Get parameter addresses from KIM object */
   KIM_API_getm_data(pkim, &ier, 18*3,
		"PARAM_FREE_cutoff", &param_cutoff, 1,
		"PARAM_FREE_A", &param_A, 1,
		"PARAM_FREE_B", &param_B, 1,
		"PARAM_FREE_rh", &param_rh, 1,
		"PARAM_FREE_sig", &param_sig, 1,
		"PARAM_FREE_lam", &param_lam, 1,
		"PARAM_FREE_gam", &param_gam, 1,
		"PARAM_FREE_b", &param_b, 1,
		"PARAM_FREE_c", &param_c, 1,
		"PARAM_FREE_mu", &param_mu, 1,
		"PARAM_FREE_Qo", &param_Qo, 1,
		"PARAM_FREE_eta", &param_eta, 1,
                "PARAM_FREE_bet", &param_bet, 1,
		"PARAM_FREE_alp", &param_alp, 1,
		"PARAM_FREE_u1", &param_u1, 1,
		"PARAM_FREE_u2", &param_u2, 1,
		"PARAM_FREE_u3", &param_u3, 1,
		"PARAM_FREE_u4", &param_u4, 1
		);
   if (KIM_STATUS_OK > ier){
   KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", ier);
   return ier;}



   /* Free the memory at the parameter addresses */
   free(param_cutoff);
   free(param_A);
   free(param_B);
   free(param_rh);
   free(param_sig);
   free(param_lam);
   free(param_gam);
   free(param_b);
   free(param_c);
   free(param_mu);
   free(param_Qo);
   free(param_eta);
   free(param_bet);
   free(param_alp);
   free(param_u1);
   free(param_u2);
   free(param_u3);
   free(param_u4);

   ier = KIM_STATUS_OK;
return ier;
}
