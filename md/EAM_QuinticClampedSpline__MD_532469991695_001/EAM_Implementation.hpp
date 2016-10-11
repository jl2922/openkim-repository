//
// CDDL HEADER START
//
// The contents of this file are subject to the terms of the Common Development
// and Distribution License Version 1.0 (the "License").
//
// You can obtain a copy of the license at
// http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
// specific language governing permissions and limitations under the License.
//
// When distributing Covered Code, include this CDDL HEADER in each file and
// include the License file in a prominent location with the name LICENSE.CDDL.
// If applicable, add the following below this CDDL HEADER, with the fields
// enclosed by brackets "[]" replaced with your own identifying information:
//
// Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
//
// CDDL HEADER END
//

//
// Copyright (c) 2013--2015, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Ryan S. Elliott
//    Stephen M. Whalen
//


#ifndef EAM_IMPLEMENTATION_HPP_
#define EAM_IMPLEMENTATION_HPP_

#include <map>
#include <cmath>
#include "EAM.hpp"
#include "KIM_API_status.h"

#define MAXLINE 1024
#define DIMENSION 3
#define ONE 1.0
#define TWO 2.0
#define HALF 0.5

#define MAX_PARAMETER_FILES 20
#define NUMBER_SETFL_COMMENT_LINES 3

#include "EAM_Spline.hpp"

//==============================================================================
//
// Type definitions, enumerations, and helper function prototypes
//
//==============================================================================

// type declaration for get neighbor functions
typedef int (GetNeighborFunction)(void**, int*, int*, int*, int*, int**,
                                  double**);
// type declaration for vector of constant dimension
typedef double VectorOfSizeDIM[DIMENSION];
// type declaration for funcfl data
struct SetOfFuncflData
{
  int numberRhoPoints[MAX_PARAMETER_FILES];
  double deltaRho[MAX_PARAMETER_FILES];
  int numberRPoints[MAX_PARAMETER_FILES];
  double deltaR[MAX_PARAMETER_FILES];
  double cutoff[MAX_PARAMETER_FILES];
  double* embeddingData[MAX_PARAMETER_FILES];
  double* densityData[MAX_PARAMETER_FILES];
  double* ZData[MAX_PARAMETER_FILES];
};
// type declaration for deferred neighbors
typedef struct
{
  int index;
  VectorOfSizeDIM rij;
} neighbor;
// type declaration for iterating over deferred neighbor lists
typedef std::multimap<int, neighbor>::iterator deferredNeighborIterator;

// enumeration for the different types of NBC's
enum NBCTypeEnum {Neigh_Rvec, Neigh_Pure, Mi_Opbc, Cluster};
// enumeration for the different methods of computing rij
enum RijEnum {Coordinates, RVec, MI_OPBC};
// enumeration for EAMFileType
enum EAMFileType {Setfl, Funcfl, FinnisSinclair, Error};

//==============================================================================
//
// Helper class definitions
//
//==============================================================================

// Iterator object for CLUSTER NBC
class ClusterIterator
{
 private:
  int* list_;
  int baseconvert_;
  int const cachedNumberContributingParticles_;
  int request_;
 public:
  ClusterIterator(KIM_API_model* const pkim,
                  GetNeighborFunction* const get_neigh,
                  int const baseconvert,
                  int const cachedNumberContributingParticles,
                  int* const i,
                  int* const numnei,
                  int** const n1atom,
                  double** const pRij)
      : baseconvert_(baseconvert),
        cachedNumberContributingParticles_(cachedNumberContributingParticles),
        request_(0)
  {
    // allocate memory for neighbor list
    list_ = new int[cachedNumberContributingParticles_];
    for (int k = 0; k < cachedNumberContributingParticles_; ++k)
      list_[k] = k - baseconvert_;

    *i = request_;
    // CLUSTER always uses half-list behavior
    *numnei = cachedNumberContributingParticles_ - request_ - 1;
    *n1atom = &(list_[request_ + 1]);
    *pRij = NULL;
  }
  ~ClusterIterator()
  {
    delete [] list_;
  }
  bool done() const
  {
    return !(request_ < cachedNumberContributingParticles_);
  }
  int next(int* const i, int* const numnei, int** const n1atom,
           double** const pRij)
  {
    ++request_;

    *i = request_;
    *numnei = cachedNumberContributingParticles_ - request_ - 1;
    *n1atom = &(list_[request_ + 1]);
    *pRij = NULL;
    return KIM_STATUS_OK;
  }
};

// Iterator object for Locator mode access to neighbor list
class LocatorIterator
{
 private:
  KIM_API_model* const pkim_;
  GetNeighborFunction* const get_neigh_;
  int const baseconvert_;
  int const cachedNumberContributingParticles_;
  int request_;
  int const mode_;
 public:
  LocatorIterator(KIM_API_model* const pkim,
                  GetNeighborFunction* const get_neigh,
                  int const baseconvert,
                  int const cachedNumberContributingParticles,
                  int* const i,
                  int* const numnei,
                  int** const n1atom,
                  double** const pRij)
      : pkim_(pkim),
        get_neigh_(get_neigh),
        baseconvert_(baseconvert),
        cachedNumberContributingParticles_(cachedNumberContributingParticles),
        request_(-baseconvert_),  // set to first value (test-based indexing)
    mode_(1)  // locator mode
  {
    next(i, numnei, n1atom, pRij);
  }
  bool done() const
  {
    return !(request_ + baseconvert_ <= cachedNumberContributingParticles_);
  }
  int next(int* const i, int* const numnei, int** const n1atom,
           double** const pRij)
  {
    int ier;
    // Allow for request_ to be incremented to one more than contributing
    // without causing an error/warning from the openkim-api
    int req = std::min(request_,
                       cachedNumberContributingParticles_-baseconvert_-1);
    ier = (*get_neigh_)(
        reinterpret_cast<void**>(const_cast<KIM_API_model**>(&pkim_)),
        (int*) &mode_,
        &req,
        (int*) i,
        (int*) numnei,
        (int**) n1atom,
        (double**) pRij);
    *i += baseconvert_;  // adjust index of current particle

    ++request_;
    return ier;
  }
};

// Iterator object for Iterator mode access to neighbor list
class IteratorIterator
{
 private:
  KIM_API_model* const pkim_;
  GetNeighborFunction* const get_neigh_;
  int const baseconvert_;
  int request_;
  int const mode_;
  int ier_;
 public:
  IteratorIterator(KIM_API_model* const pkim,
                   GetNeighborFunction* const get_neigh,
                   int const baseconvert,
                   int const cachedNumberContributingParticles,
                   int* const i,
                   int* const numnei,
                   int** const n1atom,
                   double** const pRij)
      : pkim_(pkim),
        get_neigh_(get_neigh),
        baseconvert_(baseconvert),
        request_(0),  // reset iterator
        mode_(0),     // iterator mode
        ier_(KIM_STATUS_FAIL)
  {
    ier_ = (*get_neigh_)(
        reinterpret_cast<void**>(const_cast<KIM_API_model**>(&pkim_)),
        (int*) &mode_, &request_, i, numnei, n1atom, pRij);
    if (ier_ != KIM_STATUS_NEIGH_ITER_INIT_OK)
    {
      pkim->report_error(__LINE__, __FILE__, "iterator init failed", ier_);
      ier_ = KIM_STATUS_FAIL;
    }
    request_ = 1;  // increment iterator
    // return initial set of neighbors
    next(i, numnei, n1atom, pRij);
  }
  bool done() const
  {
    return ((ier_ == KIM_STATUS_FAIL) ||
            (ier_ == KIM_STATUS_NEIGH_ITER_PAST_END));
  }
  int next(int* const i, int* const numnei, int** const n1atom,
           double** const pRij)
  {
    ier_ = (*get_neigh_)(
        reinterpret_cast<void**>(const_cast<KIM_API_model**>(&pkim_)),
        (int*) &mode_,
        &request_,
        (int*) i,
        (int*) numnei,
        (int**) n1atom,
        (double**) pRij);
    *i += baseconvert_;  // adjust index of current particle

    return ier_;
  }
};

// helper routine declarations
void AllocateAndInitialize2DArray(double**& arrayPtr, int const extentZero,
                                  int const extentOne);
void Deallocate2DArray(double**& arrayPtr);
void AllocateAndInitialize3DArray(double***& arrayPtr, int const extentZero,
                                  int const extentOne, int const extentTwo);
void Deallocate3DArray(double***& arrayPtr);

//==============================================================================
//
// Declaration of EAM_Implementation class
//
//==============================================================================

//******************************************************************************
class EAM_Implementation
{
 public:
  EAM_Implementation(KIM_API_model* const pkim,
                           char const* const parameterFileNames,
                           int const parameterFileNameLength,
                           int const numberParameterFiles,
                           int* const ier);
  ~EAM_Implementation();  // no explicit Destroy() needed here

  int Reinit(KIM_API_model* pkim);
  int Compute(KIM_API_model* pkim);

 private:
  // Constant values that never change
  //   Set in constructor (via SetConstantValues)
  //
  //
  // KIM API: Conventions
  int baseconvert_;
  NBCTypeEnum NBCType_;
  bool isHalf_;
  bool isLocatorMode_;
  //
  // KIM API: Model Input indices
  int numberOfSpeciesIndex_;
  int numberOfParticlesIndex_;
  int numberContributingParticlesIndex_;
  int particleSpeciesIndex_;
  int coordinatesIndex_;
  int boxSideLengthsIndex_;
  int get_neighIndex_;
  int process_dEdrIndex_;
  int process_d2Edr2Index_;
  //
  // KIM API: Model Output indices
  int cutoffIndex_;
  int energyIndex_;
  int forcesIndex_;
  int particleEnergyIndex_;
  //
  // EAM_Implementation: constants
  int numberModelSpecies_;
  int numberUniqueSpeciesPairs_;
  EAMFileType eamFileType_;


  // Constant values that are read from the input files and never change
  //   Set in constructor (via functions listed below)
  //
  //
  // KIM API: Model Fixed Parameters
  //   Memory allocated in   AllocateFixedParameterMemory()
  //   Memory deallocated in ~EAM_Implementation()
  //   Data set in ReadParameterFile routines
  char* comments_ptr_[MAX_PARAMETER_FILES];
  char comments_[MAX_PARAMETER_FILES][MAXLINE];
  char particleNames_[MAXLINE];
  int* particleNumber_;
  double* particleMass_;
  double* latticeConstant_;
  char** latticeType_;
  int numberRhoPoints_;
  int numberRPoints_;
  //
  // KIM API: Model Free Parameters whose (pointer) values never change
  //   Memory allocated in   AllocateFreeParameterMemory() (from constructor)
  //   Memory deallocated in ~Eam_Implementation()
  //   Data set in ReadParameterFile routines OR by KIM Simulator
  double** embeddingData_;
  double*** densityData_;
  double*** rPhiData_;

  // Free Parameter pointers to be published without repeat data
  double** publishDensityData_;
  double** publish_rPhiData_;

  // Mutable values that only change when reinit() executes
  //   Set in Reinit (via SetReinitMutableValues)
  //
  //
  // KIM API: Model Free Parameters (can be changed directly by KIM Simulator)
  double cutoffParameter_;
  double deltaR_;
  double deltaRho_;
  //
  // EAM_Implementation: values (changed only by Reinit())
  double cutoffSq_;
  double oneByDr_;
  double oneByDrho_;
  //   Memory allocated once by AllocateFreeParameterMemory()
  double** embeddingCoeff_;
  double*** densityCoeff_;
  double*** rPhiCoeff_;


  // Mutable values that can change with each call to Reinit() and Compute()
  //   Memory may be reallocated on each call
  //
  //
  // EAM_Implementation: values that change
  int cachedNumberOfParticles_;
  int cachedNumberContributingParticles_;
  double* densityValue_;
  double* embeddingDerivativeValue_;
  double* embeddingSecondDerivativeValue_;


  // Helper methods
  //
  //
  // Related to constructor
  int SetConstantValues(KIM_API_model* const pkim);
  int DetermineNBCTypeAndHalf(KIM_API_model* const pkim);
  void AllocateFixedParameterMemory();
  static int OpenParameterFiles(
      KIM_API_model* const pkim,
      char const* const parameterFileNames,
      int const parameterFileNameLength,
      int const numberParameterFiles,
      FILE* parameterFilePointers[MAX_PARAMETER_FILES]);
  static EAMFileType DetermineParameterFileTypes(
      KIM_API_model* const pkim,
      int const numberModelSpecies,
      FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
      int const numberParameterFiles);
  static EAMFileType IsFuncflOrSetfl(FILE* const fptr);
  static EAMFileType IsSetflOrFinnisSinclair(KIM_API_model* const pkim,
                                             FILE* const fptr);
  int ProcessParameterFileHeaders(
      KIM_API_model* const pkim,
      EAMFileType const eamFileType,
      FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
      int const numberParameterFiles, SetOfFuncflData& funcflData);
  int ReadSetflHeader(KIM_API_model* const pkim, FILE* const fptr);
  int ReadFuncflHeader(KIM_API_model* const pkim, FILE* const fptr,
                       int const fileIndex, int& numberRhoPoints,
                       double& deltaRho, int& numberRPoints,
                       double& deltaR, double& cutoffParameter);
  int SetParticleNamesForFuncflModels(KIM_API_model* const pkim);
  void AllocateFreeParameterMemory();
  int ProcessParameterFileData(
      KIM_API_model* const pkim,
      EAMFileType const eamFileType,
      FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
      int const numberParameterFiles, SetOfFuncflData& funcflData);
  int ReadSetflData(KIM_API_model* const pkim, FILE* const fptr);
  int ReadFinnisSinclairData(KIM_API_model* const pkim, FILE* const fptr);
  static int ReadFuncflData(KIM_API_model* const pkim, FILE* const fptr,
                            int const fileIndex,
                            SetOfFuncflData& funcflData);
  static int GrabData(KIM_API_model* const pkim, FILE* const fptr, int const n,
                      double* const list);
  void ReinterpolateAndMix(SetOfFuncflData const& funcflData);
  static void CloseParameterFiles(
      FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
      int const numberParameterFiles);
  int ConvertUnits(KIM_API_model* const pkim);
  int RegisterKIMParameters(KIM_API_model* const pkim,
                            EAMFileType const eamFileType) const;
  int RegisterKIMFunctions(KIM_API_model* const pkim) const;
  //
  // Related to Reinit()
  int SetReinitMutableValues(KIM_API_model* const pkim);
  void SplineInterpolateAllData();
  static void SplineInterpolate(double const* const dat,
                                double const delta, int const n,
                                double* const coe);
  //
  // Related to Compute()
  int SetComputeMutableValues(KIM_API_model* const pkim,
                              bool& isComputeProcess_dEdr,
                              bool& isComputeProcess_d2Edr2,
                              bool& isComputeEnergy,
                              bool& isComputeForces,
                              bool& isComputeParticleEnergy,
                              int const*& particleSpecies,
                              GetNeighborFunction *& get_neigh,
                              double const*& boxSideLengths,
                              VectorOfSizeDIM const*& coordinates,
                              double*& energy,
                              double*& particleEnergy,
                              VectorOfSizeDIM*& forces);
  int CheckParticleSpecies(KIM_API_model* const pkim,
                           int const* const particleSpecies) const;
  int GetComputeIndex(const bool& isComputeProcess_dEdr,
                      const bool& isComputeProcess_d2Edr2,
                      const bool& isComputeEnergy,
                      const bool& isComputeForces,
                      const bool& isComputeParticleEnergy) const;
  static void ApplyMIOPBC(double const* const boxSideLengths, double* const dx);

  // compute functions
  template< class Iter, bool isHalf, RijEnum rijMethod,
            bool isComputeProcess_dEdr, bool isComputeProcess_d2Edr2,
            bool isComputeEnergy, bool isComputeForces,
            bool isComputeParticleEnergy >
  int Compute(KIM_API_model* const pkim,
              const int* const particleSpecies,
              GetNeighborFunction* const get_neigh,
              const double* const boxSideLengths,
              const VectorOfSizeDIM* const coordinates,
              double* const energy,
              VectorOfSizeDIM* const forces,
              double* const particleEnergy);
};

//==============================================================================
//
// Definition of EAM_Implementation::Compute functions
//
// NOTE: Here we rely on the compiler optimizations to prune dead code
//       after the template expansions.  This provides high efficiency
//       and easy maintenance.
//
//==============================================================================

//******************************************************************************
template< class Iter, bool isHalf, RijEnum rijMethod,
          bool isComputeProcess_dEdr, bool isComputeProcess_d2Edr2,
          bool isComputeEnergy, bool isComputeForces,
          bool isComputeParticleEnergy >
int EAM_Implementation::Compute(
    KIM_API_model* const pkim,
    const int* const particleSpecies,
    GetNeighborFunction* const get_neigh,
    const double* const boxSideLengths,
    const VectorOfSizeDIM* const coordinates,
    double* const energy,
    VectorOfSizeDIM* const forces,
    double* const particleEnergy)
{
  int ier = KIM_STATUS_OK;

  // initialize electron density for each contributing particle
  for (int i = 0; i < cachedNumberContributingParticles_; ++i)
  {
    densityValue_[i] = 0.0;
    // no need to initialize embeddingDerivativeValue_
  }
  if (isComputeEnergy == true)
  {
    *energy = 0.0;
  }
  if (isComputeForces == true)
  {
    for (int i = 0; i < cachedNumberOfParticles_; ++i)
    {
      for (int j = 0; j < DIMENSION; ++j)
        forces[i][j] = 0.0;
    }
  }

  // compute electron density
  // Setup loop over contributing particles
  int i = 0;
  int numnei = 0;
  int* n1atom = 0;
  double* pRij = 0;
  for (Iter iterator(pkim, get_neigh, baseconvert_,
                     cachedNumberContributingParticles_, &i, &numnei, &n1atom,
                     &pRij);
       iterator.done() == false;
       ier = iterator.next(&i, &numnei, &n1atom, &pRij))
  {
    if (ier < KIM_STATUS_OK) // check that iterator.next was successful
    {
      pkim->report_error(__LINE__, __FILE__, "get_neigh", ier);
      return ier;
    }

    // Setup loop over neighbors of current particle
    for (int jj = 0; jj < numnei; ++jj)
    {  // adjust index of particle neighbor
      int const j = n1atom[jj] + baseconvert_;
      double* r_ij;
      double r_ijValue[DIMENSION];
      // Compute r_ij appropriately
      switch (rijMethod)
      {
        case Coordinates:
        {
          r_ij = r_ijValue;
          for (int k = 0; k < DIMENSION; ++k)
            r_ij[k] = coordinates[j][k] - coordinates[i][k];
          break;
        }
        case MI_OPBC:
        {
          r_ij = r_ijValue;
          for (int k = 0; k < DIMENSION; ++k)
            r_ij[k] = coordinates[j][k] - coordinates[i][k];
          // apply minimum image convention
          ApplyMIOPBC(boxSideLengths, r_ij);
          break;
        }
        case RVec:
        {
          r_ij = &pRij[jj * DIMENSION];
          break;
        }
      }

      // compute distance squared
      double rij2 = 0.0;
      for (int k = 0; k < DIMENSION; ++k)
        rij2 += r_ij[k] * r_ij[k];

      if (rij2 <= cutoffSq_)
      { // compute contribution to electron density
        double rijOffset;
        int rijIndex;
        double const rij = sqrt(rij2);

        // compute rijOffset and rijIndex
        GET_DELTAX_AND_INDEX(rij, oneByDr_, numberRPoints_, rijOffset,
                             rijIndex);

        // interpolate value of rho_beta(r_ij)
        double densityBetaValue;
        double const* const densityBetaCoeff
            = densityCoeff_[particleSpecies[i]][particleSpecies[j]];
        INTERPOLATE_F(densityBetaCoeff, rijOffset, rijIndex, densityBetaValue);
        densityValue_[i] += densityBetaValue;

        // compute ED contribution from neighbor if half list
        if (isHalf == true)
        {
          if (j < cachedNumberContributingParticles_)
          { // if using half list and j contributes, add its contribution.
            // interpolate value of rho_alpha(r_ij)
            double densityAlphaValue;
            double const* const densityAlphaCoeff
                = densityCoeff_[particleSpecies[j]][particleSpecies[i]];
            INTERPOLATE_F(densityAlphaCoeff, rijOffset, rijIndex,
                          densityAlphaValue);
            densityValue_[j] += densityAlphaValue;
          }
        }
      }
    }  // end of loop over neighbors
  }  // end of loop over contributing particles

  // calculate embedding function and its derivative
  for (Iter iterator(pkim, get_neigh, baseconvert_,
                     cachedNumberContributingParticles_, &i, &numnei, &n1atom,
                     &pRij);
       iterator.done() == false;
       ier = iterator.next(&i, &numnei, &n1atom, &pRij))
  {
    if (ier < KIM_STATUS_OK) // check that iterator.next was successful
    {
      pkim->report_error(__LINE__, __FILE__, "get_neigh", ier);
      return ier;
    }

    double densityOffset;
    int densityIndex;
    // compute densityOffset and densityIndex
    GET_DELTAX_AND_INDEX(densityValue_[i], oneByDrho_, numberRhoPoints_,
                         densityOffset, densityIndex);
    double const* const embeddingAlphaCoeff
        = embeddingCoeff_[particleSpecies[i]];
    if (0 < numnei)
    {
      // interpolate F_i(rho_i)
      double embeddingValue;
      INTERPOLATE_F(embeddingAlphaCoeff, densityOffset, densityIndex,
                    embeddingValue);
      // Contribute embedding term to Energy
      if (isComputeEnergy == true)
      {
        *energy += embeddingValue;
      }
      // Contribute embedding term to ParticleEnergy
      if (isComputeParticleEnergy == true)
      {
        particleEnergy[i] = embeddingValue;
      }
    }
    // Compute embedding derivative
    if ((isComputeForces == true) || (isComputeProcess_dEdr == true) ||
        (isComputeProcess_d2Edr2 == true))
    {
      // interpolate dF_i(rho_i)/d(rho_i)
      INTERPOLATE_DF(embeddingAlphaCoeff, densityOffset, densityIndex,
                     embeddingDerivativeValue_[i]);
    }
    // Compute embedding second derivative
    if (isComputeProcess_d2Edr2 == true)
    {
      // interpolate d^2F_i(rho_i)/d(rho_i)^2
      INTERPOLATE_D2F(embeddingAlphaCoeff, densityOffset, densityIndex,
                      embeddingSecondDerivativeValue_[i]);
    }
  }

  // calculate contribution from electron density to the force part
  // and from pair function
  //
  // Setup loop over contributing particles
  for (Iter iterator(pkim, get_neigh, baseconvert_,
                     cachedNumberContributingParticles_, &i, &numnei, &n1atom,
                     &pRij);
       iterator.done() == false;
       iterator.next(&i, &numnei, &n1atom, &pRij))
  {
    // Setup loop over neighbors of current particle
    for (int jj = 0; jj < numnei; ++jj)
    { // adjust index of particle neighbor
      int const j = n1atom[jj] + baseconvert_;
      double* r_ij;
      double r_ijValue[DIMENSION];
      // Compute r_ij appropriately
      switch (rijMethod)
      {
        case Coordinates:
        {
          r_ij = r_ijValue;
          for (int k = 0; k < DIMENSION; ++k)
            r_ij[k] = coordinates[j][k] - coordinates[i][k];
          break;
        }
        case MI_OPBC:
        {
          r_ij = r_ijValue;
          for (int k = 0; k < DIMENSION; ++k)
            r_ij[k] = coordinates[j][k] - coordinates[i][k];
          // apply minimum image convention
          ApplyMIOPBC(boxSideLengths, r_ij);
          break;
        }
        case RVec:
        {
          r_ij = &pRij[jj * DIMENSION];
          break;
        }
      }

      // compute distance squared
      double rij2 = 0.0;
      for (int k = 0; k < DIMENSION; ++k)
        rij2 += r_ij[k] * r_ij[k];

      if (rij2 <= cutoffSq_)
      { // compute contribution to energy and force
        double rijOffset;
        int rijIndex;
        double const rij = sqrt(rij2);

        // compute rijOffset and rijIndex
        GET_DELTAX_AND_INDEX(rij, oneByDr_, numberRPoints_, rijOffset,
                             rijIndex);

        // interpolate r_ij*phi(r_ij)
        double rijPhiValue;
        double const* const rijPhiAlphaBetaCoeff
            = rPhiCoeff_[particleSpecies[i]][particleSpecies[j]];
        INTERPOLATE_F(rijPhiAlphaBetaCoeff, rijOffset, rijIndex, rijPhiValue);

        // find phi(r_ij)
        double const oneByRij = ONE / rij;
        double const pairPotentialValue = rijPhiValue * oneByRij;

        // Contribute pair term to Energy as half or full
        if (isComputeEnergy == true)
        {
          if (isHalf == true)
          {
            if (j < cachedNumberContributingParticles_)
            {
              *energy += pairPotentialValue;
            }
            else
            {
              *energy += HALF * pairPotentialValue;
            }
          }
          else
          {
            *energy += HALF * pairPotentialValue;
          }
        }

        // Contribute pair term to Particle Energy
        if (isComputeParticleEnergy == true)
        {
          if (isHalf == true)
          {
            if (i < cachedNumberContributingParticles_)
            {
              particleEnergy[i] += HALF * pairPotentialValue;
            }
            if (j < cachedNumberContributingParticles_)
            {
              particleEnergy[j] += HALF * pairPotentialValue;
            }
          }
          else
          {
            particleEnergy[i] += HALF * pairPotentialValue;
          }
        }

        // Compute dEdrByR terms as half or full
        double dEdrByRij = 0.0;
        if ((isComputeForces == true) || (isComputeProcess_dEdr == true))
        {
          if (isHalf == true)
          { // interpolate derivative of r_ij*phi(r_ij) function
            double rijPhiDerivativeValue;
            INTERPOLATE_DF(rijPhiAlphaBetaCoeff, rijOffset, rijIndex,
                           rijPhiDerivativeValue);

            // interpolate derivative of rho_beta(r_ij)
            double densityBetaDerivativeValue;
            double const* const densityBetaCoeff
                = densityCoeff_[particleSpecies[i]][particleSpecies[j]];
            INTERPOLATE_DF(densityBetaCoeff, rijOffset, rijIndex,
                           densityBetaDerivativeValue);

            // compute dEdr contribution
            if (j < cachedNumberContributingParticles_)
            { // interpolate derivative of rho_alpha(r_ij)
              double densityAlphaDerivativeValue;
              double const* const densityAlphaCoeff
                  = densityCoeff_[particleSpecies[j]][particleSpecies[i]];
              INTERPOLATE_DF(densityAlphaCoeff, rijOffset, rijIndex,
                             densityAlphaDerivativeValue);

              // embedding contribution to dEdr
              double const embeddingContribution
                  = ((embeddingDerivativeValue_[i] * densityBetaDerivativeValue)
                     +
                     (embeddingDerivativeValue_[j] *
                      densityAlphaDerivativeValue));

              // pair potential contribution to dEdr
              double const pairPotentialContribution
                  = (rijPhiDerivativeValue - pairPotentialValue) * oneByRij;

              // divide by r so we can multiply by r_ij below
              dEdrByRij = (embeddingContribution + pairPotentialContribution)
                  * oneByRij;
            }
            else
            { // embedding contribution to dEdr
              double const embeddingContribution
                  = (embeddingDerivativeValue_[i] * densityBetaDerivativeValue);

              // pair potential contribution
              double const pairPotentialContribution
                  = HALF * (rijPhiDerivativeValue - pairPotentialValue)
                  * oneByRij;

              // divide by r so we can multiply by r_ij below
              dEdrByRij = (embeddingContribution + pairPotentialContribution)
                  * oneByRij;
            }
          }
          else
          { // interpolate derivative of r_ij*phi(r_ij) function
            double rijPhiDerivativeValue;
            INTERPOLATE_DF(rijPhiAlphaBetaCoeff, rijOffset, rijIndex,
                           rijPhiDerivativeValue);

            // interpolate derivative of rho_beta(r_ij)
            double densityBetaDerivativeValue;
            double const* const densityBetaCoeff
                = densityCoeff_[particleSpecies[i]][particleSpecies[j]];
            INTERPOLATE_DF(densityBetaCoeff, rijOffset, rijIndex,
                           densityBetaDerivativeValue);

            // compute dEdr contribution
            // embedding contribution to dEdr
            double const embeddingContribution
                = (embeddingDerivativeValue_[i] * densityBetaDerivativeValue);

            // pair potential contribution
            double const pairPotentialContribution
                = HALF * (rijPhiDerivativeValue - pairPotentialValue)
                * oneByRij;

            // divide by r so we can multiply by r_ij below
            dEdrByRij = (embeddingContribution + pairPotentialContribution)
                * oneByRij;
          }
        }

        // Contribute dEdrByR to forces
        if (isComputeForces == true)
        {
          for (int k = 0; k < DIMENSION; ++k)
          {
            forces[i][k] += dEdrByRij * r_ij[k];
            forces[j][k] -= dEdrByRij * r_ij[k];
          }
        }

        // Call process_dEdr
        if (isComputeProcess_dEdr == true)
        {
          double const dEdr = dEdrByRij * rij;
          ier = pkim->process_dEdr(const_cast<KIM_API_model**>(&pkim),
                                   const_cast<double*>(&dEdr),
                                   const_cast<double*>(&rij),
                                   &r_ij,
                                   &i,
                                   const_cast<int*>(&j));
          if (ier < KIM_STATUS_OK)
          {
            pkim->report_error(__LINE__, __FILE__, "process_dEdr", ier);
            return ier;
          }
        }
      }  // if particles i and j interact
    }  // end of first neighbor loop
  }  // end of loop over contributing particles

  // Separate loop nest for process_d2Edr2
  if (isComputeProcess_d2Edr2 == true)
  { // Separate implementations for half and full cases
    //
    // For this potential, the second derivative
    // d^2E/dr_{ij}dr_{kl} is nonzero only if i,j,k,l are
    // not distinct. As such, we need only address all
    // derivatives d^2E/dr_{ij}dr_{ik}.
    //
    // If we have full lists, we can simply iterate over
    // particle i's neighbor list in a doubly-nested
    // triangular loop.
    //
    // If we have half lists, then the presence of j in
    // i's neighbor list implies that i is also a neighbor
    // of j, so that we should be processing all derivs
    // d^2E/dr_{ij}dr_{ik} and d^2E/dr_{ji}dr_{jl}. The
    // hurdle to overcome here is the need to access the
    // neighbor list of particle j (and also particle k).
    //
    // We avoid this problem by taking a different
    // approach. For each neighbor (j or k) of i, we save
    // i in a "deferred neighbor list" for each of j and
    // k. For each particle, we process neighbors both in
    // its neighbor list and its deferred neighbor list.
    // (In effect, we are reconstructing full lists for
    // each particle.)

    // Containers for deferred neighbor lists for process_d2Edr2
    std::multimap<int, neighbor> deferredNeighborMap;

    // Iterators over deferred neighbor lists
    std::pair<deferredNeighborIterator, deferredNeighborIterator>
        deferredNeighborRange;
    deferredNeighborIterator deferredNeighborPtrJ;

    // Setup loop over contributing particles
    for (Iter iterator(pkim, get_neigh, baseconvert_,
                       cachedNumberContributingParticles_, &i, &numnei, &n1atom,
                       &pRij);
         iterator.done() == false;
         iterator.next(&i, &numnei, &n1atom, &pRij))
    { // Number of neighbors to visit. May be greater than
      // numnei if previous particles inserted deferred
      // neighbors to this particle, or if this particle
      // is its own neighbor.
      int extendedNumnei;
      extendedNumnei = numnei;

      // Begin by deferring i to each of its neighbors.
      // It is tempting to try to include this in the
      // subsequent loop that also does the calculations,
      // but this cannot be done. The complete deferred
      // list must be available the first time we go
      // through the inner, second-neighbor, loop.
      if (isHalf == true)
      {
        for (int jj = 0; jj < numnei; ++jj)
        {
          neighbor reflectedNeighJ;

          // Index of particle neighbor
          int const atomj = n1atom[jj];
          // Adjust index of particle neighbor
          int const j = atomj + baseconvert_;

          if (j < cachedNumberContributingParticles_)
          {
            // Save "unadjusted" particle index
            reflectedNeighJ.index = i - baseconvert_;

            // Save displacement if necessary
            if (RVec == rijMethod)
            {
              for (int k = 0; k < DIMENSION; ++k)
                reflectedNeighJ.rij[k] = -pRij[jj * DIMENSION + k];
            }

            deferredNeighborMap.insert(std::make_pair(j, reflectedNeighJ));
          }
        }

        // No further particles will contribute to i's neighbor list.
        // Safe to finalize the neighbor count.
        extendedNumnei += deferredNeighborMap.count(i);
        deferredNeighborRange = deferredNeighborMap.equal_range(i);
        deferredNeighborPtrJ  = deferredNeighborRange.first;
      }

      // Setup loop over neighbors of current particle
      for (int jj = 0; jj < extendedNumnei; ++jj)
      {
        int atomj;

        if (isHalf == true && jj >= numnei)
        { // Look in deferred neighbor list
          atomj = deferredNeighborPtrJ->second.index;
        }
        else
        { // Not a deferred neighbor
          atomj = n1atom[jj];
        }

        // adjust index of particle neighbor
        int const j = atomj + baseconvert_;

        // Declare enough space to hold two displacements,
        // one for each neighbor. Ensuring that the two
        // displacements are contiguous in memory allows
        // the use of r_ijValue in the call to process_d2Edr2.
        double r_ijValue[2*DIMENSION];

        // Pointer to the first displacement.
        // Perhaps confusingly, this might be set to point
        // into the RVEC displacements list (pRij), rather
        // than always pointing into r_ijValue.
        double* r_ij;

        // Compute r_ij appropriately
        switch (rijMethod)
        {
          case Coordinates:
          {
            r_ij = r_ijValue;
            for (int k = 0; k < DIMENSION; ++k)
              r_ij[k] = coordinates[j][k] - coordinates[i][k];
            break;
          }
          case MI_OPBC:
          {
            r_ij = r_ijValue;
            for (int k = 0; k < DIMENSION; ++k)
              r_ij[k] = coordinates[j][k] - coordinates[i][k];
            // apply minimum image convention
            ApplyMIOPBC(boxSideLengths, r_ij);
            break;
          }
          case RVec:
          {
            if (isHalf == true && jj >= numnei)
            { // Look in deferred neighbor list
              r_ij = r_ijValue;
              for (int k = 0; k < DIMENSION; ++k)
                r_ij[k] = deferredNeighborPtrJ->second.rij[k];
            }
            else
            {
              r_ij = &pRij[jj * DIMENSION];
            }
            break;
          }
        }

        // compute distance squared
        double rij2 = 0.0;
        for (int k = 0; k < DIMENSION; ++k)
          rij2 += r_ij[k] * r_ij[k];

        if (rij2 <= cutoffSq_)
        {
          double rijOffset;
          int rijIndex;
          double rij[2] = { sqrt(rij2), 0.0 };  // rij[1] used for ik pair

          // compute rijOffset and rijIndex
          GET_DELTAX_AND_INDEX(rij[0], oneByDr_, numberRPoints_, rijOffset,
                               rijIndex);

          // interpolate r_ij*phi(r_ij)
          double rijPhiValue;
          double const* const rijPhiAlphaBetaCoeff
              = rPhiCoeff_[particleSpecies[i]][particleSpecies[j]];
          INTERPOLATE_F(rijPhiAlphaBetaCoeff, rijOffset, rijIndex, rijPhiValue);

          // find phi(r_ij)
          double const oneByRij = ONE / rij[0];
          double const pairPotentialValue = rijPhiValue * oneByRij;

          // interpolate derivative of r_ij*phi(r_ij) function
          double rijPhiDerivativeValue;
          INTERPOLATE_DF(rijPhiAlphaBetaCoeff, rijOffset, rijIndex,
                         rijPhiDerivativeValue);

          // find derivative of phi(r_ij)
          double const pairPotentialDerivativeValue
              = (rijPhiDerivativeValue - pairPotentialValue) * oneByRij;

          // interpolate derivative of rho_beta(r_ij)
          double densityBetaDerivativeValue;
          double const* const densityBetaCoeff
              = densityCoeff_[particleSpecies[i]][particleSpecies[j]];
          INTERPOLATE_DF(densityBetaCoeff, rijOffset, rijIndex,
                         densityBetaDerivativeValue);

          // Setup second loop over neighbors of current particle
          deferredNeighborIterator
              deferredNeighborPtrK = deferredNeighborPtrJ;
          for (int kk = jj; kk < extendedNumnei; ++kk)
          {
            int atomk;
            if (isHalf == true && kk >= numnei)
            { // Look in deferred neighbor list
              atomk = deferredNeighborPtrK->second.index;
            }
            else
            {
              atomk = n1atom[kk];
            }

            // adjust index of particle neighbor
            int const k = atomk + baseconvert_;
            double* r_ik;

            // Pointer to the second displacement vector
            r_ik = &(r_ijValue[DIMENSION]);
            switch (rijMethod)
            {
              case Coordinates:
              {
                for (int d = 0; d < DIMENSION; ++d)
                  r_ik[d] = coordinates[k][d] - coordinates[i][d];
                break;
              }
              case MI_OPBC:
              {
                for (int d = 0; d < DIMENSION; ++d)
                  r_ik[d] = coordinates[k][d] - coordinates[i][d];
                // apply minimum image convention
                ApplyMIOPBC(boxSideLengths, r_ik);
                break;
              }
              case RVec:
              { // In this case, we could have earlier set r_ij
                // to point into pRij, and not into r_ijValue.
                // In other words, r_ijValue might not yet
                // contain the components of the first neighbor's
                // displacement. We need them to be there.
                for (int d = 0; d < DIMENSION; ++d)
                  r_ijValue[d] = r_ij[d];

                r_ij = r_ijValue;

                if (isHalf == true && kk >= numnei)
                { // Look in deferred neighbor list
                  for (int d = 0; d < DIMENSION; ++d)
                    r_ik[d] = deferredNeighborPtrK->second.rij[d];
                }
                else
                {
                  for (int d = 0; d < DIMENSION; ++d)
                    r_ik[d] = pRij[kk * DIMENSION + d];
                }

                break;
              }
            }

            // compute distance squared
            double rik2 = 0.0;
            for (int d = 0; d < DIMENSION; ++d)
              rik2 += r_ik[d] * r_ik[d];

            if (rik2 <= cutoffSq_)
            {
              double rikOffset;
              int rikIndex;
              rij[1] = sqrt(rik2);

              // compute rikOffset and rikIndex
              GET_DELTAX_AND_INDEX(rij[1], oneByDr_, numberRPoints_, rikOffset,
                                   rikIndex);

              // interpolate derivative of rho_gamma(r_ik)
              double densityGammaDerivativeValue;
              double const* const densityGammaCoeff
                  = densityCoeff_[particleSpecies[i]][particleSpecies[k]];
              INTERPOLATE_DF(densityGammaCoeff, rikOffset, rikIndex,
                             densityGammaDerivativeValue);

              // mixed-index embedding contribution to d2Edr2
              double const mixedEmbeddingContribution
                  = embeddingSecondDerivativeValue_[i]
                  * densityBetaDerivativeValue
                  * densityGammaDerivativeValue;

              double d2Edr2 = mixedEmbeddingContribution;

              if (kk == jj)
              { // interpolate second derivative of r_ij*phi(r_ij) function
                double rijPhiSecondDerivativeValue;
                INTERPOLATE_D2F(rijPhiAlphaBetaCoeff, rijOffset, rijIndex,
                                rijPhiSecondDerivativeValue);

                // interpolate second derivative of rho_beta(r_ij)
                double densityBetaSecondDerivativeValue;
                double const* const densityBetaCoeff
                    = densityCoeff_[particleSpecies[i]][particleSpecies[j]];
                INTERPOLATE_D2F(densityBetaCoeff, rijOffset, rijIndex,
                                densityBetaSecondDerivativeValue);

                // pair potential contribution
                double const pairPotentialSecondDerivativeValue
                    = (rijPhiSecondDerivativeValue
                       - TWO * pairPotentialDerivativeValue) * oneByRij;
                double const pairPotentialContribution
                    = HALF * pairPotentialSecondDerivativeValue;

                // second derivative of embedding contribution
                double const embeddingContribution
                    = embeddingDerivativeValue_[i]
                    * densityBetaSecondDerivativeValue;

                d2Edr2 += pairPotentialContribution + embeddingContribution;
              }
              else
              { // Process transpose pair
                const int iis[2] = { i, i };
                const int jjs[2] = { k, j };
                int const* const piis = &iis[0];
                int const* const pjjs = &jjs[0];
                const double rikrij[2] = { rij[1], rij[0] };
                double const* const prikrij = &rikrij[0];
                const double r_ikr_ij[6] = { r_ij[3], r_ij[4], r_ij[5],
                                             r_ij[0], r_ij[1], r_ij[2] };
                double const* const pr_ikr_ij = &r_ikr_ij[0];
                ier = pkim->process_d2Edr2(
                    const_cast<KIM_API_model**>(&pkim),
                    &d2Edr2,
                    const_cast<double**>(&prikrij),
                    const_cast<double**>(&pr_ikr_ij),
                    const_cast<int**>(&piis),
                    const_cast<int**>(&pjjs));
                if (ier < KIM_STATUS_OK)
                {
                  pkim->report_error(__LINE__, __FILE__, "process_d2Edr2", ier);
                  return ier;
                }
              }

              int const iis[2] = { i, i };
              int const jjs[2] = { j, k };
              int const* const piis = &iis[0];
              int const* const pjjs = &jjs[0];
              double const* const prij = &rij[0];
              ier = pkim->process_d2Edr2(
                  const_cast<KIM_API_model**>(&pkim),
                  &d2Edr2,
                  const_cast<double**>(&prij),
                  &r_ij,
                  const_cast<int**>(&piis),
                  const_cast<int**>(&pjjs));
              if (ier < KIM_STATUS_OK)
              {
                pkim->report_error(__LINE__, __FILE__, "process_d2Edr2", ier);
                return ier;
              }
            }  // if particles i and k interact

            if (isHalf == true && kk >= numnei)
            {
              ++deferredNeighborPtrK;
            }
          }  // end of second neighbor loop
        }  // if particles i and j interact

        if (isHalf == true && jj >= numnei)
        {
          ++deferredNeighborPtrJ;
        }
      }  // end of first neighbor loop

      deferredNeighborMap.erase(i);
    }  // end of loop over contributing particles
  }  // if (isComputeProcess_d2Edr2 == true)

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

#endif  // EAM_IMPLEMENTATION_HPP_
