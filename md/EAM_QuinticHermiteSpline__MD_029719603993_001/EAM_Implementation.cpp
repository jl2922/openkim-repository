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
// Copyright (c) 2013--2014, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Ryan S. Elliott
//    Stephen M. Whalen
//


#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

#include "EAM.hpp"
#include "EAM_Implementation.hpp"
#include "KIM_API_status.h"

#define IGNORE_RESULT(fn) if(fn){}


//==============================================================================
//
// Implementation of EAM_Implementation public member functions
//
//==============================================================================

//******************************************************************************
EAM_Implementation::EAM_Implementation(
    KIM_API_model* const pkim,
    char const* const  parameterFileNames,
    int const parameterFileNameLength,
    int const numberParameterFiles,
    int* const ier)
    : numberOfSpeciesIndex_(-1),  // initizlize index, pointer, and cached
      numberOfParticlesIndex_(-1),    // member variables
      numberContributingParticlesIndex_(-1),
      particleSpeciesIndex_(-1),
      coordinatesIndex_(-1),
      boxSideLengthsIndex_(-1),
      get_neighIndex_(-1),
      process_dEdrIndex_(-1),
      process_d2Edr2Index_(-1),
      cutoffIndex_(-1),
      energyIndex_(-1),
      forcesIndex_(-1),
      particleEnergyIndex_(-1),
      particleNumber_(0),
      particleMass_(0),
      latticeConstant_(0),
      latticeType_(0),
      embeddingData_(0),
      densityData_(0),
      rPhiData_(0),
      publishDensityData_(0),
      publish_rPhiData_(0),
      embeddingCoeff_(0),
      densityCoeff_(0),
      rPhiCoeff_(0),
      cachedNumberOfParticles_(0),
      cachedNumberContributingParticles_(0),
      densityValue_(0),
      embeddingDerivativeValue_(0),
      embeddingSecondDerivativeValue_(0)
{ // initialize comments to null strings and set pointers for comment fields
  for (int i = 0; i < MAX_PARAMETER_FILES; ++i)
  {
    comments_[i][0] = 0;
    comments_ptr_[i] = comments_[i];
  }

  // set particleNames to null string
  particleNames_[0] = 0;

  *ier = SetConstantValues(pkim);
  if (*ier < KIM_STATUS_OK) return;

  AllocateFixedParameterMemory();

  FILE* parameterFilePointers[MAX_PARAMETER_FILES];
  *ier = OpenParameterFiles(pkim, parameterFileNames, parameterFileNameLength,
                            numberParameterFiles, parameterFilePointers);
  if (*ier < KIM_STATUS_OK) return;

  eamFileType_ = DetermineParameterFileTypes(pkim,
                                             numberModelSpecies_,
                                             parameterFilePointers,
                                             numberParameterFiles);
  if (eamFileType_ == Error)
  {
    *ier = KIM_STATUS_FAIL;
    return;
  }

  SetOfFuncflData funcflData;
  *ier = ProcessParameterFileHeaders(pkim, eamFileType_, parameterFilePointers,
                                     numberParameterFiles, funcflData);
  if (*ier < KIM_STATUS_OK)
  {
    CloseParameterFiles(parameterFilePointers, numberParameterFiles);
    return;
  }

  AllocateFreeParameterMemory();

  *ier = ProcessParameterFileData(pkim, eamFileType_, parameterFilePointers,
                                  numberParameterFiles, funcflData);
  CloseParameterFiles(parameterFilePointers, numberParameterFiles);
  if (*ier < KIM_STATUS_OK) return;

  *ier = ConvertUnits(pkim);
  if (*ier < KIM_STATUS_OK) return;

  *ier = SetReinitMutableValues(pkim);
  if (*ier < KIM_STATUS_OK) return;

  *ier = RegisterKIMParameters(pkim, eamFileType_);
  if (*ier < KIM_STATUS_OK) return;

  *ier = RegisterKIMFunctions(pkim);
  if (*ier < KIM_STATUS_OK) return;

  // everything is good
  *ier = KIM_STATUS_OK;
  return;
}

//******************************************************************************
EAM_Implementation::~EAM_Implementation()
{ // note: it is ok to delete a null pointer and we have ensured that
  // everything is initialized to null
  delete [] particleNumber_;
  delete [] particleMass_;
  delete [] latticeConstant_;
  if (latticeType_ != 0) delete [] latticeType_[0];
  delete [] latticeType_;

  Deallocate2DArray(embeddingData_);
  Deallocate3DArray(densityData_);
  Deallocate3DArray(rPhiData_);
  Deallocate2DArray(embeddingCoeff_);
  Deallocate3DArray(densityCoeff_);
  Deallocate3DArray(rPhiCoeff_);
  Deallocate2DArray(publishDensityData_);
  Deallocate2DArray(publish_rPhiData_);

  delete [] densityValue_;
  delete [] embeddingDerivativeValue_;
  delete [] embeddingSecondDerivativeValue_;
}

//******************************************************************************
int EAM_Implementation::Reinit(KIM_API_model* const pkim)
{
  int ier;

  for (int i = 0; i < numberModelSpecies_ ; i++)
  {
    for(int j = i; j < numberModelSpecies_; j++)
    {
      int indxPhi = i*numberModelSpecies_ + j - (i*i + i)/2;
      for(int k = 0; k < numberRPoints_; k++)
      {
        rPhiData_[i][j][k] = rPhiData_[j][i][k] = publish_rPhiData_[indxPhi][k];
      }
    }

    for(int j = 0; j < numberModelSpecies_; j++)
    {
      int indxDensity = (eamFileType_ == FinnisSinclair) ?
          i*numberModelSpecies_ + j : i;
      for(int k = 0; k < numberRPoints_; k++)
      {
        densityData_[i][j][k] = publishDensityData_[indxDensity][k];
      }
    }
  }

  ier = SetReinitMutableValues(pkim);
  if (ier < KIM_STATUS_OK) return ier;

  // nothing else to do for this case

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int EAM_Implementation::Compute(KIM_API_model* const pkim)
{
  int ier;

  // KIM API Model Input compute flags
  bool isComputeProcess_dEdr;
  bool isComputeProcess_d2Edr2;
  //
  // KIM API Model Output compute flags
  bool isComputeEnergy;
  bool isComputeForces;
  bool isComputeParticleEnergy;
  //
  // KIM API Model Input
  int const* particleSpecies = 0;
  GetNeighborFunction * get_neigh = 0;
  double const* boxSideLengths = 0;
  VectorOfSizeDIM const* coordinates = 0;
  //
  // KIM API Model Output
  double* energy = 0;
  double* particleEnergy = 0;
  VectorOfSizeDIM* forces = 0;
  ier = SetComputeMutableValues(pkim, isComputeProcess_dEdr,
                                isComputeProcess_d2Edr2, isComputeEnergy,
                                isComputeForces, isComputeParticleEnergy,
                                particleSpecies, get_neigh, boxSideLengths,
                                coordinates, energy, particleEnergy, forces);
  if (ier < KIM_STATUS_OK) return ier;

  ier = CheckParticleSpecies(pkim, particleSpecies);
  if (ier < KIM_STATUS_OK) return ier;

#include "EAM_ImplementationComputeDispatch.cpp"
  return ier;
}

//==============================================================================
//
// Implementation of EAM_Implementation private member functions
//
//==============================================================================

//******************************************************************************
int EAM_Implementation::SetConstantValues(KIM_API_model* const pkim)
{
  int ier = KIM_STATUS_FAIL;

  // get baseconvert value from KIM API object
  baseconvert_ = pkim->get_model_index_shift();

  ier = DetermineNBCTypeAndHalf(pkim);
  if (ier < KIM_STATUS_OK) return ier;

  // set isLocator
  if (NBCType_ != Cluster)
  { // all NBC cases except CLUSTER
    isLocatorMode_ = pkim->get_neigh_mode(&ier) == 2;  // Locator mode
    if (ier < KIM_STATUS_OK)
    {
      pkim->report_error(__LINE__, __FILE__, "get_neigh_mode", ier);
      return ier;
    }
  }
  else
  {  // use true as a default value for CLUSTER
    isLocatorMode_ = true;
  }

  // obtain indices for various KIM API Object arguments
  pkim->getm_index(
      &ier, 3 * 13,
      "numberOfSpecies",             &numberOfSpeciesIndex_,             1,
      "numberOfParticles",           &numberOfParticlesIndex_,           1,
      "numberContributingParticles", &numberContributingParticlesIndex_, 1,
      "particleSpecies",             &particleSpeciesIndex_,             1,
      "coordinates",                 &coordinatesIndex_,                 1,
      "boxSideLengths",              &boxSideLengthsIndex_,              1,
      "get_neigh",                   &get_neighIndex_,                   1,
      "process_dEdr",                &process_dEdrIndex_,                1,
      "process_d2Edr2",              &process_d2Edr2Index_,              1,
      "cutoff",                      &cutoffIndex_,                      1,
      "energy",                      &energyIndex_,                      1,
      "forces",                      &forcesIndex_,                      1,
      "particleEnergy",              &particleEnergyIndex_,              1);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "getm_index", ier);
    return ier;
  }

  // set numberModelSpecies
  int dummy;
  ier = pkim->get_num_model_species(&numberModelSpecies_, &dummy);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "get_num_model_species", ier);
    return ier;
  }
  numberUniqueSpeciesPairs_ = ((numberModelSpecies_+1)*numberModelSpecies_)/2;

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int EAM_Implementation::DetermineNBCTypeAndHalf(KIM_API_model* const pkim)
{
  int ier;
  const char* NBC;

  // record the active NBC method
  ier = pkim->get_NBC_method(&NBC);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "get_NBC_method", ier);
    return ier;
  }
  if (strcmp(NBC, "NEIGH_RVEC_H") == 0)
  {
    NBCType_ = Neigh_Rvec;
    isHalf_ = true;
  }
  else if (strcmp(NBC, "NEIGH_PURE_H") == 0)
  {
    NBCType_ = Neigh_Pure;
    isHalf_ = true;
  }
  else if (strcmp(NBC, "NEIGH_RVEC_F") == 0)
  {
    NBCType_ = Neigh_Rvec;
    isHalf_ = false;
  }
  else if (strcmp(NBC, "NEIGH_PURE_F") == 0)
  {
    NBCType_ = Neigh_Pure;
    isHalf_ = false;
  }
  else if (strcmp(NBC, "MI_OPBC_H") == 0)
  {
    NBCType_ = Mi_Opbc;
    isHalf_ = true;
  }
  else if (strcmp(NBC, "MI_OPBC_F") == 0)
  {
    NBCType_ = Mi_Opbc;
    isHalf_ = false;
  }
  else if (strcmp(NBC, "CLUSTER") == 0)
  {
    NBCType_ = Cluster;
    isHalf_ = true;
  }
  else
  {
    ier = KIM_STATUS_FAIL;
    pkim->report_error(__LINE__, __FILE__,
                       "(eam_init_) unknown NBC type", ier);
    return ier;
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
void EAM_Implementation::AllocateFixedParameterMemory()
{ // allocate memory for particle number, mass, lattice constant, and lattice
  // type
  particleNumber_ = new int[numberModelSpecies_];
  particleMass_ = new double[numberModelSpecies_];
  latticeConstant_ = new double[numberModelSpecies_];
  latticeType_ = new char*[numberModelSpecies_];
  latticeType_[0] = new char[numberModelSpecies_ * MAXLINE];
  for (int i = 1; i < numberModelSpecies_; ++i)
  {
    latticeType_[i] = latticeType_[i-1] + MAXLINE;
  }
}

//******************************************************************************
int EAM_Implementation::OpenParameterFiles(
    KIM_API_model* const pkim,
    char const* const parameterFileNames,
    int const parameterFileNameLength,
    int const numberParameterFiles,
    FILE* parameterFilePointers[MAX_PARAMETER_FILES])
{
  int ier;

  for (int i = 0; i < numberParameterFiles; ++i)
  {
    parameterFilePointers[i]
        = fopen(&parameterFileNames[i * parameterFileNameLength], "r");
    if (parameterFilePointers[i] == 0)
    {
      char message[MAXLINE];
      sprintf(message,
              "EAM parameter file number %d cannot be opened",
              i);
      ier = KIM_STATUS_FAIL;
      pkim->report_error(__LINE__, __FILE__, message, ier);
      for (int j = i - 1; i <= 0; --i)
      {
        fclose(parameterFilePointers[j]);
      }
      return ier;
    }
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
EAMFileType EAM_Implementation::DetermineParameterFileTypes(
    KIM_API_model* const pkim,
    int const numberModelSpecies,
    FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
    int const numberParameterFiles)
{
  if ((numberParameterFiles > 1) &&
      (numberParameterFiles <= MAX_PARAMETER_FILES))
  { // should be a set of Funcfl files
    for (int i = 0; i < numberParameterFiles; ++i)
    {
      if (IsFuncflOrSetfl(parameterFilePointers[i]) != Funcfl)
      {
        char message[MAXLINE];
        sprintf(message, "EAM parameter file number %d is not a"
                " funcfl file", i);
        pkim->report_error(__LINE__, __FILE__, message, KIM_STATUS_FAIL);
        return Error;
      }
    }
    // all files are Funcfl as needed, but also check that we have the
    // right nubmer of files
    if (numberModelSpecies != numberParameterFiles)
    {
      pkim->report_error(__LINE__, __FILE__, "Wrong number of parameter"
                         " files in EAM", KIM_STATUS_FAIL);
      return Error;
    }

    // everything is good
    return Funcfl;
  }
  else if (numberParameterFiles == 1)
  {
    EAMFileType eamFileType = IsFuncflOrSetfl(parameterFilePointers[0]);

    // check that we have the right number of files
    if ((eamFileType == Funcfl) && (numberModelSpecies != 1))
    {
      pkim->report_error(__LINE__, __FILE__, "Wrong number of parameter"
                         " files in EAM", KIM_STATUS_FAIL);
      return Error;
    }

    if (eamFileType == Error)
    {
      pkim->report_error(__LINE__, __FILE__, "Unable to determine parameter"
                         " file type in EAM", KIM_STATUS_FAIL);
    }

    // distinguish between setfl and Finnis-Sinclair files
    if (eamFileType == Setfl)
    {
      eamFileType = IsSetflOrFinnisSinclair(pkim, parameterFilePointers[0]);
    }

    return eamFileType;
  }
  else
  {
    char message[MAXLINE];
    sprintf(message, "Invalid number (%d) of parameter files in EAM_",
            numberParameterFiles);
    pkim->report_error(__LINE__, __FILE__, message, KIM_STATUS_FAIL);
    return Error;
  }
}

//******************************************************************************
EAMFileType EAM_Implementation::IsFuncflOrSetfl(FILE* const fptr)
{
  int const numberOfLinesToRead = 8;
  // use 1-based counting for line numbers
  bool isInteger[numberOfLinesToRead + 1];
  bool isDouble[numberOfLinesToRead + 1];
  int intValue[numberOfLinesToRead + 1];

  // discard 1st line.  It is a comment in both funcfl and setfl
  char line[MAXLINE];
  char const* cer = fgets(line, MAXLINE, fptr);
  if (cer == 0) return Error;

  // loop over the lines 2--6
  for (int i = 2; i <= numberOfLinesToRead; ++i)
  { // get next line
    cer = fgets(line, MAXLINE, fptr);
    if (cer == 0) return Error;

    // get first token on line
    char const* const tok = strtok(line, " ,\t\n");
    if (tok == 0)
    { // nothing on the line
      isInteger[i] = false;
      isDouble[i] = false;

      // skip to next line
      continue;
    }

    char* endptr;
    intValue[i] = strtol(tok, &endptr, 10);
    if (*endptr == 0)  // entire string used up, thus, tok was an int
    {
      isInteger[i] = true;
      isDouble[i] = false;
    }
    else
    {
      IGNORE_RESULT(strtod(tok, &endptr));
      if (*endptr == 0)  // entire string used up, thus, tok was a double
      {
        isInteger[i] = false;
        isDouble[i] = true;
      }
      else
      {
        isInteger[i] = false;
        isDouble[i] = false;
      }
    }
  }

  // done with the file rewind it for later use
  rewind(fptr);

  bool const isFuncfl = (  // line 2 starts with "ielem" and <= 118
      (isInteger[2] && (intValue[2] <= 118)) &&
      // line 3 starts with "nrho"
      isInteger[3] &&
      // line 4 starts with embedding data
      isDouble[4] &&
      // line 5 is double data
      isDouble[5] &&
      // line 6 is double data
      isDouble[6] &&
      // line 7 is double data
      isDouble[7] &&
      // line 8 is double data
      isDouble[8]);

  if (isFuncfl) return Funcfl;

  bool const isSetfl = (  // line 4 starts with "ntypes"
      isInteger[4] &&
      // line 5 starts with "nrho"
      isInteger[5] &&
      // line 6 starts with "ielem(1)" and <= 118
      (isInteger[6] && (intValue[6] <= 118)) &&
      // line 7 starts with embedding data
      isDouble[7] &&
      // line 8 is double data
      isDouble[8]);

  if (isSetfl) return Setfl;

  return Error;
}

EAMFileType EAM_Implementation::IsSetflOrFinnisSinclair(
    KIM_API_model* const pkim, FILE* const fptr)
{ // We are free to assume the file format (of the header, at least)
  // is conforms to the Setfl format

  char line[MAXLINE];
  // discard first three lines.  They are comments
  for (int i = 0; i < 3; ++i)
  {
    char const* cer = fgets(line, MAXLINE, fptr);
    if (cer == 0) return Error;
  }

  // 4th line; read number of elements
  int Nelements;
  {
    char const* const cer = fgets(line, MAXLINE, fptr);
    if (cer == 0) return Error;
    int ier = sscanf(line, "%d", &Nelements);
    if (ier != 1) return Error;
  }

  // 5th line; read Nrho and Nr
  int Nrho;
  int Nr;
  {
    char const* const cer = fgets(line, MAXLINE, fptr);
    if (cer == 0) return Error;
    double dummy;
    int const ier = sscanf(line, "%d %lg %d", &Nrho, &dummy, &Nr);
    if (ier != 3) return Error;
  }

  // 6th line; discard
  {
    char const* cer = fgets(line, MAXLINE, fptr);
    if (cer == 0) return Error;
  }

  // Read element1 data set
  {
    double* const dummy = new double[(Nrho>Nr)? Nrho : Nr];
    int ier = GrabData(pkim, fptr, Nrho, dummy);
    if (ier < KIM_STATUS_OK)
    {
      delete[] dummy;
      return Error;
    }

    ier = GrabData(pkim, fptr, Nr, dummy);
    if (ier < KIM_STATUS_OK)
    {
      delete[] dummy;
      return Error;
    }

    delete[] dummy;
  }

  // Read next line (Setfl - element header; FinnisSinclair - rho(r) values)
  bool isSetfl;
  {
    char const* cer = fgets(line, MAXLINE, fptr);
    if (cer == 0) return Error;

    // get first token on line
    char const* const tok = strtok(line, " ,\t\n");
    if (tok == 0)
    { // nothing on the line
      return Error;
    }

    char* endptr;
    IGNORE_RESULT(strtol(tok, &endptr, 10));
    if (*endptr == 0)  // entire string used up, thus, tok was an int
    {
      isSetfl = true;
    }
    else
    {
      IGNORE_RESULT(strtod(tok, &endptr));
      if (*endptr == 0)  // entire string used up, thus, tok was a double
      {
        isSetfl = false;
      }
      else
      {
        return Error;
      }
    }
  }

  // done with the file rewind it for later use
  rewind(fptr);

  if (isSetfl)
    return Setfl;
  else
    return FinnisSinclair;
}

//******************************************************************************
int EAM_Implementation::ProcessParameterFileHeaders(
    KIM_API_model* const pkim,
    EAMFileType const eamFileType,
    FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
    int const numberParameterFiles,
    SetOfFuncflData& funcflData)
{
  int ier;

  switch (eamFileType)
  {
    case FinnisSinclair:
    { // FinnisSinclair has same header information as setfl
      ier = ReadSetflHeader(pkim, parameterFilePointers[0]);
      if (ier < KIM_STATUS_OK)
      {
        pkim->report_error(__LINE__, __FILE__, "ReadSetflHeader", ier);
        return ier;
      }
      break;
    }
    case Setfl:
    {
      ier = ReadSetflHeader(pkim, parameterFilePointers[0]);
      if (ier < KIM_STATUS_OK)
      {
        pkim->report_error(__LINE__, __FILE__, "ReadSetflHeader", ier);
        return ier;
      }
      break;
    }
    case Funcfl:
    { // initialize grid values
      deltaRho_ = 0.0;
      deltaR_ = 0.0;
      cutoffParameter_ = 0.0;
      double rhoMax = 0.0;
      double rMax = 0.0;

      for (int i = 0; i < numberParameterFiles; ++i)
      {
        ier = ReadFuncflHeader(pkim,
                               parameterFilePointers[i],
                               i,
                               funcflData.numberRhoPoints[i],
                               funcflData.deltaRho[i],
                               funcflData.numberRPoints[i],
                               funcflData.deltaR[i],
                               funcflData.cutoff[i]);
        if (ier < KIM_STATUS_OK)
        {
          pkim->report_error(__LINE__, __FILE__, "ReadFuncflHeader", ier);
          return ier;
        }

        // update grid values (use maximums)
        deltaRho_ = std::max(deltaRho_, funcflData.deltaRho[i]);
        deltaR_ = std::max(deltaR_, funcflData.deltaR[i]);
        cutoffParameter_ = std::max(cutoffParameter_, funcflData.cutoff[i]);
        rhoMax = std::max(rhoMax, ((funcflData.numberRhoPoints[i] - 1) *
                                   funcflData.deltaRho[i]));
        rMax = std::max(rMax, ((funcflData.numberRPoints[i] - 1) *
                               funcflData.deltaR[i]));
      }

      // determine number of rho and r points in grid.
      // add 1 to account for point at zero
      numberRhoPoints_ = int(rhoMax / deltaRho_ + 0.5) + 1;
      numberRPoints_ = int(rMax / deltaR_ + 0.5) + 1;

      // set partileNames_
      ier = SetParticleNamesForFuncflModels(pkim);
      if (ier < KIM_STATUS_OK) return ier;
      break;
    }
    default:  // should never get here
    {
      ier = KIM_STATUS_FAIL;
      pkim->report_error(__LINE__, __FILE__, "Bad parameter file(s)", ier);
      return ier;
      break;
    }
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int EAM_Implementation::ReadSetflHeader(KIM_API_model* const pkim,
                                        FILE* const fptr)
{
  int ier;
  char const* cer;
  char line[MAXLINE];

  // read lines 1, 2, and 3 (comment lines)
  for (int i = 0; i < NUMBER_SETFL_COMMENT_LINES; ++i)
  {
    cer = fgets(&comments_[i][0], MAXLINE, fptr);
    if (cer == 0)
    {
      ier = KIM_STATUS_FAIL;
      pkim->report_error(__LINE__, __FILE__,
                         "error reading comment lines of setfl file",
                         ier);
      return ier;
    }
    int const cmntlength = strlen(&comments_[i][0]);
    if (comments_[i][cmntlength-1] == '\n') comments_[i][cmntlength-1] = 0;
  }

  // read 4th line (Nelements Element1 Element2 ... ElementN)
  cer = fgets(particleNames_, MAXLINE, fptr);
  int const nameslength = strlen(particleNames_);
  if (particleNames_[nameslength-1] == '\n') particleNames_[nameslength-1] = 0;
  // parse number of particle species
  int number_of_species;
  ier = sscanf(particleNames_, "%d", &number_of_species);
  if ((cer == 0) || (ier != 1))
  {
    ier = KIM_STATUS_FAIL;
    pkim->report_error(__LINE__, __FILE__,
                       "error reading 4th line of setfl file", ier);
    return ier;
  }
  // check consistency of particle species
  if (numberModelSpecies_ != number_of_species)
  {
    ier = KIM_STATUS_FAIL;
    pkim->report_error(__LINE__, __FILE__, "setfl file and .kim file"
                       " are inconsistent: number of particles don't match",
                       ier);
    return ier;
  }
  // parse the remainder of 4th line for particle specie names
  char* const tmpnames = new char[strlen(particleNames_)+1];
  strcpy(tmpnames, particleNames_);
  char** const elems = new char*[numberModelSpecies_];
  char* tmpstring = strtok(tmpnames, " ,\t");  // ignore first token
  int counter = 0;
  while (tmpstring != 0)
  {
    tmpstring = strtok(0, " ,\t\n");
    elems[counter] = tmpstring;
    ++counter;
    if (counter >= numberModelSpecies_)
    {
      break;
    }
  }
  // check to ensure that the atom species match with the KIM descriptor file
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    const int acode = pkim->get_species_code(elems[i], &ier);
    if (ier < KIM_STATUS_OK)
    { // the atom specie is not listed in kim file or other error
      pkim->report_error(__LINE__, __FILE__, "setfl file and .kim"
                         " file are inconsistent: elements (particles) don't"
                         " match", ier);
      delete [] elems;
      delete [] tmpnames;
      return ier;
    }
    if (acode != i)
    { // order of atom species not consistent with kim file list
      ier = KIM_STATUS_FAIL;
      pkim->report_error(__LINE__, __FILE__,
                         "setfl file and .kim file are inconsistent:"
                         "elements (particles) are out of order", ier);
      delete [] elems;
      delete [] tmpnames;
      return ier;
    }
  }
  delete [] elems;
  delete [] tmpnames;

  // read 5th line (Nrho, deltaRho, Nr, deltaR, cutoff)
  cer = fgets(line, MAXLINE, fptr);
  ier = sscanf(line, "%d %lg %d %lg %lg", &numberRhoPoints_, &deltaRho_,
               &numberRPoints_, &deltaR_, &cutoffParameter_);
  if ((cer == 0) || (ier != 5))
  {
    ier = KIM_STATUS_FAIL;
    pkim->report_error(__LINE__, __FILE__,
                       "error reading 5th lines of setfl file", ier);
    return ier;
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int EAM_Implementation::ReadFuncflHeader(KIM_API_model* const pkim,
                                         FILE* const fptr,
                                         int const fileIndex,
                                         int& numberRhoPoints,
                                         double& deltaRho,
                                         int& numberRPoints,
                                         double& deltaR,
                                         double& cutoffParameter)
{
  int ier;
  char const* cer;
  char line[MAXLINE];

  // read 1st line (comment line)
  cer = fgets(&comments_[fileIndex][0], MAXLINE, fptr);
  if (cer == 0)
  {
    ier = KIM_STATUS_FAIL;
    pkim->report_error(__LINE__, __FILE__,
                       "error reading 1st line of funcfl file", ier);
    return ier;
  }
  int const cmntlength = strlen(&comments_[fileIndex][0]);
  if (comments_[fileIndex][cmntlength-1] == '\n')
    comments_[fileIndex][cmntlength-1] = 0;

  // read 2nd line (particle number, mass, lattice constant, and lattice type)
  cer = fgets(line, MAXLINE, fptr);
  ier = sscanf(line, "%d %lg %lg %s", &(particleNumber_[fileIndex]),
               &(particleMass_[fileIndex]), &(latticeConstant_[fileIndex]),
               latticeType_[fileIndex]);
  if ((cer == 0) || (ier != 4))
  {
    ier = KIM_STATUS_FAIL;
    pkim->report_error(__LINE__, __FILE__,
                       "error 2nd line of funcfl file", ier);
    return ier;
  }

  // read 3rd line (Nrho, deltaRho, Nr, deltaR, cutoff)
  cer = fgets(line, MAXLINE, fptr);
  ier = sscanf(line, "%d %lg %d %lg %lg", &numberRhoPoints, &deltaRho,
               &numberRPoints, &deltaR, &cutoffParameter);
  if ((cer == 0) || (ier != 5))
  {
    ier = KIM_STATUS_FAIL;
    pkim->report_error(__LINE__, __FILE__,
                       "error reading 3rd line of funcfl file", ier);
    return ier;
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int EAM_Implementation::SetParticleNamesForFuncflModels(
    KIM_API_model* const pkim)
{
  int ier;

  // get and correctly order the particle names
  const char** const particleNames = new const char*[numberModelSpecies_];
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    const char* kimModelParticleSpecies;
    ier = pkim->get_model_species(i, &kimModelParticleSpecies);
    if (ier < KIM_STATUS_OK)
    {
      pkim->report_error(__LINE__, __FILE__, "get_model_species", ier);
      delete [] particleNames;
      return ier;
    }
    int const index = pkim->get_species_code(kimModelParticleSpecies, &ier);
    particleNames[index] = kimModelParticleSpecies;
  }

  // write particleNames_ string
  sprintf(particleNames_, "%d ", numberModelSpecies_);
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    strcat(particleNames_, particleNames[i]);
    strcat(particleNames_, " ");
  }
  int const nmlength = strlen(particleNames_);
  particleNames_[nmlength - 1] = 0;
  delete [] particleNames;

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
void EAM_Implementation::AllocateFreeParameterMemory()
{
  // allocate memory for data
  AllocateAndInitialize2DArray(embeddingData_, numberModelSpecies_,
                               numberRhoPoints_);
  AllocateAndInitialize3DArray(densityData_, numberModelSpecies_,
                               numberModelSpecies_, numberRPoints_);
  AllocateAndInitialize3DArray(rPhiData_, numberModelSpecies_,
                               numberModelSpecies_, numberRPoints_);
  // allocate memory for non-repeat data
  AllocateAndInitialize2DArray(
      publishDensityData_,
      numberModelSpecies_
      *((eamFileType_ == FinnisSinclair) ? numberModelSpecies_ : 1),
      numberRPoints_);
  AllocateAndInitialize2DArray(publish_rPhiData_, numberUniqueSpeciesPairs_,
                               numberRPoints_);

  // allocate memory for coefficients
  AllocateAndInitialize2DArray(embeddingCoeff_, numberModelSpecies_,
                               numberRhoPoints_ * NUMBER_SPLINE_COEFF);
  AllocateAndInitialize3DArray(densityCoeff_, numberModelSpecies_,
                               numberModelSpecies_,
                               numberRPoints_ * NUMBER_SPLINE_COEFF);
  AllocateAndInitialize3DArray(rPhiCoeff_, numberModelSpecies_,
                               numberModelSpecies_,
                               numberRPoints_ * NUMBER_SPLINE_COEFF);
}

//******************************************************************************
int EAM_Implementation::ProcessParameterFileData(
    KIM_API_model* const pkim,
    EAMFileType const eamFileType,
    FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
    int const numberParameterFiles,
    SetOfFuncflData& funcflData)
{
  int ier;

  // read data file based on indicated type
  switch (eamFileType)
  {
    case FinnisSinclair:
    {
      ier = ReadFinnisSinclairData(pkim, parameterFilePointers[0]);
      if (ier < KIM_STATUS_OK)
      {
        pkim->report_error(__LINE__, __FILE__, "ReadSetflData", ier);
        return ier;
      }
      break;
    }
    case Setfl:
    {
      ier = ReadSetflData(pkim, parameterFilePointers[0]);
      if (ier < KIM_STATUS_OK)
      {
        pkim->report_error(__LINE__, __FILE__, "ReadSetflData", ier);
        return ier;
      }
      break;
    }
    case Funcfl:
    {
      for (int i = 0; i < numberParameterFiles; ++i)
      { // allocate memory for Funcfl data
        funcflData.embeddingData[i] = new double[funcflData.numberRhoPoints[i]];
        funcflData.densityData[i] = new double[funcflData.numberRPoints[i]];
        funcflData.ZData[i] = new double[funcflData.numberRPoints[i]];

        ier = ReadFuncflData(pkim, parameterFilePointers[i], i,
                                   funcflData);
        if (ier < KIM_STATUS_OK)
        {
          pkim->report_error(__LINE__, __FILE__, "ReadFuncflData", ier);
          for (int j = 0; j <= i; ++j)
          {
            delete [] funcflData.embeddingData[i];
            delete [] funcflData.densityData[i];
            delete [] funcflData.ZData[i];
          }
          return ier;
        }
      }

      ReinterpolateAndMix(funcflData);
      for (int i = 0; i < numberParameterFiles; ++i)
      {
        delete [] funcflData.embeddingData[i];
        delete [] funcflData.densityData[i];
        delete [] funcflData.ZData[i];
      }

      break;
    }
    default:  // should never get here
    {
      ier = KIM_STATUS_FAIL;
      pkim->report_error(__LINE__, __FILE__, "Bad parameter files(s)", ier);
      return ier;
    }
    break;
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int EAM_Implementation::ReadSetflData(KIM_API_model* const pkim,
                                      FILE* const fptr)
{
  int ier;
  char const* cer;
  char line[MAXLINE];

  // loop over each atom type in the data file
  for (int i = 0; i < numberModelSpecies_; ++i)
  { // read header line (partcle number, mass, lattice constant, lattice type)
    cer = fgets(line, MAXLINE, fptr);
    ier = sscanf(line, "%d %lg %lg %s", &(particleNumber_[i]),
                 &(particleMass_[i]), &(latticeConstant_[i]), latticeType_[i]);
    if ((cer == 0) || (ier != 4))
    {
      ier = KIM_STATUS_FAIL;
      pkim->report_error(__LINE__, __FILE__,
                         "error reading lines of setfl file", ier);
      return ier;
    }

    // read "embed_dat"
    ier = GrabData(pkim, fptr, numberRhoPoints_, embeddingData_[i]);
    if (ier < KIM_STATUS_OK)
    {
      pkim->report_error(__LINE__, __FILE__,
                         "error reading embeddingData lines of setfl"
                         " file", ier);
      return ier;
    }

    // read "densityData"
    ier = GrabData(pkim, fptr, numberRPoints_, densityData_[i][0]);
    if (ier < KIM_STATUS_OK)
    {
      pkim->report_error(__LINE__, __FILE__,
                         "error reading densityData lines of setfl"
                         " file", ier);
      return ier;
    }
    // fill in remaining columns
    for (int j = 1; j < numberModelSpecies_; ++j)
    {
      for (int k = 0; k < numberRPoints_; ++k)
      {
        densityData_[i][j][k] = densityData_[i][0][k];
      }
    }
  }

  // read "rPhiData"
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      ier = GrabData(pkim, fptr, numberRPoints_, rPhiData_[i][j]);
      if (ier < KIM_STATUS_OK)
      {
        pkim->report_error(__LINE__, __FILE__,
                           "error reading rPhiData lines of setfl"
                           " file", ier);
        return ier;
      }
    }
  }

  // filling in upper-triangular part of rPhiData
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    for (int j = i + 1; j < numberModelSpecies_; ++j)
    {
      for (int k = 0; k < numberRPoints_; ++k)
      {
        rPhiData_[i][j][k] = rPhiData_[j][i][k];
      }
    }
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int EAM_Implementation::ReadFinnisSinclairData(
    KIM_API_model* const pkim,
    FILE* const fptr)
{
  int ier;
  char const* cer;
  char line[MAXLINE];

  // loop over each atom species in the data file
  for (int i = 0; i < numberModelSpecies_; ++i)
  { // read header line (partcle number, mass, lattice constant, lattice type)
    cer = fgets(line, MAXLINE, fptr);
    ier = sscanf(line, "%d %lg %lg %s", &(particleNumber_[i]),
                 &(particleMass_[i]), &(latticeConstant_[i]), latticeType_[i]);
    if ((cer == 0) || (ier != 4))
    {
      ier = KIM_STATUS_FAIL;
      pkim->report_error(__LINE__, __FILE__,
                         "error reading lines of setfl file", ier);
      return ier;
    }

    // read "embed_dat"
    ier = GrabData(pkim, fptr, numberRhoPoints_, embeddingData_[i]);
    if (ier < KIM_STATUS_OK)
    {
      pkim->report_error(__LINE__, __FILE__,
                         "error reading embeddingData lines of setfl"
                         " file", ier);
      return ier;
    }

    // read "densityData"
    for (int j = 0; j < numberModelSpecies_; ++j)
    {
      ier = GrabData(pkim, fptr, numberRPoints_, densityData_[i][j]);
      if (ier < KIM_STATUS_OK)
      {
        pkim->report_error(__LINE__, __FILE__,
                           "error reading densityData lines of setfl"
                           " file", ier);
        return ier;
      }
    }
  }

  // read "rPhiData"
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      ier = GrabData(pkim, fptr, numberRPoints_, rPhiData_[i][j]);
      if (ier < KIM_STATUS_OK)
      {
        pkim->report_error(__LINE__, __FILE__,
                           "error reading rPhiData lines of setfl"
                           " file", ier);
        return ier;
      }
    }
  }

  // filling in upper-triangular part of rPhiData
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    for (int j = i + 1; j < numberModelSpecies_; ++j)
    {
      for (int k = 0; k < numberRPoints_; ++k)
      {
        rPhiData_[i][j][k] = rPhiData_[j][i][k];
      }
    }
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int EAM_Implementation::ReadFuncflData(KIM_API_model* const pkim,
                                       FILE* const fptr,
                                       int const fileIndex,
                                       SetOfFuncflData& funcflData)
{
  int ier;

  // read "embed_dat"
  ier = GrabData(pkim, fptr, funcflData.numberRhoPoints[fileIndex],
                 funcflData.embeddingData[fileIndex]);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__,
                       "error reading embeddingData lines of funcfl"
                       " file", ier);
    return ier;
  }

  // read "Z_dat"
  ier = GrabData(pkim, fptr, funcflData.numberRPoints[fileIndex],
                 funcflData.ZData[fileIndex]);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__,
                       "error reading Z_dat lines of funcfl file",
                       ier);
    return ier;
  }

  // read "densityData"
  ier = GrabData(pkim, fptr, funcflData.numberRPoints[fileIndex],
                 funcflData.densityData[fileIndex]);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__,
                       "error reading densityData lines of funcfl file"
                       , ier);
    return ier;
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int EAM_Implementation::GrabData(KIM_API_model* const pkim,
                                 FILE* const fptr, int const n,
                                 double* const list)
{ // This function originally obtained under CDDL from Steve Plimpton
  int ier;
  char const* cer;
  char const* ptr;
  char line[MAXLINE];

  int i = 0;
  while (i < n)
  {
    cer = fgets(line, MAXLINE, fptr);
    if (cer == 0)
    {
      ier = KIM_STATUS_FAIL;
      pkim->report_error(__LINE__, __FILE__,
                         "error reading data from file", ier);
      return ier;
    }

    ptr = strtok(line, " \t\n\r\f");
    list[i] = atof(ptr);
    ++i;
    while ((ptr = strtok(0, " \t\n\r\f")))
    {
      list[i++] = atof(ptr);
    }
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
void EAM_Implementation::ReinterpolateAndMix(
    SetOfFuncflData const& funcflData)
{
  // conversion constants
  double const Hartree = 27.2;  // ev
  double const Bohr = 0.529;    // Angstroms

  // reinterpolate
  double const oneByDrho = ONE / deltaRho_;
  double const oneByDr = ONE / deltaR_;
  if (numberModelSpecies_ > 1)
  {
    for (int i = 0; i < numberModelSpecies_; ++i)
    {
      double* const embeddingCoeff
          = new double[funcflData.numberRhoPoints[i] * NUMBER_SPLINE_COEFF];
      double* const densityCoeff = new double[funcflData.numberRPoints[i] *
                                              NUMBER_SPLINE_COEFF];
      double* const rPhiCoeff = new double[funcflData.numberRPoints[i] *
                                           NUMBER_SPLINE_COEFF];

      SplineInterpolate(funcflData.embeddingData[i],
                        funcflData.deltaRho[i],
                        funcflData.numberRhoPoints[i], embeddingCoeff);
      SplineInterpolate(funcflData.densityData[i], funcflData.deltaR[i],
                        funcflData.numberRPoints[i], densityCoeff);
      SplineInterpolate(funcflData.ZData[i], funcflData.deltaR[i],
                        funcflData.numberRPoints[i], rPhiCoeff);

      for (int j = 0; j < numberRhoPoints_; ++j)
      {
        double densityOffset;
        int densityIndex;
        // compute densityOffset and densityIndex
        double const densityValue = j * deltaRho_;
        GET_DELTAX_AND_INDEX(densityValue, oneByDrho, numberRhoPoints_,
                             densityOffset, densityIndex);
        // interpolate value of embeddingData_[i][j]
        INTERPOLATE_F(embeddingCoeff, densityOffset, densityIndex,
                      embeddingData_[i][j]);
      }

      for (int j = 0; j < numberRPoints_; ++j)
      {
        double rOffset;
        int rIndex;
        double const r = j * deltaR_;
        // compute rOffset and rIndex
        GET_DELTAX_AND_INDEX(r, oneByDr, numberRPoints_, rOffset, rIndex);
        // interpolate value of densityData_[i][j]
        INTERPOLATE_F(densityCoeff, rOffset, rIndex, densityData_[i][0][j]);
        for (int k = 1; k < numberModelSpecies_; ++k)
        {
          densityData_[i][k][j] = densityData_[i][0][j];
        }
        // interpolate value of rPhiData_[i][i][j]
        INTERPOLATE_F(rPhiCoeff, rOffset, rIndex, rPhiData_[i][i][j]);
      }

      delete [] embeddingCoeff;
      delete [] densityCoeff;
      delete [] rPhiCoeff;
    }

    // convert "Z_dat" to r*phi, mix, and store in rPhiData_
    for (int i = 0; i < numberModelSpecies_; ++i)
    {
      for (int j = numberModelSpecies_ - 1; j > i; --j)
      {
        for (int k = 0; k < numberRPoints_; ++k)
          rPhiData_[j][i][k] = rPhiData_[i][j][k]
              = (rPhiData_[i][i][k] * rPhiData_[j][j][k]) * Hartree * Bohr;
      }
      for (int k = 0; k < numberRPoints_; ++k)
        rPhiData_[i][i][k] = (rPhiData_[i][i][k] * rPhiData_[i][i][k])
            * Hartree * Bohr;
    }
  }
  else
  { // if numberModelSpecies_ == 1, don't reinterpolate
    for (int i = 0; i < numberRhoPoints_; ++i)
      embeddingData_[0][i] = funcflData.embeddingData[0][i];
    for (int i = 0; i < numberRPoints_; ++i)
    {
      densityData_[0][0][i] = funcflData.densityData[0][i];
      rPhiData_[0][0][i] = funcflData.ZData[0][i] * funcflData.ZData[0][i]
          * Hartree * Bohr;
    }
  }
}

//******************************************************************************
void EAM_Implementation::CloseParameterFiles(
    FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
    int const numberParameterFiles)
{
  for (int i = 0; i < numberParameterFiles; ++i)
    fclose(parameterFilePointers[i]);
}

//******************************************************************************
int EAM_Implementation::ConvertUnits(KIM_API_model* const pkim)
{
  int ier;

  // define default base units
  char length[] = "A";
  char energy[] = "eV";
  char charge[] = "e";
  char temperature[] = "K";
  char time[] = "ps";

  // changing units of particle mass and lattice constant
  double const
      convertMass = pkim->convert_to_act_unit("m", "kJ/mol", "e", "K", "s",
                                              -2.0, 1.0, 0.0, 0.0, 2.0, &ier);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "convert_to_act_unit", ier);
    return ier;
  }
  double const convertLength
      = pkim->convert_to_act_unit(length, energy, charge, temperature, time,
                                  1.0, 0.0, 0.0, 0.0, 0.0, &ier);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "convert_to_act_unit", ier);
    return ier;
  }
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    particleMass_[i] *= 1.0e-6;            // convert to (kJ/mole)*(s^2)/(m^2)
    particleMass_[i] *= convertMass;       // convert to active units
    latticeConstant_[i] *= convertLength;  // convert to active units
  }

  // changing units of embedding function values
  // don't convert the density units (argument of embedding function)
  double const convertEnergy
      = pkim->convert_to_act_unit(length, energy, charge, temperature, time,
                                  0.0, 1.0, 0.0, 0.0, 0.0, &ier);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "convert_to_act_unit", ier);
    return ier;
  }
  if (convertEnergy != ONE)
  {
    for (int i = 0; i < numberModelSpecies_; ++i)
    {
      for (int j = 0; j < numberRhoPoints_; ++j)
      {
        embeddingData_[i][j] *= convertEnergy;
      }
    }
  }

  //
  // don't convert units of density (rho) --- they are ambiguous
  //

  // changing units of r*phi (stored in rPhiData) function values
  double const convertRPhi
      = pkim->convert_to_act_unit(length, energy, charge, temperature, time,
                                  1.0, 1.0, 0.0, 0.0, 0.0, &ier);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "convert_to_act_unit", ier);
    return ier;
  }
  if (convertRPhi != ONE)
  {
    for (int i = 0; i < numberModelSpecies_; ++i)
    {
      for (int j = 0; j < numberModelSpecies_; ++j)
      {
        for (int k = 0; k < numberRPoints_; ++k)
        {
          rPhiData_[i][j][k] *= convertRPhi;
        }
      }
    }
  }

  // changing units of cutoff radius and deltaR
  if (convertLength != ONE)
  {
    cutoffParameter_ *= convertLength;
    deltaR_ *= convertLength;

    //
    // don't convert units of deltaRho --- they are ambiguous
    //
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int EAM_Implementation::RegisterKIMParameters(
    KIM_API_model* const pkim,
    EAMFileType const eamFileType) const
{
  int ier;

  // publish parameters
  int const numberCommentLines = (eamFileType == Funcfl) ?
      numberModelSpecies_ : NUMBER_SETFL_COMMENT_LINES;
  pkim->setm_data(&ier, 11 * 4,
                  "PARAM_FIXED_comments",
                  numberCommentLines, (void*)comments_ptr_, 1,
                  "PARAM_FIXED_particleNames",
                  1, (void*)particleNames_, 1,
                  "PARAM_FIXED_particleNumber",
                  numberModelSpecies_, (void*)particleNumber_, 1,
                  "PARAM_FIXED_particleMass",
                  numberModelSpecies_, (void*)particleMass_, 1,
                  "PARAM_FIXED_latticeConstant",
                  numberModelSpecies_, (void*)latticeConstant_, 1,
                  "PARAM_FIXED_latticeType",
                  numberModelSpecies_, (void*)latticeType_, 1,
                  "PARAM_FIXED_numberRhoPoints",
                  1, (void*)&(numberRhoPoints_), 1,
                  "PARAM_FIXED_numberRPoints",
                  1, (void*)&(numberRPoints_), 1,
                  "PARAM_FREE_cutoff",
                  1, (void*)&(cutoffParameter_), 1,
                  "PARAM_FREE_deltaRho",
                  1, (void*)&(deltaRho_), 1,
                  "PARAM_FREE_deltaR",
                  1, (void*)&(deltaR_), 1);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "setm_data", ier);
    return ier;
  }

  for (int i = 0; i < numberModelSpecies_ ; i++)
  {
    for(int j = i; j < numberModelSpecies_; j++)
    {
      int indxPhi = i*numberModelSpecies_ + j - (i*i + i)/2;
      for(int k = 0; k < numberRPoints_; k++)
      {
        publish_rPhiData_[indxPhi][k] = rPhiData_[i][j][k];
      }
    }
  }

  int publishDensityLength = (eamFileType_ == FinnisSinclair) ?
      numberModelSpecies_*numberModelSpecies_ :
      numberModelSpecies_;

  for (int i = 0; i < numberModelSpecies_ ; i++)
  {
    for(int j = 0; j < numberModelSpecies_; j++)
    {
      int indxDensity = (eamFileType_ == FinnisSinclair) ?
          i*numberModelSpecies_ + j : i;
      for(int k = 0; k < numberRPoints_; k++)
      {
        publishDensityData_[indxDensity][k] = densityData_[i][j][k];
      }
      if (eamFileType_ != FinnisSinclair) break;
    }
  }
  int shape_embed[2] = {numberModelSpecies_, numberRhoPoints_};
  int shape_r_phi[2]
      = {numberUniqueSpeciesPairs_, numberRPoints_};
  int shape_density[2]
      = {publishDensityLength, numberRPoints_};
  pkim->set_shape("PARAM_FREE_embeddingData", shape_embed, 2, &ier);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "set_shape", ier);
    return ier;
  }
  pkim->set_shape("PARAM_FREE_densityData", shape_density, 2, &ier);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "set_shape", ier);
    return ier;
  }
  pkim->set_shape("PARAM_FREE_rPhiData", shape_r_phi, 2, &ier);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "set_shape", ier);
    return ier;
  }
  pkim->setm_data(&ier, 3 * 4,
                  "PARAM_FREE_embeddingData",
                  numberModelSpecies_ * numberRhoPoints_,
                  (void*) embeddingData_[0],
                  1,
                  //
                  "PARAM_FREE_densityData",
                  publishDensityLength * numberRPoints_,
                  (void*) publishDensityData_[0],
                  1,
                  //
                  "PARAM_FREE_rPhiData",
                  numberUniqueSpeciesPairs_ * numberRPoints_,
                  (void*) publish_rPhiData_[0],
                  1);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "setm_data", ier);
    return ier;
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int EAM_Implementation::RegisterKIMFunctions(KIM_API_model* const pkim)
    const
{
  int ier;

  // register the destroy() and reinit() functions
  pkim->setm_method(&ier, 3 * 4,
                    "destroy", 1, (func_ptr) &(EAM::Destroy), 1,
                    "reinit",  1, (func_ptr) &(EAM::Reinit),  1,
                    "compute", 1, (func_ptr) &(EAM::Compute), 1);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "setm_method", ier);
    return ier;
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int EAM_Implementation::SetReinitMutableValues(KIM_API_model* const pkim)
{ // use (possibly) new values of free parameters to compute other quantities
  int ier;

  // get cutoff pointer
  double* const cutoff
      = static_cast<double*>(pkim->get_data_by_index(cutoffIndex_, &ier));
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "get_data_by_index", ier);
    return ier;
  }

  // update cutoff value in KIM API object
  *cutoff = cutoffParameter_;

  // update EAM_Implementation values
  cutoffSq_ = cutoffParameter_ * cutoffParameter_;
  oneByDr_ = ONE / deltaR_;
  oneByDrho_ = ONE / deltaRho_;

  // calculate spline interpolating coefficients
  SplineInterpolateAllData();

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
void EAM_Implementation::SplineInterpolateAllData()
{
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    SplineInterpolate(embeddingData_[i], deltaRho_, numberRhoPoints_,
                      embeddingCoeff_[i]);
    for (int j = 0; j < numberModelSpecies_; ++j)
    {
      SplineInterpolate(densityData_[i][j], deltaR_, numberRPoints_,
                        densityCoeff_[i][j]);
      SplineInterpolate(rPhiData_[i][j], deltaR_, numberRPoints_,
                        rPhiCoeff_[i][j]);
    }
  }
}

//******************************************************************************
int EAM_Implementation::SetComputeMutableValues(
    KIM_API_model* const pkim,
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
    VectorOfSizeDIM*& forces)
{
  int ier = KIM_STATUS_FAIL;

  // get compute flags
  int compEnergy;
  int compForces;
  int compParticleEnergy;
  int compProcess_dEdr;
  int compProcess_d2Edr2;
  pkim->getm_compute_by_index(&ier, 3 * 5,
                              energyIndex_,         &compEnergy,         1,
                              forcesIndex_,         &compForces,         1,
                              particleEnergyIndex_, &compParticleEnergy, 1,
                              process_dEdrIndex_,   &compProcess_dEdr,   1,
                              process_d2Edr2Index_, &compProcess_d2Edr2, 1);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "getm_compute_by_index", ier);
    return ier;
  }

  isComputeEnergy = (compEnergy == KIM_COMPUTE_TRUE);
  isComputeForces = (compForces == KIM_COMPUTE_TRUE);
  isComputeParticleEnergy = (compParticleEnergy == KIM_COMPUTE_TRUE);
  isComputeProcess_dEdr = (compProcess_dEdr == KIM_COMPUTE_TRUE);
  isComputeProcess_d2Edr2 = (compProcess_d2Edr2 == KIM_COMPUTE_TRUE);

  // extract pointers based on compute flags
  //
  // double const* cutoff;            // currently unused
  // int const* numberOfSpecies;  // currently unused
  int const* numberOfParticles;
  int const* numberContributingParticles;
  pkim->getm_data_by_index(
      &ier, 3 * 8,
      // cutoffIndex_, &cutoff, 1,
      // numberOfSpeciesIndex_, &numberOfSpecies, 1,
      numberOfParticlesIndex_, &numberOfParticles, 1,
      numberContributingParticlesIndex_, &numberContributingParticles, isHalf_,
      particleSpeciesIndex_, &particleSpecies, 1,
      boxSideLengthsIndex_, &boxSideLengths, (NBCType_ == Mi_Opbc),
      coordinatesIndex_, &coordinates, 1,
      energyIndex_, &energy, compEnergy,
      particleEnergyIndex_, &particleEnergy, compParticleEnergy,
      forcesIndex_, &forces, compForces);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "getm_data_by_index", ier);
    return ier;
  }
  if (NBCType_ != Cluster)
  {
    get_neigh = (GetNeighborFunction *)
        pkim->get_method_by_index(get_neighIndex_, &ier);
    if (ier < KIM_STATUS_OK)
    {
      pkim->report_error(__LINE__, __FILE__, "get_method_by_index", ier);
      return ier;
    }
  }

  // allocate memory if needed
  if (cachedNumberOfParticles_ < *numberOfParticles)
  {
    delete [] densityValue_;  // ok to delete null pointer
    densityValue_ = new double[*numberOfParticles];

    delete [] embeddingDerivativeValue_;  // ok to delete null pointer
    embeddingDerivativeValue_ = new double[*numberOfParticles];

    delete [] embeddingSecondDerivativeValue_;  // ok to delete null pointer
    embeddingSecondDerivativeValue_ = new double[*numberOfParticles];
  }

  // update values
  cachedNumberOfParticles_ = *numberOfParticles;

  // set cachedNumberContributingParticles based on half/full neighbor list
  if ((isHalf_) && (NBCType_ != Cluster))
  {
    cachedNumberContributingParticles_ = *numberContributingParticles;
  }
  else
  { // set so that it can be used even with a full neighbor list
    cachedNumberContributingParticles_ = *numberOfParticles;
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int EAM_Implementation::CheckParticleSpecies(
    KIM_API_model* const pkim,
    int const* const particleSpecies)
    const
{
  int ier;
  for (int i = 0; i < cachedNumberOfParticles_; ++i)
  {
    if ((particleSpecies[i] < 0) || (particleSpecies[i] >= numberModelSpecies_))
    {
      ier = KIM_STATUS_FAIL;
      pkim->report_error(__LINE__, __FILE__,
                         "unsupported particle species detected", ier);
      return ier;
    }
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int EAM_Implementation::GetComputeIndex(
    const bool& isComputeProcess_dEdr,
    const bool& isComputeProcess_d2Edr2,
    const bool& isComputeEnergy,
    const bool& isComputeForces,
    const bool& isComputeParticleEnergy) const
{
  // const int iter = 3;
  const int half = 2;
  const int rij = 3;
  const int processdE = 2;
  const int processd2E = 2;
  const int energy = 2;
  const int force = 2;
  const int particleEnergy = 2;


  int index = 0;

  // iter
  //  Cluster  = 0
  //        (for EAM when NBCType_==Cluster then isLocatorMode_==true)
  //  Locator  = 1
  //  Iterator = 2
  index += (int(NBCType_ != Cluster) + int(!isLocatorMode_))
      * half * rij * processdE * processd2E * energy * force * particleEnergy;

  // half
  index += (int(isHalf_))
      * rij * processdE * processd2E * energy * force * particleEnergy;

  // rij
  //  Coordinates = 0
  //       (NBCType_ != Neigh_Rvec and != Mi_OPbc)
  //  RVec        = 1
  //  Mi_Opbc     = 2
  index += (int(NBCType_ == Neigh_Rvec) + 2*int(NBCType_ == Mi_Opbc))
      * processdE * processd2E * energy * force * particleEnergy;

  // processdE
  index += (int(isComputeProcess_dEdr))
      * processd2E * energy * force * particleEnergy;

  // processd2E
  index += (int(isComputeProcess_d2Edr2))
      * energy * force * particleEnergy;

  // energy
  index += (int(isComputeEnergy))
      * force * particleEnergy;

  // force
  index += (int(isComputeForces))
      * particleEnergy;

  // particleEnergy
  index += (int(isComputeParticleEnergy));

  return index;
}

//******************************************************************************
void EAM_Implementation::ApplyMIOPBC(double const* const boxSideLengths,
                                     double* const dx)
{
  double sign;
  for (int i = 0; i < DIMENSION; ++i)
  {
    sign = dx[i] > 0 ? ONE : -ONE;
    dx[i] = (abs(dx[i]) > HALF * boxSideLengths[i]) ? dx[i]
        - sign * boxSideLengths[i] : dx[i];
  }
}

//==============================================================================
//
// Implementation of helper functions
//
//==============================================================================

//******************************************************************************
void AllocateAndInitialize2DArray(double**& arrayPtr, int const extentZero,
                                  int const extentOne)
{ // allocate memory and set pointers
  arrayPtr = new double*[extentZero];
  arrayPtr[0] = new double[extentZero * extentOne];
  for (int i = 1; i < extentZero; ++i)
  {
    arrayPtr[i] = arrayPtr[i-1] + extentOne;
  }

  // initialize
  for (int i = 0; i < extentZero; ++i)
  {
    for (int j = 0; j < extentOne; ++j)
    {
      arrayPtr[i][j] = 0.0;
    }
  }
}

//******************************************************************************
void Deallocate2DArray(double**& arrayPtr)
{ // deallocate memory
  if (arrayPtr != 0) delete [] arrayPtr[0];
  delete [] arrayPtr;

  // nullify pointer
  arrayPtr = 0;
}

//******************************************************************************
void AllocateAndInitialize3DArray(double***& arrayPtr, int const extentZero,
                                  int const extentOne, int const extentTwo)
{ // allocate memory and set pointers
  arrayPtr = new double**[extentZero];
  arrayPtr[0] = new double*[extentZero * extentOne];
  arrayPtr[0][0] = new double[extentZero * extentOne * extentTwo];
  for (int j = 1; j < extentZero; ++j)
  {
    arrayPtr[j] = arrayPtr[j-1] + extentOne;
    arrayPtr[0][j] = arrayPtr[0][j-1] + extentTwo;
  }
  for (int i = 1; i < extentZero; ++i)
  {
    arrayPtr[i][0] = arrayPtr[i-1][extentOne-1] + extentTwo;
    for (int j = 1; j < extentOne; ++j)
    {
      arrayPtr[i][j] = arrayPtr[i][j-1] + extentTwo;
    }
  }

  // initialize
  for (int i = 0; i < extentZero; ++i)
  {
    for (int j = 0; j < extentOne; ++j)
    {
      for (int k = 0; k < extentTwo; ++k)
      {
        arrayPtr[i][j][k] = 0.0;
      }
    }
  }
}

//******************************************************************************
void Deallocate3DArray(double***& arrayPtr)
{ // deallocate memory
  if (arrayPtr != 0)
  {
    delete [] arrayPtr[0][0];
    delete [] arrayPtr[0];
  }
  delete [] arrayPtr;

  // nullify pointer
  arrayPtr = 0;
}
