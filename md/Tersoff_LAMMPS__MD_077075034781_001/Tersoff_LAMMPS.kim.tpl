# -*- conf -*-
#
# Copyright (c) 2012,2013 Tobias Brink
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


KIM_API_Version := 1.6.0

## Units ###############################################################
Unit_Handling    := flexible
Unit_length      := A
Unit_energy      := eV
Unit_charge      := e
Unit_temperature := K
Unit_time        := ps

## Supported atom types ################################################
# Add more species if needed, no other changes needed.  A model which
# uses this driver MUST have continuous numbering for the species,
# i.e. the model driver expects the righmost column to have all
# numbers in [0;n_spec).  Failure to do this WILL lead to an error.
PARTICLE_SPECIES:
# Symbol/name               Type                    code
SPECIES_001_NAME_STR        spec                    0
SPECIES_002_NAME_STR        spec                    1
SPECIES_003_NAME_STR        spec                    2
SPECIES_004_NAME_STR        spec                    3
SPECIES_005_NAME_STR        spec                    4
SPECIES_006_NAME_STR        spec                    5
SPECIES_007_NAME_STR        spec                    6
SPECIES_008_NAME_STR        spec                    7
SPECIES_009_NAME_STR        spec                    8
SPECIES_010_NAME_STR        spec                    9
SPECIES_011_NAME_STR        spec                   10
SPECIES_012_NAME_STR        spec                   11
SPECIES_013_NAME_STR        spec                   12
SPECIES_014_NAME_STR        spec                   13
SPECIES_015_NAME_STR        spec                   14
SPECIES_016_NAME_STR        spec                   15

## Supported neighbor list/iteration mode ##############################
CONVENTIONS:

ZeroBasedLists              flag
Neigh_IterAccess            flag
Neigh_LocaAccess            flag
NEIGH_RVEC_F                flag
NEIGH_PURE_F                flag
MI_OPBC_F                   flag
CLUSTER                     flag

# Cannot support half lists due to functional form which would then
# need neighbors of i and j at once in an inner loop.  This is not
# possible to implement efficiently.


## Input spec ##########################################################
MODEL_INPUT:

# variable           type      unit      dimensions            requirements
numberOfParticles    integer   none      []
numberOfSpecies  integer   none      []
particleSpecies        integer   none      [numberOfParticles]
coordinates          double    length    [numberOfParticles,3]
get_neigh            method    none      []                    optional
neighObject          pointer   none      []                    optional
boxSideLengths       double    length    [3]                   optional
# This one is unused in the model.
numberContributingParticles integer none []                    optional

## Output spec #########################################################
MODEL_OUTPUT:

# variable           type      unit      dimensions            requirements
destroy              method    none      []
compute              method    none      []
#reinit              method    none      []                    optional
cutoff               double    length    []
energy               double    energy    []                    optional
forces               double    force     [numberOfParticles,3] optional
particleEnergy       double    energy    [numberOfParticles]   optional
virial               double    energy    [6]                   optional
particleVirial       double    energy    [numberOfParticles,6] optional


## Parameters ##########################################################
MODEL_PARAMETERS:

# variable           type      unit      dimensions
#PARAM_FREE_A        double    energy    [numberOfSpecies,numberOfSpecies]
#PARAM_FREE_B        double    energy    [numberOfSpecies,numberOfSpecies]
#PARAM_FREE_lambda1  double    length^-1 [numberOfSpecies,numberOfSpecies]
#PARAM_FREE_lambda2  double    length^-1 [numberOfSpecies,numberOfSpecies]
#PARAM_FREE_lambda3  double    length^-1 [numberOfSpecies,numberOfSpecies]
#PARAM_FREE_beta     double    none      [numberOfSpecies,numberOfSpecies]
#PARAM_FREE_n        double    none      [numberOfSpecies,numberOfSpecies]
#PARAM_FREE_m        integer   none      [numberOfSpecies,numberOfSpecies]
#PARAM_FREE_gamma    double    none      [numberOfSpecies,numberOfSpecies]
#PARAM_FREE_c        double    none      [numberOfSpecies,numberOfSpecies]
#PARAM_FREE_d        double    none      [numberOfSpecies,numberOfSpecies]
#PARAM_FREE_h        double    none      [numberOfSpecies,numberOfSpecies]
#PARAM_FREE_Rc       double    length    [numberOfSpecies,numberOfSpecies]
#PARAM_FREE_Dc       double    length    [numberOfSpecies,numberOfSpecies]

# Maximum cutoff = Rc + Dc.
#PARAM_FIXED_cut     double    length    [numberOfSpecies,numberOfSpecies]
# Square of maximum cutoff.
#PARAM_FIXED_cutsq   double    length^2  [numberOfSpecies,numberOfSpecies]
# Value of cutoff function is 1.0 below this value, i.e. Rc - Dc.
#PARAM_FIXED_cutmin  double    length    [numberOfSpecies,numberOfSpecies]
# These are used to check if certain numerical approximations may be used.
#PARAM_FIXED_n_precomp double  none      [numberOfSpecies,numberOfSpecies,4]
# c^2
#PARAM_FIXED_c2      double    none      [numberOfSpecies,numberOfSpecies]
# d^2
#PARAM_FIXED_d2      double    none      [numberOfSpecies,numberOfSpecies]
# c^2 / d^2
#PARAM_FIXED_c2_d2   double    none      [numberOfSpecies,numberOfSpecies]
