#
# CDDL HEADER START
#
# The contents of this file are subject to the terms of the Common Development
# and Distribution License Version 1.0 (the "License").
#
# You can obtain a copy of the license at
# http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
# specific language governing permissions and limitations under the License.
#
# When distributing Covered Code, include this CDDL HEADER in each file and
# include the License file in a prominent location with the name LICENSE.CDDL.
# If applicable, add the following below this CDDL HEADER, with the fields
# enclosed by brackets "[]" replaced with your own identifying information:
#
# Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
#
# CDDL HEADER END
#

#
# Copyright (c) 2012, Regents of the University of Minnesota.  All rights reserved.
#
# Contributors:
#    Ryan S. Elliott
#    Ellad B. Tadmor
#    Valeriu Smirichinski
#


KIM_API_Version := 1.6.0

Unit_Handling    := fixed
Unit_length      := A
Unit_energy      := eV
Unit_charge      := e
Unit_temperature := K
Unit_time        := ps


#######################################################################################################
PARTICLE_SPECIES:
# Symbol/name               Type                    code

SPECIES_001_NAME_STR         spec                    -1
SPECIES_002_NAME_STR         spec                    -1
SPECIES_003_NAME_STR         spec                    -1
SPECIES_004_NAME_STR         spec                    -1
SPECIES_005_NAME_STR         spec                    -1
SPECIES_006_NAME_STR         spec                    -1
SPECIES_007_NAME_STR         spec                    -1
SPECIES_008_NAME_STR         spec                    -1
SPECIES_009_NAME_STR         spec                    -1
SPECIES_010_NAME_STR         spec                    -1
SPECIES_011_NAME_STR         spec                    -1
SPECIES_012_NAME_STR         spec                    -1
SPECIES_013_NAME_STR         spec                    -1
SPECIES_014_NAME_STR         spec                    -1
SPECIES_015_NAME_STR         spec                    -1
SPECIES_016_NAME_STR         spec                    -1
SPECIES_017_NAME_STR         spec                    -1
SPECIES_018_NAME_STR         spec                    -1
SPECIES_019_NAME_STR         spec                    -1
SPECIES_020_NAME_STR         spec                    -1
SPECIES_021_NAME_STR         spec                    -1
SPECIES_022_NAME_STR         spec                    -1
SPECIES_023_NAME_STR         spec                    -1
SPECIES_024_NAME_STR         spec                    -1
SPECIES_025_NAME_STR         spec                    -1
SPECIES_026_NAME_STR         spec                    -1
SPECIES_027_NAME_STR         spec                    -1
SPECIES_028_NAME_STR         spec                    -1
SPECIES_029_NAME_STR         spec                    -1
SPECIES_030_NAME_STR         spec                    -1
SPECIES_031_NAME_STR         spec                    -1
SPECIES_032_NAME_STR         spec                    -1
SPECIES_033_NAME_STR         spec                    -1
SPECIES_034_NAME_STR         spec                    -1
SPECIES_035_NAME_STR         spec                    -1
SPECIES_036_NAME_STR         spec                    -1
SPECIES_037_NAME_STR         spec                    -1
SPECIES_038_NAME_STR         spec                    -1
SPECIES_039_NAME_STR         spec                    -1
SPECIES_040_NAME_STR         spec                    -1
SPECIES_041_NAME_STR         spec                    -1
SPECIES_042_NAME_STR         spec                    -1
SPECIES_043_NAME_STR         spec                    -1
SPECIES_044_NAME_STR         spec                    -1
SPECIES_045_NAME_STR         spec                    -1
SPECIES_046_NAME_STR         spec                    -1
SPECIES_047_NAME_STR         spec                    -1
SPECIES_048_NAME_STR         spec                    -1
SPECIES_049_NAME_STR         spec                    -1
SPECIES_050_NAME_STR         spec                    -1


#######################################################################################################
CONVENTIONS:
# Name                      Type

OneBasedLists               flag

Neigh_LocaAccess            flag

MI_OPBC_F                   flag


#######################################################################################################
MODEL_INPUT:
# Name                      Type         Unit                Shape              Requirements

numberOfParticles           integer      none                []

numberContributingParticles integer      none                []                 optional

numberOfSpecies         integer      none                []

particleSpecies               integer      none                [numberOfParticles]

coordinates                 double       length              [numberOfParticles,3]

boxSideLengths              double       length              [3]                optional

get_neigh                   method       none                []                 optional

neighObject                 pointer      none                []                 optional


#######################################################################################################
MODEL_OUTPUT:
# Name                      Type         Unit                Shape              Requirements

destroy                     method       none                []                 optional

compute                     method       none                []

reinit                      method       none                []                 optional

cutoff                      double       length              []

energy                      double       energy              []                 optional

forces                      double       force               [numberOfParticles,3]  optional

particleEnergy              double       energy              [numberOfParticles]    optional

virial                      double       energy              [6]                optional


#######################################################################################################
MODEL_PARAMETERS:
# Name                      Type         Unit                Shape              Requirements

# --- No published parameters ---
