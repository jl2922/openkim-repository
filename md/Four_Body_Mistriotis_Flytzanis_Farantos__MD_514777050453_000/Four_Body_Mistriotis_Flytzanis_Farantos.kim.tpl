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
#    Amit Singh

KIM_API_Version := 1.7.2

Unit_Handling    := flexible
Unit_length      := A
Unit_energy      := eV
Unit_charge      := e
Unit_temperature := K
Unit_time        := ps

#######################################################################################################
SUPPORTED_ATOM/PARTICLES_TYPES:
# Symbol/name               Type                    code

SPECIES_001_NAME_STR        spec                    1
SPECIES_002_NAME_STR        spec                    2

#######################################################################################################
CONVENTIONS:
# Name                      Type

ZeroBasedLists              flag

Neigh_IterAccess            flag

Neigh_LocaAccess            flag

NEIGH_PURE_F                flag

NEIGH_RVEC_F                flag

MI_OPBC_F                   flag

CLUSTER                     flag

#######################################################################################################
MODEL_INPUT:
# Name                      Type         Unit                Shape              Requirements

numberOfParticles           integer      none                []

numberContributingParticles integer      none                []                  optional

numberOfSpecies              integer      none                []

particleSpecies               integer      none                [numberOfParticles]

coordinates                 double       length              [numberOfParticles,3]

boxSideLengths              double       length              [3]                  optional

get_neigh                   method       none                []                   optional

neighObject                 pointer      none                []                   optional

#######################################################################################################
MODEL_OUTPUT:
# Name                      Type         Unit                Shape              Requirements

destroy                     method       none                []

compute                     method       none                []

reinit                      method       none                []                    optional

cutoff                      double       length              []

energy                      double       energy              []                    optional

forces                      double       force               [numberOfParticles,3] optional

particleEnergy              double       energy              [numberOfParticles]   optional

process_dEdr                method       none                []                    optional

process_d2Edr2              method       none                []                    optional

#######################################################################################################
MODEL_PARAMETERS:
# Name                      Type         Unit                Shape              Requirements

PARAM_FREE_cutoff           double       length                 []                optional

PARAM_FREE_A                double       none                   []

PARAM_FREE_B                double       none                   []

PARAM_FREE_p                double       none                   []

PARAM_FREE_q                double       none                   []

PARAM_FREE_a                double       none                   []

PARAM_FREE_lambda           double       none                   []

PARAM_FREE_lambda_2           double       none                   []

PARAM_FREE_gamma            double       none                   []

PARAM_FREE_sigma            double       length                 []

PARAM_FREE_epsilon          double       energy                 []

PARAM_FREE_Q                double       none                 []

PARAM_FREE_costhetat        double       none                   []

PARAM_FIXED_cutsq           double       length                 []

