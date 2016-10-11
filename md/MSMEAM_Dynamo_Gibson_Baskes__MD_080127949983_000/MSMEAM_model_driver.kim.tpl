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
# Copyright (c) 2013--2014, Regents of the University of Minnesota.
# All rights reserved.
#
# Contributors:
#    Ryan S. Elliott
#    Ellad B. Tadmor
#    Valeriu Smirichinski
#    Stephen M. Whalen
#
#
#
# Portions Copyright (c) 2015, Joshua S. Gibson.
# All rights reserved.
#
# Contributors:
#    J.S. Gibson
#    S.G  Srivilliputhur
#    M.I. Baskes
#    A.K. Wilson
#    R.E. Miller
#

#######################################################################################################
#
# See src/standard.kim for documentation about this file
#
#######################################################################################################


KIM_API_Version := 1.6.3

Unit_Handling    := fixed
Unit_length      := A
Unit_energy      := eV
Unit_charge      := e
Unit_temperature := K
Unit_time        := ps


#######################################################################################################
PARTICLE_SPECIES:
# Symbol/name               Type                    code

SPECIES_001_NAME_STR        spec                    1


#######################################################################################################
CONVENTIONS:
# Name                      Type

OneBasedLists               flag

Neigh_IterAccess            flag

Neigh_LocaAccess            flag

NEIGH_PURE_F                flag

#######################################################################################################
MODEL_INPUT:
# Name                      Type         Unit                Shape              Requirements

numberOfParticles           integer      none                []

numberContributingParticles integer      none                []                 optional

numberOfSpecies             integer      none                []

particleSpecies             integer      none                [numberOfParticles]

coordinates                 double       length              [numberOfParticles,3]

boxSideLengths              double       length              [3]                optional

get_neigh                   method       none                []                 optional

neighObject                 pointer      none                []                 optional

#######################################################################################################
MODEL_OUTPUT:
# Name                      Type         Unit                Shape              Requirements

destroy                     method       none                []

compute                     method       none                []

reinit                      method       none                []                 optional

cutoff                      double       length              []

energy                      double       energy              []                 optional

forces                      double       force               [numberOfParticles,3]  optional

particleEnergy              double       energy              [numberOfParticles]    optional

particleVirial              double       energy              [numberOfParticles,6]  optional

virial                      double       energy              [6]                    optional

#######################################################################################################
MODEL_PARAMETERS:
# Name                      Type         Unit                Shape              Requirements

# --- No published parameters ---
