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
#    Mingjian Wen


# load all basic KIM make configuration
include ../Makefile.KIM_Config

# set model driver specific details
MODEL_DRIVER_NAME := Three_Body_Stillinger_Weber__MD_335816936951_002
MODEL_NAME := Three_Body_Stillinger_Weber_CdTeZnSeHgS__MO_503261197030_000
SPECIES_001_NAME    := Cd
SPECIES_002_NAME    := Te
SPECIES_003_NAME    := Zn
SPECIES_004_NAME    := Se
SPECIES_005_NAME    := Hg
SPECIES_006_NAME    := S
PARAM_FILE_001_NAME := Three_Body_Stillinger_Weber_CdTeZnSeHgS.params

# APPEND to compiler option flag lists
#FFLAGS   +=
#CFLAGS   +=
#CXXFLAGS +=
#LDFLAGS  +=

# load remaining KIM make configuration
include $(KIM_DIR)/$(builddir)/Makefile.ParameterizedModel
