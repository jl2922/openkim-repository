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
#    Nikhil Chandra Admal
#    Ellad B. Tadmor
#



# load all basic KIM make configuration
ifeq ($(wildcard ../Makefile.KIM_Config),)
  $(error ../Makefile.KIM_Config does not exist.  Something is wrong with your KIM API package setup)
endif
include ../Makefile.KIM_Config

# set model driver specific details
MODEL_DRIVER_NAME := Pair_LJ_Smooth__MD_716364606252_000
MODEL_NAME        := Pair_LJ_Smooth_Bernardes_Ar__MO_764178710049_000
SPECIES_001_NAME  := Ar
PARAM_FILE_001_NAME := Pair_LJ_Smooth.params

# APPEND to compiler option flag lists
#FFLAGS   +=
#CCFLAGS  +=
#CXXFLAGS +=
#LDFLAGS  +=

# load remaining KIM make configuration
include $(KIM_DIR)/$(builddir)/Makefile.ParameterizedModel
