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
# Copyright (c) 2012,   Institute for Theoretical and Applied Physics
#      			University of Stuttgart, D-70550 Stuttgart, Germany.
# 			All rights reserved.

#
# Contributors:
#    Daniel Schopf
#

# load all basic KIM make configuration
include ../Makefile.KIM_Config

# set model driver specific details
MODEL_DRIVER_NAME := IMD_EAM__MD_113599595631_002
MODEL_NAME := IMD_EAM_Schopf_CaCd__MO_145183423516_002
SPECIES_001_NAME    := Ca
SPECIES_002_NAME    := Cd
PARAM_FILE_001_NAME := cacd.222_gga.23170_new_phi.imd.pt
PARAM_FILE_002_NAME := cacd.222_gga.23170_new_rho.imd.pt
PARAM_FILE_003_NAME := cacd.222_gga.23170_new_F.imd.pt

# APPEND to compiler option flag lists
#FFLAGS   +=
#CFLAGS   +=
#CXXFLAGS +=
#LDFLAGS  +=

# load remaining KIM make configuration
include $(KIM_DIR)/$(builddir)/Makefile.ParameterizedModel
