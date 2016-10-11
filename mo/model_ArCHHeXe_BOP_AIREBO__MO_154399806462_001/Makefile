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


# load all basic KIM make configuration
include ../Makefile.KIM_Config

# set model specific details
MODEL_NAME := model_ArCHHeXe_BOP_AIREBO__MO_154399806462_001
MODEL_KIM_FILE := model_ArCHHeXe_BOP_AIREBO.kim
MODEL_INIT_FUNCTION_NAME := model_init

MODELSRCS = blkdat.f blkdtw.f caguts.f calcforce.f dgesv.f die.f dihedrals.f \
            init.f initpots.f Jstuff.f ljguts.f ljparam.f \
            mathstuff.f mtable.f pairlist.f pairwise.f param.f \
            pibond.f qstuff.f screen.f setin.f splines.f spoof.f

LOCALOBJ = model_ArCHHeXe_BOP_AIREBO.o ${MODELSRCS:.f=.o}

LOCALCLEAN = dgesv.f blkdtw.f flush.f perror.f configured

all: configured

configured:
	./configure_kim "$(FC)" "$(FFLAGS)"

blkdtw.f: configured
dgesv.f: configured
flush.f: configured
perror.f: configured


# APPEND to compiler option flag lists
#FFLAGS   +=
#CFLAGS   +=
#CXXFLAGS +=
#LDFLAGS  +=

# load remaining KIM make configuration
include $(KIM_DIR)/$(builddir)/Makefile.Model
