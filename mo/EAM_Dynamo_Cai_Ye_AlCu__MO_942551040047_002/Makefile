# Copyright 2013 Sang-Pil Kim
# 
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

# load all basic KIM make configuration
include ../Makefile.KIM_Config

# set model driver specific details
MODEL_DRIVER_NAME := EAM_Dynamo__MD_120291908751_002
MODEL_NAME := EAM_Dynamo_Cai_Ye_AlCu__MO_942551040047_002
SPECIES_001_NAME    := Al
SPECIES_002_NAME    := Cu
PARAM_FILE_001_NAME := AlCu.eam.alloy

# APPEND to compiler option flag lists
#FFLAGS   +=
#CFLAGS   +=
#CXXFLAGS +=
#LDFLAGS  +=

# load remaining KIM make configuration
include $(KIM_DIR)/$(builddir)/Makefile.ParameterizedModel
