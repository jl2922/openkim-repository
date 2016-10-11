# load all basic KIM make configuration
ifeq ($(wildcard ../Makefile.KIM_Config),)
  $(error ../Makefile.KIM_Config does not exist.  Something is wrong with your KIM API package setup)
endif
include ../Makefile.KIM_Config

# set model driver specific details
MODEL_DRIVER_NAME   := EAM_Dynamo__MD_120291908751_002
MODEL_NAME          := EAM_Dynamo_Foiles_Baskes_Universal3_Cu__MO_666348409573_001
SPECIES_001_NAME    := Cu
PARAM_FILE_001_NAME := Cu_u3.eam

# load remaining KIM make configuration
include $(KIM_DIR)/$(builddir)/Makefile.ParameterizedModel
