# load all basic KIM make configuration
ifeq ($(wildcard ../Makefile.KIM_Config),)
  $(error ../Makefile.KIM_Config does not exist.  Something is wrong with your KIM API package setup)
endif
include ../Makefile.KIM_Config

# set model driver specific details
MODEL_DRIVER_NAME   := Four_Body_Mistriotis_Flytzanis_Farantos__MD_514777050453_000
MODEL_NAME          := Four_Body_Mistriotis_Flytzanis_Farantos_Si__MO_080526771943_000
SPECIES_001_NAME    := Si
PARAM_FILE_001_NAME := Four_Body_Mistriotis_Flytzanis_Farantos_Si.params

# load remaining KIM make configuration
include $(KIM_DIR)/$(builddir)/Makefile.ParameterizedModel
