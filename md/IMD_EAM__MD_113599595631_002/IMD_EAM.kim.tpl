
KIM_API_Version := 1.6.0

Unit_Handling    := fixed
Unit_length      := A
Unit_energy      := eV
Unit_charge      := e
Unit_temperature := K
Unit_time        := ps

PARTICLE_SPECIES:
# Symbol/name               Type                    code
SPECIES_001_NAME_STR        spec                    0
SPECIES_002_NAME_STR        spec                    1
SPECIES_003_NAME_STR        spec                    2
SPECIES_004_NAME_STR        spec                    3
SPECIES_005_NAME_STR        spec                    4
SPECIES_006_NAME_STR        spec                    5
SPECIES_007_NAME_STR        spec                    6
SPECIES_008_NAME_STR        spec                    7
SPECIES_009_NAME_STR        spec                    8

CONVENTIONS:
# Name                      Type
ZeroBasedLists              flag
Neigh_LocaAccess 	    flag
Neigh_IterAccess            flag
NEIGH_RVEC_H                flag 	optional
NEIGH_RVEC_F                flag 	optional
NEIGH_PURE_H                flag 	optional
NEIGH_PURE_F                flag 	optional
CLUSTER                     flag 	optional

MODEL_INPUT:
# Name                      Type         Unit                Shape              Requirements
get_neigh                   method       none                []                 optional
neighObject                 pointer      none                []                 optional
numberContributingParticles integer      none 		     []
numberOfParticles           integer      none                []
numberOfSpecies         integer      none                []
particleSpecies               integer      none                [numberOfParticles]
coordinates                 double       length              [numberOfParticles,3]
process_dEdr 		    method       none 		     [] 		optional

MODEL_OUTPUT:
# Name                      Type         Unit                Shape              Requirements
compute                     method       none                []
destroy 		    method       none 	 	     [] 		optional
reinit 			    method       none                [] 		optional
cutoff                      double       length              []
energy                      double       energy              [] 		optional
forces 			    double       force 		     [numberOfParticles,3] optional
particleEnergy 		    double 	 energy 	     [numberOfParticles] optional

MODEL_PARAMETERS:
# Name                      Type         Unit                Shape              Requirements
