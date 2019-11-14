# Create variables for source and target raw data files to copy from Source drive to lustre

# Change depending on source and target drives
SOURCE_ROOT = $(GROUPS_ROOT)
TARGET_ROOT = $(LUSTRE_ROOT)

# data directory from source drive (vars taken from main makefile)
SOURCE_DATA := $(addsuffix /data, $(SOURCE_ROOT))

# Target directory to copy source data to
TARGET_DATA := $(addsuffix /data, $(TARGET_ROOT))

# RUN_DIR (sequencing run directories)
## These tissue samples are multiplexed M2 = mouse 2
M2T1 := 2015-11-02_M2_Loading_Titration1
M2T2 := 2015-11-06_M2_Loading_Titration2
M2F5 := 2015-11-17_M2_5-10kb_FullRun
M2F3 := 2015-11-20_M2_3-6kb_FullRun
E7_8CELL := 2016-06-29_EDL_7-15kb_8cell
CE7_16CELL := 2016-07-07_Cardiac_EDL_7-15kb_16cell
SE7_16CELL := 2016-07-13_Soleus_EDL_7-15kb_16cell

######################
###  Pooled 3-6kb  ###
######################

# P3A = Pooled, 3-6kb, RunFolder A (2015-11-02_M2_Loading_Titration1)
# B = 2015-11-06_M2_Loading_Titration2
# C = 2015-11-20_M2_3-6kb_FullRun

# Subfolders containing cells with 3-6kb size fractions
P3A_CELL := A01_1 B01_1 C01_1 D01_1
P3B_CELL := A01_1 B01_1 C01_1 D01_1
# B01_1 cell failed; excluded
P3C_CELL := A01_1 C01_1 D01_1 E01_1 F01_1 G01_1 H01_1 A02_1 B02_1 C02_1 D02_1 E02_1 F02_1 G02_1 H02_1

P3A_SOURCE_CELL_PATH := $(addprefix $(addprefix $(SOURCE_DATA)/, $(M2T1))/, $(P3A_CELL))
P3B_SOURCE_CELL_PATH := $(addprefix $(addprefix $(SOURCE_DATA)/, $(M2T2))/, $(P3B_CELL))
P3C_SOURCE_CELL_PATH := $(addprefix $(addprefix $(SOURCE_DATA)/, $(M2F3))/, $(P3C_CELL))

#######################
###  Pooled 5-10kb  ###
#######################

P5A_CELL := E01_1 F01_1 G01_1 H01_1
P5B_CELL := E01_1 F01_1 G01_1 H01_1
P5C_CELL := A01_1 B01_1 C01_1 D01_1 E01_1 F01_1 G01_1 H01_1 A02_1 B02_1 C02_1 D02_1 E02_1 F02_1 G02_1 H02_1

P5A_SOURCE_CELL_PATH := $(addprefix $(addprefix $(SOURCE_DATA)/, $(M2T1))/, $(P5A_CELL))
P5B_SOURCE_CELL_PATH := $(addprefix $(addprefix $(SOURCE_DATA)/, $(M2T2))/, $(P5B_CELL))
P5C_SOURCE_CELL_PATH := $(addprefix $(addprefix $(SOURCE_DATA)/, $(M2F5))/, $(P5C_CELL))

################
###  7-15kb  ###
################

# Slightly different here because samples were not multiplexed and just ran on
# separate cells.

## SOLEUS
S7_CELL := A01_1 B01_1 C01_1 D01_1 E01_1 F01_1 G01_1 H01_1 A02_1 B02_1 C02_1 D02_1
S7_SOURCE_CELL_PATH := $(addprefix $(addprefix $(SOURCE_DATA)/, $(SE7_16CELL))/, $(S7_CELL))

## EDL CELL_DIR for E7_8CELL, CE7_16CELL and SE7_16CELL
E7A_CELL := A01_1 B01_1 C01_1 D01_1 E01_1 F01_1 G01_1 H01_1
E7B_CELL := E02_1 F02_1 G02_1 H02_1
E7C_CELL := E02_1 F02_1 G02_1 H02_1

E7A_SOURCE_CELL_PATH := $(addprefix $(addprefix $(SOURCE_DATA)/, $(E7_8CELL))/, $(E7A_CELL))
E7B_SOURCE_CELL_PATH := $(addprefix $(addprefix $(SOURCE_DATA)/, $(CE7_16CELL))/, $(E7B_CELL))
E7C_SOURCE_CELL_PATH := $(addprefix $(addprefix $(SOURCE_DATA)/, $(SE7_16CELL))/, $(E7C_CELL))

## CARDIAC  (cell F01_1 failed; excluded)
C7_CELL := A01_1 B01_1 C01_1 D01_1 E01_1 G01_1 H01_1 A02_1 B02_1 C02_1 D02_1
C7_SOURCE_CELL_PATH := $(addprefix $(addprefix $(SOURCE_DATA)/, $(CE7_16CELL))/, $(C7_CELL))
C7_TARGET_PATH := $(addprefix $(TARGET_DATA)/, $(CE7_16CELL))

####################
### TARGET PATHS ###
####################

M2T1_TARGET_PATH := $(addprefix $(TARGET_DATA)/, $(M2T1))
M2T2_TARGET_PATH := $(addprefix $(TARGET_DATA)/, $(M2T2))
M2F3_TARGET_PATH := $(addprefix $(TARGET_DATA)/, $(M2F3))
M2F5_TARGET_PATH := $(addprefix $(TARGET_DATA)/, $(M2F5))
E7_TARGET_PATH := $(addprefix $(TARGET_DATA)/, $(E7_8CELL))
CE7_TARGET_PATH := $(addprefix $(TARGET_DATA)/, $(CE7_16CELL))
SE7_TARGET_PATH := $(addprefix $(TARGET_DATA)/, $(SE7_16CELL))

#################
### FUNCTIONS ###
#################

# There are 3 run folders containing 3-6kb data 
.PHONY : copy_P3A copy_P3B copy_P3C copy_P5A copy_P5B copy_P5C
copy_P3A : $(P3A_SOURCE_CELL_PATH) $(M2T1_TARGET_PATH) 
	rsync -rhvupE $^

copy_P3B : $(P3B_SOURCE_CELL_PATH) $(M2T2_TARGET_PATH) 
	rsync -rhvupE $^

copy_P3C : $(P3C_SOURCE_CELL_PATH) $(M2F3_TARGET_PATH)
	rsync -rhvupE $^

# 3 run folders containing 5-10kb data
.PHONY : copy_P5A copy_P5B copy_P5C
copy_P5A : $(P5A_SOURCE_CELL_PATH) $(M2T1_TARGET_PATH)
	rsync -rhvupE $^

copy_P5B : $(P5B_SOURCE_CELL_PATH) $(M2T2_TARGET_PATH)
	rsync -rhvupE $^

copy_P5C : $(P5C_SOURCE_CELL_PATH) $(M2F5_TARGET_PATH)
	rsync -rhvupE $^

# 5 different run folders for 7-15kb data
.PHONY : copy_S7 copy_E7A copy_E7B copy_E7C copy_C7
copy_S7 : $(S7_SOURCE_CELL_PATH) $(SE7_TARGET_PATH)
	rsync -rhvupE $^

copy_E7A : $(E7A_SOURCE_CELL_PATH) $(E7_TARGET_PATH)
	rsync -rhvupE $^

copy_E7B : $(E7B_SOURCE_CELL_PATH) $(CE7_TARGET_PATH)
	rsync -rhvupE $^

copy_E7C : $(E7C_SOURCE_CELL_PATH) $(SE7_TARGET_PATH)
	rsync -rhvupE $^

copy_C7 : $(C7_SOURCE_CELL_PATH) $(CE7_TARGET_PATH)
	rsync -rhvupE $^

# Create directories for copying the size fractions to in lustre drive
$(M2T1_TARGET_PATH) :
	mkdir -p $@

$(M2T2_TARGET_PATH) :
	mkdir -p $@

$(M2F3_TARGET_PATH) :
	mkdir -p $@

$(M2F5_TARGET_PATH) :
	mkdir -p $@

$(E7_TARGET_PATH) :
	mkdir -p $@

$(SE7_TARGET_PATH) :
	mkdir -p $@

$(CE7_TARGET_PATH) :
	mkdir -p $@


############
### TEST ###
############

TE7_CELL := A01_1 
TE7_SOURCE_CELL_PATH := $(addprefix $(addprefix $(SOURCE_DATA)/, $(E7_8CELL))/, $(TE7_CELL))

copy_test : $(TE7_SOURCE_CELL_PATH) $(E7_TARGET_PATH)
	rsync -nrhvupE $^

.PHONY : test_vars
test_vars :
	@echo "TE7_CELL : $(TE7_CELL)" > test_results.txt
	@echo "TE7_SOURCE_CELL_PATH : $(TE7_SOURCE_CELL_PATH)" >> test_results.txt
	@echo "TE7_TARGET_PATH : $(TE7_TARGET_PATH)" >> test_results.txt