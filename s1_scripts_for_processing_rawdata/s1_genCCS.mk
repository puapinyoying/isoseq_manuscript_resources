
######################################################################
### Creating submit scripts for generating CCS reads from raw data ###
######################################################################

# Generate submit scripts for generating circular consensus data from rawdata
# Many path and file variables are defined in the copy_rawdata.mk file

# Use for python script --dataRoot
LUSTRE_DATA_ROOT := $(TARGET_ROOT)
LUSTRE_RESULTS_ROOT := $(addprefix $(LUSTRE_ROOT)/, results)

STEP1 = s1_ccs
CELL_DESC = indv_cell

SIZE1 = 3-6kb
SIZE2 = 5-10kb
SIZE3 = 7-15kb

TISSUE1 = soleus
TISSUE2 = edl
TISSUE3 = cardiac
TISSUE4 = pooled

# Generate submit scripts for generating ccs reads by feeding directory lists
# to an external python script. The script takes 7 arguments: '--runFolder',
# '--cellList', '--dataRoot', '--resultsRoot', '--tissue',
# '--shortFolderName' and '--size'

# Python script to call
CCS_SCRIPT_SRC := s1_genCCS.py
CCS_SCRIPT_EXE := python $(CCS_SCRIPT_SRC)

# 3kb size fraction
.PHONY : P3A_genCCS P3B_genCCS P3C_genCCS

P3A_genCCS :  $(CCS_SCRIPT_SRC)
	$(CCS_SCRIPT_EXE) --size $(SIZE1) --runFolder $(M2T1) --cellList $(P3A_CELL) \
	--dataRoot $(LUSTRE_DATA_ROOT) --resultsRoot $(LUSTRE_RESULTS_ROOT) \
	--fofnRoot $(FOFN_ROOT) --scriptsRoot $(SCRIPTS_ROOT) --logRoot $(LOG_ROOT) \
	--tissue $(TISSUE4) --shortFolderName M2T1 

P3B_genCCS : $(CCS_SCRIPT_SRC)
	$(CCS_SCRIPT_EXE) --size $(SIZE1) --runFolder $(M2T2) --cellList $(P3B_CELL) \
	--dataRoot $(LUSTRE_DATA_ROOT) --resultsRoot $(LUSTRE_RESULTS_ROOT) \
	--fofnRoot $(FOFN_ROOT) --scriptsRoot $(SCRIPTS_ROOT) --logRoot $(LOG_ROOT) \
	--tissue $(TISSUE4) --shortFolderName M2T2

P3C_genCCS : $(CCS_SCRIPT_SRC)
	$(CCS_SCRIPT_EXE) --size $(SIZE1) --runFolder $(M2F3) --cellList $(P3C_CELL) \
	--dataRoot $(LUSTRE_DATA_ROOT) --resultsRoot $(LUSTRE_RESULTS_ROOT) \
	--fofnRoot $(FOFN_ROOT) --scriptsRoot $(SCRIPTS_ROOT) --logRoot $(LOG_ROOT) \
	--tissue $(TISSUE4) --shortFolderName M2F3

# 5-10kb size fractions
.PHONY :  P5A_genCCS P5B_genCCS P5C_genCCS
P5A_genCCS : $(CCS_SCRIPT_SRC)
	$(CCS_SCRIPT_EXE) --size $(SIZE2) --runFolder $(M2T1) --cellList $(P5A_CELL) \
	--dataRoot $(LUSTRE_DATA_ROOT) --resultsRoot $(LUSTRE_RESULTS_ROOT) \
	--fofnRoot $(FOFN_ROOT) --scriptsRoot $(SCRIPTS_ROOT) --logRoot $(LOG_ROOT) \
	--tissue $(TISSUE4) --shortFolderName M2T1

P5B_genCCS : $(CCS_SCRIPT_SRC)
	$(CCS_SCRIPT_EXE) --size $(SIZE2) --runFolder $(M2T2) --cellList $(P5B_CELL) \
	--dataRoot $(LUSTRE_DATA_ROOT) --resultsRoot $(LUSTRE_RESULTS_ROOT) \
	--fofnRoot $(FOFN_ROOT) --scriptsRoot $(SCRIPTS_ROOT) --logRoot $(LOG_ROOT) \
	--tissue $(TISSUE4) --shortFolderName M2T2

P5C_genCCS : $(CCS_SCRIPT_SRC)
	$(CCS_SCRIPT_EXE) --size $(SIZE2) --runFolder $(M2F5) --cellList $(P5C_CELL) \
	--dataRoot $(LUSTRE_DATA_ROOT) --resultsRoot $(LUSTRE_RESULTS_ROOT) \
	--fofnRoot $(FOFN_ROOT) --scriptsRoot $(SCRIPTS_ROOT) --logRoot $(LOG_ROOT) \
	--tissue $(TISSUE4) --shortFolderName M2F5

# 7-15kb size fractions
.PHONY : S7_genCCS E7A_genCCS E7B_genCCS E7C_genCCS C7_genCCS
S7_genCCS : $(CCS_SCRIPT_SRC)
	$(CCS_SCRIPT_EXE) --size $(SIZE3) --runFolder $(SE7_16CELL) --cellList $(S7_CELL) \
	--dataRoot $(LUSTRE_DATA_ROOT) --resultsRoot $(LUSTRE_RESULTS_ROOT) \
	--fofnRoot $(FOFN_ROOT) --scriptsRoot $(SCRIPTS_ROOT) --logRoot $(LOG_ROOT) \
	--tissue $(TISSUE1) --shortFolderName SE7

E7A_genCCS : $(CCS_SCRIPT_SRC)
	$(CCS_SCRIPT_EXE) --size $(SIZE3) --runFolder $(E7_8CELL) --cellList $(E7A_CELL) \
	--dataRoot $(LUSTRE_DATA_ROOT) --resultsRoot $(LUSTRE_RESULTS_ROOT) \
	--fofnRoot $(FOFN_ROOT) --scriptsRoot $(SCRIPTS_ROOT) --logRoot $(LOG_ROOT) \
	--tissue $(TISSUE2) --shortFolderName E7_8CELL

E7B_genCCS : $(CCS_SCRIPT_SRC)
	$(CCS_SCRIPT_EXE) --size $(SIZE3) --runFolder $(CE7_16CELL) --cellList $(E7B_CELL) \
	--dataRoot $(LUSTRE_DATA_ROOT) --resultsRoot $(LUSTRE_RESULTS_ROOT) \
	--fofnRoot $(FOFN_ROOT) --scriptsRoot $(SCRIPTS_ROOT) --logRoot $(LOG_ROOT) \
	--tissue $(TISSUE2) --shortFolderName CE7_16CELL

E7C_genCCS : $(CCS_SCRIPT_SRC)
	$(CCS_SCRIPT_EXE) --size $(SIZE3) --runFolder $(SE7_16CELL) --cellList $(E7C_CELL) \
	--dataRoot $(LUSTRE_DATA_ROOT) --resultsRoot $(LUSTRE_RESULTS_ROOT) \
	--fofnRoot $(FOFN_ROOT) --scriptsRoot $(SCRIPTS_ROOT) --logRoot $(LOG_ROOT) \
	--tissue $(TISSUE2) --shortFolderName SE7_16CELL

C7_genCCS : $(CCS_SCRIPT_SRC)
	$(CCS_SCRIPT_EXE) --size $(SIZE3) --runFolder $(CE7_16CELL) --cellList $(C7_CELL) \
	--dataRoot $(LUSTRE_DATA_ROOT) --resultsRoot $(LUSTRE_RESULTS_ROOT) \
	--fofnRoot $(FOFN_ROOT) --scriptsRoot $(SCRIPTS_ROOT) --logRoot $(LOG_ROOT) \
	--tissue $(TISSUE3) --shortFolderName CE7_16CELL
