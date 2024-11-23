# Note: This took about 10 minutes to run with the main GENIE files on a personal
#.  computer.

#!/usr/bin/env bash

# input and output maf files
IMAF="/data/genomic/main_genie/mut_ready_to_annotate.txt"
OMAF="/data/genomic/main_genie/mut_onco.txt"

# The input CNA file is shaped and filtered from the raw input.
ICNA="/data/genomic/main_genie/cna_ready_to_annotate.txt"
OCNA="/data/genomic/main_genie/cna_onco.txt"

IFUS="/data/genomic/main_genie/fus_ready_to_annotate.txt"
OFUS="/data/genomic/main_genie/fus_onco.txt"

# Three environment variables needed:
# PROJ_ROOT - the location of the R project file.
# ONCO_ANNO_LOC - the location of the oncoKB annotator scripts.
# ONCOKB_KEY - the API key for oncoKB.

# For all of these you can set with commands like "export ONCOKB_KEY=''"
# check that it exists with echo ${ONCOKB_KEY}

# It worked best for me to use the full path e.g. "/Users/.../oncokb-annotator"

# OncoKB annotator is a program available from a github repo:
# https://github.com/oncokb/oncokb-annotator
# Clone this to a directory of your choosing and reference it with ONCO_ANNO_LOC above.
# The location you refer to should contain MafAnnotator.py (for example).
python ${ONCO_ANNO_LOC}/MafAnnotator.py -i "${PROJ_ROOT}${IMAF}" -o "${PROJ_ROOT}${OMAF}" -b ${ONCOKB_KEY} -r GRCh37 
# The -z flag for CNA gives us the Gain/Loss data.  This may have limited value,
#   and it's addressed in the help file for the script.
python ${ONCO_ANNO_LOC}/CnaAnnotator.py -i "${PROJ_ROOT}${ICNA}" -o "${PROJ_ROOT}${OCNA}" -b ${ONCOKB_KEY} -z
python ${ONCO_ANNO_LOC}/FusionAnnotator.py -i "${PROJ_ROOT}${IFUS}" -o "${PROJ_ROOT}${OFUS}" -b ${ONCOKB_KEY}  