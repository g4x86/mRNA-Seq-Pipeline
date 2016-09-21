#! /bin/bash

# Directories
PROJ_DIR="$HOME/LINCS/Alignment"
PROG_DIR="${PROJ_DIR}/Programs"
DATA_DIR="${PROJ_DIR}/Repo/LINCS.Dataset.20150409"
REF_DIR="${PROJ_DIR}/References/Broad_UMI"
# Reference files
SPECIES_DIR="Human_RefSeq"
REF_SEQ_FILE="refMrna_ERCC_polyAstrip.hg19.fa"
SYM2REF_FILE="refGene.hg19.sym2ref.dat"
ERCC_SEQ_FILE="ERCC92.fa"
BARCODE_FILE="barcodes_trugrade_96_set4.dat"
# Data files
SEQ_DIR="${DATA_DIR}/Seqs"
ALIGN_DIR="${DATA_DIR}/Aligns"
COUNT_DIR="${DATA_DIR}/Counts"
# Parallel parameter
THREAD_NUMBER=8

# Start converting RNA sequencing reads to gene read counts.

#python "${PROG_DIR}/run_DGE_analysis.py" "${DATA_DIR}/SampleMap.txt" "Human" "Trugrade_96_set4" "${DATA_DIR}/Aligns" "${DATA_DIR}/Counts"

# Main sample ID
SAMPLE_ID="Undetermined"

# Processing 6 lanes of FASTQ data
LANES=6
let "i = 1"
while [ $i -le $LANES ]; do
	SUBSAMPLE_ID="Lane$i"
	SEQ_FILE_R1="${SAMPLE_ID}_${SUBSAMPLE_ID}_R1.fastq.gz"
	SEQ_FILE_R2="${SAMPLE_ID}_${SUBSAMPLE_ID}_R2.fastq.gz"
	# Split input paired FASTQ files to multiple intermediate FASTQ files.
	python "${PROG_DIR}/split_and_align2.py" "${SAMPLE_ID}" "${SUBSAMPLE_ID}" "${SEQ_DIR}/${SEQ_FILE_R1}" "${SEQ_DIR}/${SEQ_FILE_R2}" "${ALIGN_DIR}/" "${REF_DIR}/${SPECIES_DIR}/${REF_SEQ_FILE}" write
	# Align individual intermediate FASTQ files.
	OLD_IFS=$IFS
	IFS=$'\n'
	for SEQ_FILE in $(find "${ALIGN_DIR}" -name "${SAMPLE_ID}\.${SUBSAMPLE_ID}*\.fastq" -printf "%p$IFS"); do
		SAM_FILE="${SEQ_FILE}.sam"
		bwa aln -l 24 -t "${THREAD_NUMBER}" "${REF_DIR}/${SPECIES_DIR}/${REF_SEQ_FILE}" "${SEQ_FILE}" | bwa samse -n 20 "${REF_DIR}/${SPECIES_DIR}/${REF_SEQ_FILE}" - "${SEQ_FILE}" > "${SAM_FILE}"
	done
	IFS=$OLD_IFS
	let "i = $i + 1"
done

# Merge all lanes of mapped reads and count total reads for each genes.
python "${PROG_DIR}/merge_and_count.py" "${SAMPLE_ID}" "${REF_DIR}/${SPECIES_DIR}/${SYM2REF_FILE}" "${REF_DIR}/${ERCC_SEQ_FILE}" "${REF_DIR}/${BARCODE_FILE}" "${ALIGN_DIR}/" "${COUNT_DIR}/" False

exit 0
