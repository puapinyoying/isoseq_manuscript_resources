#!/bin/sh

# Submit script for indexing genome files
#SBATCH -o mapFastq_job_%j.out
#SBATCH -e mapFastq_job_%j.err

# Specify debug partition
#SBATCH -p defq

# Specify number of nodes 
#SBATCH -N 1

# Specify number of cores (each node has 16 cores)
#SBATCH -n 15

# one hour timelimit:
#SBATCH --time 14-00:00:00

usage() {
    cat <<EOF
Usage: $0 [options] [--] [file...]

Arguments:

  -h, --help
    Display this usage message and exit.

  -d <val>, --genomeDir <val>, --genomeDir=<val>
    Genome directory, where to output indexed files

  -f <val>, --fastq <val>, --fastq=<val>
    Full path to mereged, quivered, high-quality '.fastq' file from s5_mergeFastq step

  -t <val>, --threads <val>, --threads=<val>
    Number of threads to use

EOF
}

# handy logging and error handling functions
log() { printf '%s\n' "$*"; }
error() { log "ERROR: $*" >&2; }
fatal() { error "$*"; exit 1; }
usage_fatal() { error "$*"; usage >&2; exit 1; }

# Check if there are any inputs
if [ $# -eq 0 ]; then
    usage
    printf "** Error: No inputs.**\n\n"

else
    # parse options
    threads=15  # default = 16 cores per node, use only 15

    while [ "$#" -gt 0 ]; do
        arg=$1
        case $1 in
            # convert "--opt=the value" to --opt "the value".
            # the quotes around the equals sign is to work around a
            # bug in emacs' syntax parsing
            --*'='*) shift; set -- "${arg%%=*}" "${arg#*=}" "$@"; continue;;
            -d|--genomeDir) shift; GENOME_DIR=$1;;
            -f|--fastq) shift; FASTQ=$1;;
            -t|--threads) shift; THREADS=$1;;
            -h|--help) usage; exit 0;;
            --) shift; break;;
            -*) usage_fatal "unknown option: '$1'";;
            *) break;; # reached the list of file names
        esac
        shift || usage_fatal "option '${arg}' requires a value"
    done

    # Load the STAR module
    module load star/2.5

    STARlong --runThreadN $THREADS \
    --genomeDir $GENOME_DIR \
    --runMode alignReads \
    --outSAMattributes NH HI NM MD \
    --readNameSeparator space \
    --outFilterMultimapScoreRange 1 \
    --outFilterMismatchNmax 2000 \
    --scoreGapNoncan -20 \
    --scoreGapGCAG -4 \
    --scoreGapATAC -8 \
    --scoreDelOpen -1 \
    --scoreDelBase -1 \
    --scoreInsOpen -1 \
    --scoreInsBase -1 \
    --alignEndsType Local \
    --seedSearchStartLmax 50 \
    --seedPerReadNmax 100000 \
    --seedPerWindowNmax 1000 \
    --alignTranscriptsPerReadNmax 100000 \
    --alignTranscriptsPerWindowNmax 10000 \
    --readFilesIn $FASTQ
fi
