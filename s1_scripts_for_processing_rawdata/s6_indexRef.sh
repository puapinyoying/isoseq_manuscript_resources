#!/bin/sh

# Submit script for indexing genome files
#SBATCH -o indexGenome_job_%j.out
#SBATCH -e indexGenome_job_%j.err

# Specify debug partition
#SBATCH -p defq

# Specify number of nodes 
#SBATCH -N 1

# Specify number of cores (each node has 16 cores)
#SBATCH -n 15

# one hour timelimit:
#SBATCH --time 06:00:00

usage() {
    cat <<EOF
Usage: $0 [options] [--] [file...]

Arguments:

  -h, --help
    Display this usage message and exit.

  -d <val>, --genomeDir <val>, --genomeDir=<val>
    Genome directory, where to output indexed files

  -f <val>, --fasta <val>, --fasta=<val>
    Reference sequence file ('.fasta') to index e.g. gencode_mouse.fasta

  -g <val>, --gtf <val>, --gtf=<val>
    Reference annotation file ('.gtf') to index e.g. gencode_mouse.gtf

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
            -d|--genomeDir) shift; genomeDir=$1;;
            -f|--fasta) shift; fasta=$1;;
            -g|--gtf) shift; gtf=$1;;
            -t|--threads) shift; threads=$1;;
            -h|--help) usage; exit 0;;
            --) shift; break;;
            -*) usage_fatal "unknown option: '$1'";;
            *) break;; # reached the list of file names
        esac
        shift || usage_fatal "option '${arg}' requires a value"
    done

    # Load the STAR module
    module load star/2.5

    STAR --runThreadN $threads --runMode genomeGenerate --genomeDir $genomeDir \
    --genomeFastaFiles $fasta --sjdbGTFfile $gtf --sjdbOverhang 100
fi