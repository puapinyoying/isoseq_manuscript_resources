#!/bin/sh

# Submit script for indexing genome files
#SBATCH -o sortCount_job_%j.out
#SBATCH -e sortCount_job_%j.err

# Specify debug partition
#SBATCH -p defq

# Specify number of nodes 
#SBATCH -N 1

# Specify number of cores (each node has 16 cores)
#SBATCH -n 15

# one day timelimit:
#SBATCH --time 1-00:00:00

usage() {
    cat <<EOF
Usage: $0 [options] [--] [file...]

Arguments:

  -h, --help
    Display this usage message and exit.

  -p <val>, --prefix <val>, --prefix=<val>
    Prefix name to call collapsed sorted SAM and htseq-count files

  -s <val>, --sam <val>, --sam=<val>
    Full path to unsorted collapsed sam file output from STARlong second alignment
    from step s8_secondMap (reused s6_starMap.sh to do this step)

  -g <val>, --gtf <val>, --gtf=<val>
    Reference annotation (.gtf) file. e.g. gencode.vM10.primary_assembly.annotation.gtf

  -o <val>, --outdir <val>, --outdir=<val>
    Full path to output directory; where sorted collapsed sam/bams & counts files will be
    No trailing slash in path!

  -t <val>, --threads <val>, --threads=<val>
    Number of threads/cores to use for samtools view & sort steps

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

    while [ "$#" -gt 0 ]; do
        arg=$1
        case $1 in
            # convert "--opt=the value" to --opt "the value".
            # the quotes around the equals sign is to work around a
            # bug in emacs' syntax parsing
            --*'='*) shift; set -- "${arg%%=*}" "${arg#*=}" "$@"; continue;;
            -p|--prefix) shift; PREFIX=$1;;
            -s|--sam) shift; SAM=$1;;
            -g|--gtf) shift; GTF=$1;;
            -o|--outdir) shift; OUTDIR=$1;;
            -t|--threads) shift; THREADS=$1;;
            -h|--help) usage; exit 0;;
            --) shift; break;;
            -*) usage_fatal "unknown option: '$1'";;
            *) break;; # reached the list of file names
        esac
        shift || usage_fatal "option '${arg}' requires a value"
    done

    # Load the SMRT_ANALYSIS module & python ver that has htseq-count
    module load samtools/1.2
    module load python/2.7.6

    cd ${OUTDIR}

    echo "Started samtools conv to bam, sort bam, index, conv to sam $(date)" > ${PREFIX}_sortCount.elog
    samtools view -@ ${THREADS} -Sb ${SAM} | samtools sort -@ ${THREADS} - ${PREFIX}_collapsed_sorted && samtools index ${PREFIX}_collapsed_sorted.bam
    samtools view -@ ${THREADS} -h ${PREFIX}_collapsed_sorted.bam > ${PREFIX}_collapsed_sorted.sam
    echo "Completed samtools steps $(date)" >> ${PREFIX}_sortCount.elog

    echo "Start htseq-count $(date)" >> ${PREFIX}_sortCount.elog
    htseq-count -i gene_name ${PREFIX}_collapsed_sorted.sam ${GTF} > ${PREFIX}_collapsed_isoform_count.txt
    echo "Completed htseq-count $(date)" >> ${PREFIX}_sortCount.elog
fi
