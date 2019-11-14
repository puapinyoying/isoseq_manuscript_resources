#!/bin/sh

# Submit script for indexing genome files
#SBATCH -o collapse_job_%j.out
#SBATCH -e collapse_job_%j.err

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
    Prefix name to call sorted SAM and collapsed isoform files

  -f <val>, --fastq <val>, --fastq=<val>
    Full path to mereged, quivered, high-quality '.fastq' file from s5_mergeFastq step

  -s <val>, --sam <val>, --sam=<val>
    Full path to unsorted sam file output from STARlong alignment / s6_starMap step

  -o <val>, --outdir <val>, --outdir=<val>
    Full path to output directory; where collapsed sam files & logs will be
    No trailing slash in path!

  -t <val>, --threads <val>, --threads=<val>
    Number of threads/cores to use for samtools view & sort steps

  -u <val>, --tofu <val>, --tofu=<val>
    Full path to virtualenv directory of installed pbtranscript-ToFU
    e.g. '/home/brianu/vtofu/VENV_TOFU'.  Should be able to find the 'bin'
    folder directly inside the provided path. No trailing slash in path!

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
            -f|--fastq) shift; FASTQ=$1;;
            -s|--sam) shift; SAM=$1;;
            -o|--outdir) shift; OUTDIR=$1;;
            -t|--threads) shift; THREADS=$1;;
            -u|--tofu) shift; TOFU=$1;;
            -h|--help) usage; exit 0;;
            --) shift; break;;
            -*) usage_fatal "unknown option: '$1'";;
            *) break;; # reached the list of file names
        esac
        shift || usage_fatal "option '${arg}' requires a value"
    done

    # House keeping: clean trailing slashes from directory paths
    #printf -v cTOFU '%s' "${TOFU}[@]%/}"
    #printf -v cOUTDIR '%s' "${OUTDIR}[@]%/}"

    # Load the SMRT_ANALYSIS module & Samtools
    module load smrt_analysis/2.3.0
    module load samtools/1.2

    # Load the enviornmental variables & Smrtshell
    source /c1/apps/smrt_analysis/2.3.0/current/etc/setup.sh
    /c1/apps/smrt_analysis/2.3.0/smrtcmds/bin/smrtshell

    # Activate ToFU - must be installed in a virtualenv first
    # See "Installing_pbtranscript-TOFU.md"
    source ${TOFU}/bin/activate

    cd ${OUTDIR}

    echo "echo Started samtools conv to bam, sort bam, index bam, conv to sam $(date)" > ${OUTDIR}/${PREFIX}_collapse.elog
    samtools view -@ ${THREADS} -Sb ${SAM} | samtools sort -@ ${THREADS} - ${PREFIX}_sorted && samtools index ${PREFIX}_sorted.bam
    samtools view -@ ${THREADS} -h ${PREFIX}_sorted.bam > ${PREFIX}_sorted.sam
    echo "Completed samtools sort $(date)" >> ${OUTDIR}/${PREFIX}_collapse.elog

    echo "Started collapse_isoforms_by_sam.py $(date)" >> ${OUTDIR}/${PREFIX}_collapse.elog
    collapse_isoforms_by_sam.py --input ${FASTQ} --fq \
    -s ${PREFIX}_sorted.sam \
    -o ${PREFIX} \
    >> ${OUTDIR}/${PREFIX}_collapse.olog 2>> ${OUTDIR}/${PREFIX}_collapse.elog
    echo "Finsihed $(date)" >> ${OUTDIR}/${PREFIX}_collapse.elog

fi

