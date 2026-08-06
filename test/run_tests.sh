#!/usr/bin/env bash
#
# Regression tests for insize.
#
# Runs the binary over test.sam (converted to BAM when samtools is available)
# and compares stdout against the files in expected/.
#
# Usage:   ./run_tests.sh            # from anywhere
#          INSIZE=/path/to/insize ./run_tests.sh
#

set -u

TESTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
INSIZE=${INSIZE:-${TESTDIR}/../src/insize}
SAMTOOLS=${SAMTOOLS:-}

if [ ! -x "${INSIZE}" ]; then
    echo "ERROR: insize binary not found at ${INSIZE}"
    echo "       build it first with:  cd src/ && make"
    exit 1
fi

# find a samtools to build the BAM with: the one we compile against, then $PATH
if [ -z "${SAMTOOLS}" ]; then
    if [ -x "${TESTDIR}/../lib/samtools/samtools" ]; then
        SAMTOOLS="${TESTDIR}/../lib/samtools/samtools"
    elif command -v samtools >/dev/null 2>&1; then
        SAMTOOLS=$( command -v samtools )
    fi
fi

WORKDIR=$( mktemp -d )
trap 'rm -rf "${WORKDIR}"' EXIT

SAM="${TESTDIR}/test.sam"
BAM="${WORKDIR}/test.bam"

if [ -n "${SAMTOOLS}" ]; then
    "${SAMTOOLS}" view -b -o "${BAM}" "${SAM}" || { echo "ERROR: could not convert ${SAM} to BAM"; exit 1; }
    INPUT="${BAM}"
    LABEL="<bam>"
    echo "input: BAM (converted with ${SAMTOOLS})"
else
    # htslib autodetects the format, so the SAM works as a fallback
    INPUT="${SAM}"
    LABEL="<sam>"
    echo "input: SAM (samtools not found, skipping the BAM round-trip)"
fi
echo

NPASS=0
NFAIL=0

# run_case <name> <expected file> <insize options...>
run_case(){
    local name=$1; shift
    local expected="${TESTDIR}/expected/$1"; shift

    local got="${WORKDIR}/${name}.out"
    "${INSIZE}" "$@" "${INPUT}" > "${got}" 2>"${WORKDIR}/${name}.err"
    local rc=$?

    if [ ${rc} -ne 0 ]; then
        echo "FAIL ${name}  (insize $* ${LABEL}): exit code ${rc}"
        sed 's/^/       /' "${WORKDIR}/${name}.err"
        NFAIL=$(( NFAIL + 1 ))
        return
    fi

    if diff -u "${expected}" "${got}" > "${WORKDIR}/${name}.diff"; then
        echo "PASS ${name}  (insize $* ${LABEL})"
        NPASS=$(( NPASS + 1 ))
    else
        echo "FAIL ${name}  (insize $* ${LABEL})"
        sed 's/^/       /' "${WORKDIR}/${name}.diff"
        NFAIL=$(( NFAIL + 1 ))
    fi
}

# all fragments: 3 paired (100,150,200), the short/long pairs (30,500)
# and the 3 single-end reads (75,60,45)
run_case default          default.txt
# -m drops the unmapped single-end read (60)
run_case mapped           mapped.txt            -m
# -p drops the pair that is not flagged as properly paired (200)
run_case properpair       properpair.txt        -p
run_case mapped_properpair mapped_properpair.txt -m -p
# -l/-L keep only fragments in [50,300]
run_case range_50_300     range_50_300.txt      -l 50 -L 300
# -v prints the histogram over [-l,-L]; only the 100 bp pair falls in this window
run_case verbose_90_110   verbose_90_110.txt    -v -l 90 -L 110
# nothing in range -> NA
run_case verbose_empty    verbose_empty.txt     -v -l 1000 -L 1100

# the SAM itself must give the same answer as the BAM
if [ "${INPUT}" != "${SAM}" ]; then
    got="${WORKDIR}/sam_input.out"
    "${INSIZE}" "${SAM}" > "${got}" 2>/dev/null
    if diff -u "${TESTDIR}/expected/default.txt" "${got}" > "${WORKDIR}/sam_input.diff"; then
        echo "PASS sam_input  (insize <sam>)"
        NPASS=$(( NPASS + 1 ))
    else
        echo "FAIL sam_input  (insize <sam>)"
        sed 's/^/       /' "${WORKDIR}/sam_input.diff"
        NFAIL=$(( NFAIL + 1 ))
    fi
fi

# no arguments must print the usage and fail
if "${INSIZE}" >/dev/null 2>&1; then
    echo "FAIL usage  (no arguments should exit non-zero)"
    NFAIL=$(( NFAIL + 1 ))
else
    echo "PASS usage  (no arguments exits non-zero)"
    NPASS=$(( NPASS + 1 ))
fi

# a missing input file must fail rather than report an empty result
if "${INSIZE}" "${WORKDIR}/does_not_exist.bam" >/dev/null 2>&1; then
    echo "FAIL missing_file  (missing input should exit non-zero)"
    NFAIL=$(( NFAIL + 1 ))
else
    echo "PASS missing_file  (missing input exits non-zero)"
    NPASS=$(( NPASS + 1 ))
fi

echo
echo "${NPASS} passed, ${NFAIL} failed"
[ ${NFAIL} -eq 0 ]
