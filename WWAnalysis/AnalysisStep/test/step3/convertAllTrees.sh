#!/usr/bin/env bash

usage()
{
if [ $# -gt 0 ] ; then a="Error: $*"; fi
cat << EOF
$a
usage: $0 options fileToConvert [inputDir=.] [outputDir.]

mwl Wrapper for CRAB

OPTIONS:
   -h      Show this message
   -v      Verbose

EOF
}

set -- `getopt "hv" "$@"`

VERBOSE=
while [ "$1" != "--" ] ; do
  case "$1" in
    -h) usage; exit 1;;
    -v) VERBOSE=1;;
  esac
  shift
done
shift


if [ $# -lt 1 ] ; then
    usage "I would appreciate something to do ... ask tar what they think."
    exit 3
fi

inputFile=$1
inDir="."
outDir="."
shift

[ $# -gt 0 ] && { inDir=$1; shift; }
[ $# -gt 0 ] && { outDir=$1; shift; }

root -b -l -q convertToLatinos.C\(\"$inDir\",\"$outDir\",\"$inputFile\"\)
hadd -f $outDir/${inputFile/tree_/latino_} $outDir/????_$inputFile


