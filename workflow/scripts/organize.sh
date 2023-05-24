#!/bin/bash

usage() { echo "Usage: $0 [-n <ont>] [-i <illumina>] [-r <reference>] [-t <table>]" 1>&2; exit 1; }

while getopts ":n:i:r:t:" o; do
    case "${o}" in
        n)
            n=${OPTARG}
            ;;
        i)
            i=${OPTARG}
            ;;
        r)
            r=${OPTARG}
            ;;
        t)
            t=${OPTARG}
            ;;
        h | *)
            usage
            exit 0
            ;;
        ?)
            usage
            exit 0
            ;;
    esac
done
shift $((OPTIND-1))

#define all base directories
script_dir=$(dirname -- "$0")
snake_base=$(readlink -f $script_dir/../..)

#reference - symlink and config
reference_base=$(readlink -f $snake_base/resources/reference)
mkdir -p $ont_base $illumina_base $reference_base
reference_path=$(readlink -f $r)
reference_name=$(basename $reference_path)
ln -sf $reference_path $reference_base/$reference_name
echo -e "reference: resources/reference/$reference_name" > config/config.yaml
echo -e "samples: samples.tsv" >> config/config.yaml

#table
table=$(readlink -f $t)
pattern_out=$(tail -n+2 $table | cut -f 3)
pattern_ont=$(tail -n+2 $table | cut -f 1)
pattern_ill=$(tail -n+2 $table | cut -f 2)

#ont - symlink
ont_path=$(readlink -f $n)
#iterate over 2 variables - nanopore
set $pattern_out
for l in $pattern_ont; do
    ont_out=$snake_base/resources/reads/"$1"/ont
    mkdir -p $ont_out
    #echo $l $ont_out
    files=$(find $ont_path -name *$l*)
    #echo $files
    for f in $files; do 
        filename=$(basename -- $f)
        ln -sf $f $ont_out/$filename
    done
    shift
done

#illumina - symlink
ill_path=$(readlink -f $i)
#iterate over 2 variables - nanopore
set $pattern_out
for l in $pattern_ill; do
    ill_out=$snake_base/resources/reads/"$1"/illumina
    mkdir -p $ill_out
    #echo $l $ill_out
    files=$(find $ill_path -name *$l*)
    #echo $files
    for f in $files; do 
        filename=$(basename -- $f)
        ln -sf $f $ill_out/$filename
    done
    shift
done

#samples
echo -e "sample\tpath" > config/samples.tsv

for l in $pattern_out; do

    echo -e $l"\t"$snake_base/resources/reads/"$l" >>config/samples.tsv

done