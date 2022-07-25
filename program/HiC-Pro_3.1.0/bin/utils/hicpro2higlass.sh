#!/bin/bash

## HiC-Pro           
## Copyleft 2017 Institut Curie
## Author(s): Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD licence

##
## First version of converter between HiCPro and higlass.
## The cooler python package should be properly installed, as well as the higlass software
##

##
## A few notes about higlass
##
## docker run will install the docker image and start it
## sudo docker run --detach  --publish 8888:80  --volume ~/hg-data:/data --volume ~/hg-tmp:/tmp --name higlass-container  gehlenborglab/higlass
## sudo docker start higlass-container
## sudo docker ps -all
##
## Once higlass is installed, you can just run it using
## sudo docker start higlass-container
## higlass will then be available at http://localhost:8888
##

###########################
## trap handler
###########################
function trap_error()
{    
    echo "Error: $1 - line $2 - exit status of last command: $?. Exit" >&2
    exit 1
}

function trap_exit()
{
    ##Since bash-4.0 $LINENO is reset to 1 when the trap is triggered
    if [ "$?" != "0" ]; then
	echo "Error: exit status detected. Exit." >&2
    fi

    if [[ ! -z ${tmp_dir} && -e ${tmp_dir} ]]; then 
	echo -e "Cleaning temporary folders ..." >&2
	/bin/rm -rf ${tmp_dir}
    fi
}

trap 'trap_error "$0" "$LINENO"' ERR
trap 'trap_exit' 0 1 2 3

set -E ## export trap to functions
set -o pipefail  ## trace ERR through pipes         

## 0 =
## 1 >
## 2 <
vercomp () {
    if [[ $1 == $2 ]]
    then
        return 0
    fi
    local IFS=.
    local i ver1=($1) ver2=($2)
    # fill empty fields in ver1 with zeros
    for ((i=${#ver1[@]}; i<${#ver2[@]}; i++))
    do
        ver1[i]=0
    done

    for ((i=0; i<${#ver1[@]}; i++))
    do
        if [[ -z ${ver2[i]} ]]
        then
            # fill empty fields in ver2 with zeros
            ver2[i]=0
        fi
        if ((10#${ver1[i]} > 10#${ver2[i]}))
        then
            echo 1
        fi
        if ((10#${ver1[i]} < 10#${ver2[i]}))
        then
            echo 2
        fi
    done
    echo 0
}

function usage {
    echo -e "usage : hicpro2higlass -i INPUT -r RESOLUTION -c CHROMSIZE [-n] [-o ODIR] [-t TEMP] [-f FLOAT] [-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "Generate Higlass input file from HiC-Pro results"
    echo "See https://github.com/hms-dbmi/higlass-website for details about Higlass"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -i|--input INPUT : allValidPairs or matrix file generated by HiC-Pro"
    echo "   -r|--res RESOLUTION : .matrix file resolution or maximum resolution to reach from the .allValidPairs input file"
    echo "   -c|--chrom CHROMSIZE : chromosome size file"
    echo "   -p|--proc NB_CPU : number of CPUs for cooler"
    echo "   [-n|--norm] : run cooler matrix balancing algorithm"
    echo "   [-f|--float] : cooler count as float for already normalized data"
    echo "   [-o|--out] : output path. Default is current path"
    echo "   [-t|--temp] TEMP : path to tmp folder. Default is current path"
    echo "   [-h|--help]: help"
    exit;
}


if [ $# -lt 1 ]
then
    usage
    exit
fi

# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--input") set -- "$@" "-i" ;;
      "--bed")   set -- "$@" "-b" ;;
      "--res")   set -- "$@" "-r" ;;
      "--chrom") set -- "$@" "-c" ;;
      "--proc") set -- "$@" "-p" ;;
      "--out") set -- "$@" "-o" ;;
      "--temp") set -- "$@" "-t" ;;
      "--norm")   set -- "$@" "-n" ;;
      "--float")   set -- "$@" "-f" ;;   
      "--help")   set -- "$@" "-h" ;;
       *)        set -- "$@" "$arg"
  esac
done

INPUT_HICPRO=""
INPUT_BED=""
NORMALIZE=0
FLOAT=0
NPROC=1
CHROMSIZES_FILE=""
RES=10000
OUT="./"
TEMP="./"

while getopts ":i:b:c:p:r:o:t:nfh" OPT
do
    case $OPT in
	i) INPUT_HICPRO=$OPTARG;;
	b) INPUT_BED=$OPTARG;;
	n) NORMALIZE=1;;
	f) FLOAT=1;;
	c) CHROMSIZES_FILE=$OPTARG;;
	p) NPROC=$OPTARG;;
	r) RES=$OPTARG;;
	o) OUT=$OPTARG;;
	t) TEMP=$OPTARG;;
	h) help ;;
	\?)
	    echo "Invalid option: -$OPTARG" >&2
	    usage
	    exit 1
	    ;;
	:)
	    echo "Option -$OPTARG requires an argument." >&2
	    usage
	    exit 1
	    ;;
    esac
done

if [[ -z $INPUT_HICPRO ]];
then
    usage
    exit
fi

if [[ ! -e $CHROMSIZES_FILE ]]; then
    echo -e "$CHROMSIZES_FILE file not found. Exit"
    exit 1
fi

## Detect input data type
DATATYPE=""
if [[ $INPUT_HICPRO == *.mat* ]]; then
    DATATYPE="MATRIX"
elif [[ $INPUT_HICPRO == *allValidPairs* || $INPUT_HICPRO == *validPairs* ]]; then
    DATATYPE="VALID"
else
    echo -e "Unknown input data type. Expect .matrix or _allValidPairs input files."
    exit 1
fi
echo -e "$DATATYPE input file detected ..."

## Check cooler version
which cooler > /dev/null;
if [ $? != "0" ]; then
    echo -e "Cooler is not installed or is not in your $PATH. See https://github.com/mirnylab/cooler for details."
    exit 1;
fi

COOLER_VERSION=$(cooler --version 2>&1 | awk '{print $NF}')
echo "Cooler version $COOLER_VERSION detected ..."
cres=$(vercomp ${COOLER_VERSION} "0.7.6")
if [[ $cres == "2" ]]; then
    echo "Cooler version must be >= 0.7.6 ! Stop."
    exit 1
fi

if [[ $DATATYPE == "VALID" ]]; then
    which pairix > /dev/null;
    if [ $? != "0" ]; then
	echo -e "Pairix is not installed or is not in your PATH. See https://github.com/4dn-dcic/pairix."
	exit 1;
    fi
fi

echo -e "\nGenerating .cool files ..."
tmp_dir=${TEMP}/_tmp$$
mkdir -p $tmp_dir

if [[ $DATATYPE == "MATRIX" ]]; then
    out=$(basename $INPUT_HICPRO | sed -e 's/.mat.*/.cool/')
    
    cooler makebins $CHROMSIZES_FILE $RES > $tmp_dir/bins.bed
    if [ $FLOAT == 1 ];
    then
	cooler load --count-as-float -f coo --one-based $tmp_dir/bins.bed $INPUT_HICPRO $tmp_dir/$out
    else
	cooler load -f coo --one-based $tmp_dir/bins.bed $INPUT_HICPRO $tmp_dir/$out
    fi
    echo -e "\nZoomify .cool file ..."
    if [[ $NORMALIZE == 1 ]]; then
	cooler zoomify --nproc ${NPROC} --balance $tmp_dir/$out
    else
	cooler zoomify --nproc ${NPROC} $tmp_dir/$out
    fi
    out=$(basename $INPUT_HICPRO | sed -e 's/.mat.*/.mcool/')
    
elif [[ $DATATYPE == "VALID" ]]; then
    out=$(basename $INPUT_HICPRO | sed -e 's/.allValidPairs.*/.cool/')

    awk '{OFS="\t";print $2,$3,$4,$5,$6,$7,1}' $INPUT_HICPRO | sed -e 's/+/1/g' -e 's/-/16/g' > $tmp_dir/contacts.txt
    cooler csort --nproc ${NPROC} -c1 1 -p1 2 -s1 3 -c2 4 -p2 5 -s2 6 \
	   -o $tmp_dir/contacts.sorted.txt.gz  \
	   $tmp_dir/contacts.txt \
	   $CHROMSIZES_FILE
    
    cooler makebins $CHROMSIZES_FILE $RES > $tmp_dir/bins.bed
    cooler cload pairix --nproc ${NPROC} $tmp_dir/bins.bed $tmp_dir/contacts.sorted.txt.gz $tmp_dir/$out

    echo -e "\nZoomify .cool file ..."
    if [[ $NORMALIZE == 1 ]]; then
	cooler zoomify --nproc ${NPROC} --balance $tmp_dir/$out
    else
	cooler zoomify --nproc ${NPROC} $tmp_dir/$out
    fi
    out=$(basename $INPUT_HICPRO | sed -e 's/.allValidPairs.*/.mcool/')
fi

## mv to out
mv $tmp_dir/*cool ${OUT}/

## clean
/bin/rm -rf $tmp_dir

echo -e "\nCooler file generated with success ..."
echo "Please copy the file $out in your Higlass input directory and run :"
echo "sudo docker exec higlass-container python higlass-server/manage.py  ingest_tileset --filename /tmp/$out --datatype matrix --filetype cooler" 



