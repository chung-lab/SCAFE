#! /bin/sh

# Copyright 2011, 2012 Martin C. Frith

# This script simplifies the output of paraclu.  It omits clusters
# that are too long, or are singletons.  Then, it removes clusters
# whose fold-increase in density is too small.  Finally, it removes
# clusters that are contained in larger clusters.

# This script assumes that the input is sorted by start position
# (ascending) then end position (descending).

maxLengthDefault=200
minDensityRiseDefault=2

maxLength=$maxLengthDefault
minDensityRise=$minDensityRiseDefault
isSingleClusterIncrease=

while getopts hl:d:s opt
do
    case $opt in
	h)  cat <<EOF
Usage: $(basename $0) [options] [paraclu-output-file(s)]
Extract a subset of clusters from the output of paraclu.

Options:
  -h  show this help message and exit
  -l  maximum cluster length (default $maxLengthDefault)
  -d  minimum density increase (default $minDensityRiseDefault)
  -s  density increase applies to single clusters, not cumulatively
EOF
	    exit
	    ;;
	l)  maxLength=$OPTARG
	    ;;
	d)  minDensityRise=$OPTARG
	    ;;
	s)  isSingleClusterIncrease=1
	    ;;
    esac
done

shift $(($OPTIND - 1))

# Remove singleton clusters.
# Remove clusters longer than maxClusterLength.
awk '
$3 < $4 && $4 - $3 <= '$maxLength'
' "$@" |

if [ -n "$isSingleClusterIncrease" ]
then
    # Remove clusters less stable than minDensityRise.
    awk '$7 * '$minDensityRise' <= $8'
else
    # Remove clusters whose density is insufficiently above the local
    # baseline density.
    awk '
{n = $1" "$2}
n != o {e=-1; o=n}
$4 > e {e=$4; s=$7}
s * '$minDensityRise' <= $8 {print}
'
fi |

# Remove clusters that are contained in larger clusters.
awk '
{n = $1" "$2}
n != o {e=-1; o=n}
$4 > e {print; e=$4}
'
