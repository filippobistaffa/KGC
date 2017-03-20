#!/usr/bin/env bash

red='\033[0;31m'			# Red
nc='\033[0m'				# No color
re='^[0-9]+$'				# Regular expression to detect natural numbers
k=5					# Default k = 5

usage() { echo -e "Usage: $0 -i <filename> [-k <max_card>] [-c]\n-i\tInput filename\n-s\tMaximum cardinality (optional, default k = 5)\n-c\tEnable CSV output (optional)" 1>&2; exit 1; }

while getopts ":i:k:c" o; do
	case "${o}" in
	i)
		i=${OPTARG}
		if [ ! -f "$i" ]
		then
			echo -e "${red}Input file \"$i\" not found!${nc}\n" 1>&2
			usage
		fi
		;;
	k)
		k=${OPTARG}
		if ! [[ $k =~ $re ]] ; then
			echo -e "${red}Parameter k must be a number!${nc}\n"
			usage
		fi
		;;
	c)
		c=1
		;;
	\?)
		echo -e "${red}-$OPTARG is not a valid option!${nc}\n" 1>&2
		usage
		;;
	esac
done
shift $((OPTIND-1))

if [ -z "${i}" ]
then
	echo -e "${red}Missing input file!${nc}\n" 1>&2
	usage
fi

e=`cat "$i" | wc -l`
n=`grep -v " " "$i" | wc -l`
e=$(( $e - $n ))

tmp=`mktemp`
echo "#define N $n" > $tmp
echo "#define E $e" >> $tmp
echo "#define K $k" >> $tmp

if [ ! -z "${c}" ]
then
	echo "#define CSV" >> $tmp
fi

if [ ! -f instance.h ]
then
	mv $tmp "instance.h"
else
	md5a=`md5sum instance.h | cut -d\  -f 1`
	md5b=`md5sum $tmp | cut -d\  -f 1`

	if [ $md5a != $md5b ]
	then
		mv $tmp "instance.h"
	else
		rm $tmp
	fi
fi

make -j
if [[ $? == 0 ]]
then
	bin=$0
	bin=${bin%???}
	$bin $i
fi
