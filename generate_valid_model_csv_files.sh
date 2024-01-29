#!/bin/bash
#
#   This script creates a list of the valid cmip5 and cmip6 models, based on the data files it finds in 
#   $cmip5_dir and $cmip6_dir
#
#   cmip5_models.csv and cmip6_models.csv each contain columns:
#       model calendar nlats nlons
#           (1 line for each model)
#
#   valid_models_cmip5.csv and valid_models_cmip6.csv contain columns:
#
#       varname model scenario ensemble [grgn] calendar nlats nlons start_year end_year
#
#           (1 line for each unique varname/model/scenario/ensemble/[grgn] 
#
#           ("grgn" is the gridding code (gn=native gridding, gr=regridded), which is only present in the cmip6 data.)
#
#   These files are read by various ARRM_V2 programs to check input against expected values.
#
#   To skip cmip5 or cmip6, add -s cmip5  or -s cmip6
#
#   To skip adding start- and end-dates, use -d  (needed when working with multiple un-joined files in the directories.
#
#   filenames which don't parse as expected (where start_year or end_year aren't all integers) are skipped, 
#   and their names written to cmip[56]_badfilenames.txt
#
#   NOTE:  "calendar" string may be modified from the actual calendar in the file:
#
#               standard        replaces "gregorian" and proleptic_gregorian"
#               365_day         replaces "noleap" and "365-day"
#
#           this is done because some of the models use a different name for the calendar in different files,
#           ARRM_V2 uses calendar_length(calendar) so this should not be an issue in ARRM_v2 code.  If you're using
#           the .csv files for some other purpose, see comments below where the calendar is extracted via ncdump

usage()
{
    prog=`basename $0`
    mypath=`readlink -f \`dirname $0 \` `
                cat << EOF
usage: $prog [-h] [-s cmip5|cmip6] [-d] [cmip5_dir  [cmip6_dir] ]   

    -h            help (this output)
    -s cmip[5|6]  skip running for cmip5 or cmip6
    -d            skip finding start and end dates.

	finds all valid model combinations for cmip5 and cmip6 model data
	by finding all netcdf files in the specified folders, parsing the filenames for
	model, scenario, ensemble, variable, etc. and extracts metadata using ncdump and grep. 

	Use -d to skip adding start and end dates when working with multiple raw, un-joined netcdf files for each set.

	cmip5_dir must be present to specify cmip5_dir.  If skipping cmip5 files, just use "skip" for cmip5_dir.

	creates two files for each:

	cmip[56]_models.csv 		with model, calendar, #lats, #lons
	valid_models_cmip[56].csv	with varname,model,scenario,ensemble,calendar,nlats,nlons,start_year,end_year
					(start_ and end_year are optional)

        NOTE:  if there are too many netcdf files, this script will fail.  "Too many"  depends on how many chars are allowed on a line by the shell.
               In that case, just join the files first, and work with the joined files.
    (run from $mypath)
EOF
        exit -1
}

myhost=$(hostname)

#..............parse options on command line
while getopts "s:hd" OPTION
do
    case $OPTION in
        h)
            usage
            exit -1
            ;;
	d)  skip_dates=1      ;;
        s)  skips=$OPTARG     ;;
        *)  
            echo "bad argument: $OPTION"
            usage
            exit -1
            ;;
    esac
done

shift $((OPTIND - 1))

cmip5_dir=$1
cmip6_dir=$2

# echo cmip5_dir $cmip5_dir
# echo cmip6_dir $cmip6_dir

	# default folders to search if not provided on command line
if [[ -z $cmip5_dir ]]
then	
if [[ $myhost == neys.ttu.edu || $myhost == icsf-kmac.ttu.edu ]]
    then
        cmip5_dir="/Volumes/lacie_1/data/gcm_cmip5/rotated_netcdf"
    else
        cmip5_dir="/lustre/research/hayhoe/cmip5_rotated"
    fi
fi

if [[ -z $cmip6_dir ]]
then	
if [[ $myhost == neys || $myhost == neys.ttu.edu || $myhost == icsf-kmac.ttu.edu  ]]
    then
#        cmip6_dir="/Volumes/lacie_1/data/gcm_cmip6/rotated_netcdf"
        cmip6_dir="/Volumes/lacie_1/data/gcm_cmip6/tllc"
    else
#        cmip6_dir="/lustre/research/hayhoe/cmip6_rotated"
        cmip6_dir="/lustre/research/hayhoe/cmip6_tllc"
    fi
fi

f5name="valid_models_cmip5.csv"
f6name="valid_models_cmip6.csv"
m5name="cmip5_models.csv"
m6name="cmip6_models.csv"
m5badnames="cmip5_badfilenames.txt"
m6badnames="cmip6_badfilenames.txt"

echo "cmip5_dir: $cmip5_dir"
echo "cmip6_dir: $cmip6_dir"
echo "f5name:    $f5name $m5name"
echo "f6name:    $f6name $m6name";
echo

#		regular expression to test if start_year or end_year are all integers.
#		if not, skip that filename.
re='^[0-9]+$'

    # get cmip5 stuff
if [[ "$skips" == "cmip5" ]]
then
    echo "skipping cmip5"
elif [[ ! -d $cmip5_dir ]]
then
    echo "can't find cmip5's llt dir"
    exit 1
else

    if [[ -f tmpfile1 ]]
    then
        rm tmpfile1
    fi
    if [[ -f tmpfile2 ]]
    then
        rm tmpfile2
    fi

    if [[ $skip_dates ]]
    then
        echo "varname,model,scenario,ensemble,calendar,nlats,nlons" > $f5name
        fns=$(find $cmip5_dir -name "*.nc" -print | fgrep -v ".US." | fgrep -v ".future." | sort )
    else
        echo "varname,model,scenario,ensemble,calendar,nlats,nlons,start_year,end_year" > $f5name
        fns=$(find $cmip5_dir -name "*.llt.nc" -print | fgrep -v ".US." | fgrep -v ".future." | sort )
    fi
    > $m5badnames 		
    nbad=0
    echo "model,calendar,nlats,nlons" > $m5name
    last_set="";
    last_file="";
    for fn in $fns
    do
        bn=$(basename $fn)
        varname=$(  echo $bn | awk 'BEGIN {FS="."; OFS=","} {print $3}' )
        model=$(    echo $bn | awk 'BEGIN {FS="."; OFS=","} {print $1}' )
        scenario=$( echo $bn | awk 'BEGIN {FS="."; OFS=","} {print $4}' )
        ensemble=$( echo $bn | awk 'BEGIN {FS="."; OFS=","} {print $2}' )
#       grgn=$(     echo $bn | awk 'BEGIN {FS="."; OFS=","} {print $2}' )  # cmip5 files don't have a grgn.
        syear=$(    echo $bn | awk 'BEGIN {FS="."; OFS=","} {print $6}' )
        eyear=$(    echo $bn | awk 'BEGIN {FS="."; OFS=","} {print $7}' )

		# skip this file if it matches last one and we're not doing dates.
	if [[ $skip_dates ]]		
	then
		cur_set="${varname},${model},${scenario},${ensemble}"
		if [[ "$cur_set" == "$last_set" ]]
		then
			echo "skipping matching set:  $bn $last_file"
			continue
		fi
		last_file=$bn
		last_set=$cur_set
	fi

	if [[ $syear =~ $re && $eyear =~ $re ]]
        then
            calendar=$(ncdump -h $fn | grep "time:calendar =" | awk '{print $3}' | tr -d '"'| sed "s/proleptic_gregorian/standard/" | sed "s/gregorian/standard/" | sed "s/noleap/365-day/" | tr "-" "_" )
                # note:  replace the line above with the line below to get the actual calendar if you need it: 
#           calendar=$(ncdump -h $fn | grep "time:calendar =" | awk '{print $3}' | tr -d '"' )
            nlats=$(   ncdump -h $fn | grep "lat = "          | awk '{print $3}'             )
            nlons=$(   ncdump -h $fn | grep "lon = "          | awk '{print $3}'             )
            if [[ "$nlons" == "UNLIMITED" ]]
            then
                nlons=$(   ncdump -h $fn | grep "lon = "        | awk '{print $6}' | tr -d '(' )
            fi

            if [[ $skip_dates ]]
            then
                echo "$varname,$model,$scenario,$ensemble,$calendar,$nlats,$nlons"              >> tmpfile1
            else
                echo "$varname,$model,$scenario,$ensemble,$calendar,$nlats,$nlons,$syear,$eyear" >> tmpfile1
            fi

            echo "$model,$calendar,$nlats,$nlons" >> tmpfile2
            printf "%-60s %-20s %-20s %10s %10s\n" $bn $model $calendar $nlats $nlons

        else
            echo "Bad cmip5 netcdf filename:  $bn"
	    echo $bn >> $m5badnames
	    nbad=$((nbad+1))
	fi



    done
    cat tmpfile1 | sort | uniq >> $f5name
 #  rm tmpfile1
    cat tmpfile2 | sort | uniq >> $m5name
 #  rm tmpfile2

    echo
    echo "cmip5 valid model sets written to $f5name"
    echo "cmip5 model list written to       $m5name"
    if (( nbad > 0 ))
    then
        echo "bad filenames written to          $m5badnames"
    fi
    echo
fi

        # get cmip6 stuff:

if [[ "$skips" == "skip_cmip6" ]]
then
    echo "skipping cmip6"
elif [[ ! -d $cmip6_dir ]]
then
    echo "can't find cmip6's tllc dir $cmip6_dir"
    exit 1
else

    if [[ -f tmpfile1 ]]
    then
        rm tmpfile1
    fi
    if [[ -f tmpfile2 ]]
    then
        rm tmpfile2
    fi

    if [[ $skip_dates ]]
    then
        echo "varname,model,scenario,ensemble,grgn,calendar,nlats,nlons" > $f6name
        fns=$(find $cmip6_dir -name "*.nc" -print | fgrep -v ".US." | fgrep -v ".future." | sort )
    else
        echo "varname,model,scenario,ensemble,grgn,calendar,nlats,nlons,start_year,end_year" > $f6name
        fns=$(find $cmip6_dir -name "*.tllc.nc" -print | fgrep -v ".US." | fgrep -v ".future." | sort )
    fi
 #  echo "varname,model,scenario,ensemble,grgn,calendar,nlats,nlons,start_year,end_year" > $f6name
    > $m6badnames
    nbad=0
    echo "model,calendar,nlats,nlons" > $m6name
    last_set="";
    last_file="";
    for fn in $fns
    do
        bn=$(basename $fn)
        varname=$(  echo $bn | tr "_" "." | awk 'BEGIN {FS="."; OFS=","} {print $1}' )
        model=$(    echo $bn | tr "_" "." | awk 'BEGIN {FS="."; OFS=","} {print $3}' )
        scenario=$( echo $bn | tr "_" "." | awk 'BEGIN {FS="."; OFS=","} {print $4}' )
        ensemble=$( echo $bn | tr "_" "." | awk 'BEGIN {FS="."; OFS=","} {print $5}' )
        grgn=$(     echo $bn | tr "_" "." | awk 'BEGIN {FS="."; OFS=","} {print $6}' )
        syear=$(    echo $bn | tr "_" "." | awk 'BEGIN {FS="."; OFS=","} {print $7}' | awk 'BEGIN {FS="-"; OFS=","} {print $1}' )
        eyear=$(    echo $bn | tr "_" "." | awk 'BEGIN {FS="."; OFS=","} {print $7}' | awk 'BEGIN {FS="-"; OFS=","} {print $2}' )

		# skip this file if it matches last one and we're not doing dates.
	if [[ $skip_dates ]]		
	then
		cur_set="${varname},${model},${scenario},${ensemble},${grgn}"
		if [[ "$cur_set" == "$last_set" ]]
		then
			echo "skipping matching set:  $bn $last_file"
			continue
		fi
		last_file=$bn
		last_set=$cur_set
	fi

	if [[ $syear =~ $re && $eyear =~ $re ]]
        then
            calendar=$(ncdump -h $fn | grep "time:calendar =" | awk '{print $3}' | tr -d '"'| sed "s/proleptic_gregorian/standard/" | sed "s/gregorian/standard/" | sed "s/noleap/365-day/" | tr "-" "_" )
                # note:  replace the line above with the line below to get the actual calendar if you need it: 
#       calendar=$(ncdump -h $fn | grep "time:calendar =" | awk '{print $3}' | tr -d '"' )
            nlats=$(   ncdump -h $fn | grep "lat = "          | awk '{print $3}'             )
            nlons=$(   ncdump -h $fn | grep "lon = "          | awk '{print $3}'             )
            if [[ "$nlons" == "UNLIMITED" ]]
            then
                nlons=$(   ncdump -h $fn | grep "lon = "        | awk '{print $6}' | tr -d '(' )
            fi
            syear=${syear:0:4}
            eyear=${eyear:0:4}

            if [[ $skip_dates ]]
            then
                echo "$varname,$model,$scenario,$ensemble,$grgn,$calendar,$nlats,$nlons"               >> tmpfile1
            else
                echo "$varname,$model,$scenario,$ensemble,$grgn,$calendar,$nlats,$nlons,$syear,$eyear" >> tmpfile1
            fi

#       echo "$varname,$model,$scenario,$ensemble,$grgn,$calendar,$nlats,$nlons,$syear,$eyear" >> tmpfile1
            echo "$model,$calendar,$nlats,$nlons" >> tmpfile2
            printf "%-80s %-20s %-20s %10s %10s\n" $bn         $model $calendar $nlats $nlons

        else
            echo "Bad cmip6 netcdf filename:  $bn"
	    echo $bn >> $m6badnames
	    nbad=$((nbad+1))
	fi


    done
    cat tmpfile1 | sort | uniq >> $f6name
    cat tmpfile2 | sort | uniq >> $m6name
    
    rm tmpfile1
    rm tmpfile2

    echo
    echo "cmip6 valid model sets written to $f6name"
    echo "cmip6 model list written to       $m6name"
    if (( nbad > 0 ))
    then
        echo "bad filenames written to          $m6badnames"
    fi
    echo
fi



