#!/bin/bash

# quakepy/pmc/pmcmap.bash
# $Id: pmcmap.bash 328 2012-03-28 15:08:30Z tyrone $
#
# The QuakePy package
# http://www.quakepy.org
#

############################################################################
#    Copyright (C) 2007-2011 by Danijel Schorlemmer & Fabian Euchner       #
#    fabian@fabian-euchner.de                                              #
#                                                                          #
#    This program is free software; you can redistribute it and#or modify  #
#    it under the terms of the GNU General Public License as published by  #
#    the Free Software Foundation; either version 2 of the License, or     #
#    (at your option) any later version.                                   #
#                                                                          #
#    This program is distributed in the hope that it will be useful,       #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#    GNU General Public License for more details.                          #
#                                                                          #
#    You should have received a copy of the GNU General Public License     #
#    along with this program; if not, write to the                         #
#    Free Software Foundation, Inc.,                                       #
#    59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             #
############################################################################

# Generates a map of probability-based completeness

function PrintHelp()
{
  echo "Creates PMC maps with GMT"
  echo "Usage: pmcmap.bash [OPTIONS]"
  echo "  Select data"
  echo "    -t <value>   Type of map"
  echo "                   m: Magnitude of completeness map (default)"
  echo "                   p: Probability map"
  echo "    -A <value>   Region for map, sets defaults for -a, -j, -B, -G"
  echo "                   SC: Southern California"
  echo "                   IT: Italy"
  echo "                   CH: Switzerland"
  echo "                   JP: Japan"
  echo "                   NZ: New Zealand"
  echo "                   NC: Northern California"
  echo "    -d <file>    Data file"
  echo "    -r <file>    Topography file (default: topo.grd)"
  echo "    -R <file>    Topography intensity file (default: topo.intens.grd)"
  echo "    -o <file>    PMC povray color table (default: pmc.pov)"
  echo "    -O <file>    Probability povray color table (default: pmc.prob.pov)"
  echo "    -m <value>   Target magnitude for probability map"
  echo "    -p <value>   Target probability for completeness magnitude map"
  echo "    -e <value>   Target depth for completeness magnitude map"
  echo "    -U <value>   Magnitude maximum cutoff (default: 4.0)"
  echo "    -u <value>   Magnitude minimum cutoff (default: 0.0)"
  echo "  Plot options"
  echo "    -a <value>   Plot area (GMT format)"
  echo "    -j <value>   Projection (GMT format)"
  echo "    -M <value>   Frame settings (GMT format)"
  echo "    -B <value>   Color bar position (GMT format)"
  echo "    -G <value>   Color bar legend position (GMT format)"
  echo "    -y <layers>  Plot layers"
  echo "                   s: stations"
  echo "                   c: contour lines (defined with -c)"
  echo "                   f: faults"
  echo "    -s <file>    Station file"
  echo "    -c <levels>  Define contour levels <value1>[A]/<value2>[A]/... (default: 0.99A/0.999A)"  
  echo "    -F <file>    Fault file"
  echo "    -b <file>    Boundary file"
  echo "    -D           Plot date label"
  echo "    -T           Do not plot topography"
  echo "  Special options"
  echo "    -f <value>   Type of figure"
  echo "                   e: EPS Figure (default)"
  echo "                   r: Raster image PNG"
  echo "                   p: Preview (PNG) for website"
  echo "    -i <file>    Input (= parameter) filename"
  echo "    -z <file>    Target filename"
  echo "    -Z           Leave target file in output directory"
  echo "    -S <dir>     Output directory"
  echo "    -x <dir>     Path of GMT executables"
  echo "    -P <file>    Absolute file name of Python interpreter (default: /usr/bin/python)"
  echo "    -C <file>    Absolute file name of convert executable (default: /usr/bin/convert)"
  echo "    -E <file>    Absolute file name of eps2eps executable (default: /usr/bin/eps2eps)"
  echo "    -L <dir>     LD_LIBRARY_PATH environment setting for netCDF library (default: /usr/local/netCDF/lib/)"
  echo "    -g           Do garbage collection after run"
  echo "    -V           Verbose mode"
  echo "    -h           Print this message and exit"
  return
}

# set data and scratch directory
DATA='.'
SKR="./$$"

# Set defaults
date_label=0
garbage_collection=0
plot_topo=1
verbose_mode=0
leave_final_file_in_skr=0

mag_cut_min=0.0
mag_cut_max=4.0
target_magnitude=2.4
target_probability=0.999
target_depth='None'

pmc_pov_file="$DATA/pmc.pov"
prob_pov_file="$DATA/pmc.prob.pov"

topo_grid_file="$DATA/topo.grd"
intens_grid_file="$DATA/topo.intens.grd"

type_figure="e"
type_map="m"

gmt_path=""
input_file='None'

#data_file='None'
data_file="$DATA/map.gmt.dat"

logfile='None'
final_file='None'

boundary_file='None'
fault_file='None'

region_code=""

stationfile="$DATA/station.gmt.dat"
contour_level='0.99A/0.999A'
plot_layers=""

python_exec="/usr/bin/python"
convert_exec="/usr/bin/convert"
eps2eps_exec="/usr/bin/eps2eps"

LIBPATH="/usr/local/netCDF/lib/"

frame_options="--PS_LINE_JOIN=miter"
white_frame_opt="--BASEMAP_FRAME_RGB=255/255/255"
black_frame_opt="--BASEMAP_FRAME_RGB=0/0/0"

#colorbar_font_options="16 0 0 CM"
colorbar_font_options="16 0 21 CM"
colorbar_line_options="--PS_LINE_CAP=square --PS_LINE_JOIN=miter --TICK_PEN=1p"

## set defaults for Southern California

area_code="-R-122/-113.5/31.5/38" # -R-122/-113.5/31.5/38 # -R5.0/11.5/45.25/48.25
projection_code="-JM6i"
frame="-Ba2f1/a2f1WeSn"

colorbar_position="-d5.23i/4i/2i/0.2i"
colorbar_position="-D5.23i/4i/2i/0.2i"
colorbar_legend_position="-114.15 37.72"
colorbar_frame_m="-Bf0.5a1";
coastline_width=2
political_boundary_width=3
datelabel_position="-121.0 37.0"
         
# Parse options
while getopts "DghTVZa:A:b:B:c:C:d:e:E:f:F:G:i:j:l:L:m:M:o:O:p:P:r:R:s:S:t:u:U:x:y:z:" flag
do
  case $flag in
    D ) date_label=1;;
    g ) garbage_collection=1;;
    h ) PrintHelp; exit;;
    T ) plot_topo=0;;
    V ) verbose_mode=1;;
    Z ) leave_final_file_in_skr=1;;
    a ) area_code=$OPTARG;;
    A ) region_code=$OPTARG;;
    b ) boundary_file=$OPTARG;;
    B ) colorbar_position=$OPTARG;;
    c ) contour_level=$OPTARG;;
    C ) convert_exec=$OPTARG;;
    d ) data_file=$OPTARG;;
    e ) target_depth=$OPTARG;;
    E ) eps2eps_exec=$OPTARG;;
    f ) type_figure=$OPTARG;;
    F ) fault_file=$OPTARG;;
    G ) colorbar_legend_position=$OPTARG;;
    i ) input_file=$OPTARG;;
    j ) projection_code=$OPTARG;;
    l ) logfile=$OPTARG;;
    m ) target_magnitude=$OPTARG;;
    M ) frame=$OPTARG;;
    L ) LIBPATH=$OPTARG;;
    o ) pmc_pov_file=$OPTARG;;
    O ) prob_pov_file=$OPTARG;;
    P ) python_exec=$OPTARG;;
    p ) target_probability=$OPTARG;;
    r ) topo_grid_file=$OPTARG;;
    R ) intens_grid_file=$OPTARG;;
    s ) stationfile=$OPTARG;;
    S ) SKR=$OPTARG;;
    t ) type_map=$OPTARG;;
    u ) mag_cut_min=$OPTARG;;
    U ) mag_cut_max=$OPTARG;;
    x ) gmt_path=$OPTARG;;
    y ) plot_layers=$OPTARG;;
    z ) final_file=$OPTARG;;
  esac
done

# read parameters from file
if [ -s $input_file ]; then

    # loop through each line of input file
    while read line
    do
        # get the parameter and value from the line
        param=`echo $line | cut -d = -f 1`
        value=`echo $line | cut -d = -f 2`

        #set the variable
        eval ${param}=\$value

    done < <(cat $input_file)
fi

# set path for libnetcdf
export LD_LIBRARY_PATH=$LIBPATH:$LD_LIBRARY_PATH

# sanity checks
if [ ! -s $data_file ]; then
    echo "pmcmap error: no valid data file at $data_file" 1>&2
    exit
fi

if [ $final_file == 'None' ]; then
    case $type_figure in
        p | r ) final_file="map.png";;
        * )     final_file="map.eps";;
    esac
fi
    
# Create output directory
if [ ! -d $SKR ]; then
  # output directory is newly created for this run, safe to remove
  mkdir $SKR
  new_skr_dir=1
else
  # output directory already exists, do not remove in garbage collection
  new_skr_dir=0
fi

# Set logfile name
if [ $logfile == 'None' ]; then
    logfile="$SKR/map.gmt.log"
fi

# create empty logfile if not already existing
if [ ! -w $logfile ]; then
    touch $logfile
fi

# start logging
echo "START LOGGING" >> $logfile

# internal file names, in output directory
output="$SKR/out.gmt.ps"
output_converted="$SKR/out.gmt.eps"
output_eps="$SKR/map.gmt.eps"
output_png="$SKR/map.gmt.png"

datafile_gmt="$SKR/datafile.gmt.dat"
datafile_gawk="$SKR/datafile.gawk.dat"

datafile_grid="$SKR/datafile.grd"
datafile_hires="$SKR/datafile.hires.grd"

contours="$SKR/contours.dat"
colormap="$SKR/colormap.cpt"

# get presets for specific regions

case $region_code in

         # defaults for Northern California
    NC ) frame="-Ba2f1/a2f1WeSn";
         map_lon_min=-127.0;
         map_lon_max=-115.0;
         map_lat_min=33.0;
         map_lat_max=44.0;
         area_code="-R$map_lon_min/$map_lon_max/$map_lat_min/$map_lat_max";
         projection_code="-JM6i";

         colorbar_position="-D5.3i/3.95i/2i/0.2i";
         colorbar_legend_position="-114.18 37.72";
         colorbar_frame_m="-Bf0.5a1";
         coastline_width=2;
         political_boundary_width=3;
         datelabel_position="-121.0 37.0";;

         # defaults for Southern California
    SC ) frame="-Ba2f1/a2f1WeSn";
         map_lon_min=-122.0;
         map_lon_max=-113.5;
         map_lat_min=31.5;
         map_lat_max=38.0;
         area_code="-R$map_lon_min/$map_lon_max/$map_lat_min/$map_lat_max";
         projection_code="-JM6i";

         colorbar_position="-D5.3i/3.95i/2i/0.2i";
         colorbar_legend_position="-114.18 37.72";
         colorbar_frame_m="-Bf0.5a1";
         coastline_width=2;
         political_boundary_width=3;
         datelabel_position="-121.0 37.0";;

         # defaults for Italy
    IT ) frame="-Ba2f1/a2f1WeSn";
         map_lon_min=4.0;
         map_lon_max=21.0;
         map_lat_min=34.0;
         map_lat_max=49.0;
         area_code="-R$map_lon_min/$map_lon_max/$map_lat_min/$map_lat_max";
         projection_code="-JM6i";

         colorbar_position="-D5.23i/5.4i/2i/0.2i";
         colorbar_legend_position="19.45 48.3";
         colorbar_frame_m="-Bf0.5a1";
         coastline_width=2;
         political_boundary_width=3;
         datelabel_position="5.0 21.0";;

         # defaults for Switzerland
    CH ) frame="-Ba2f1/a2f1WeSn";
         map_lon_min=5.0;
         map_lon_max=11.5;
         map_lat_min=45.25;
         map_lat_max=48.25;
         area_code="-R$map_lon_min/$map_lon_max/$map_lat_min/$map_lat_max";
         projection_code="-JM6i";

#         colorbar_position="-D0.27i/2.45i/1.8i/0.18i";
         colorbar_position="-D0.18i/2.45i/1.8i/0.18i";
         colorbar_legend_position="5.518 48.05";
         colorbar_frame_m="-Bf0.5a0.5";
         coastline_width=2;
         political_boundary_width=3;
         datelabel_position="5.0 21.0";;

         # defaults for Japan
    JP ) frame="-Ba8f1/a8f1WeSn";
         map_lon_min=120.0;
         map_lon_max=150.0;
         map_lat_min=20.0;         
         map_lat_max=49.0;
         area_code="-R$map_lon_min/$map_lon_max/$map_lat_min/$map_lat_max";
         projection_code="-JM6i";

         colorbar_position="-D5.2i/1.22i/2i/0.2i";
         colorbar_legend_position="147.2 31.8";
         colorbar_frame_m="-Bf0.5a1";
         coastline_width=2;
         political_boundary_width=3;
         datelabel_position="124.0 48.0";;

         # defaults for New Zealand
    NZ ) frame="-Ba5f1/a5f1WeSn";
         map_lon_min=163.0;
         map_lon_max=182.0;
         map_lat_min=-50.0;         
         map_lat_max=-33.0;
         area_code="-R$map_lon_min/$map_lon_max/$map_lat_min/$map_lat_max";
         projection_code="-JM6i";

         colorbar_position="-D0.4i/5.6i/2i/0.2i";
         colorbar_legend_position="165.0 -33.6";
         colorbar_frame_m="-Bf0.5a1";
         coastline_width=2;
         political_boundary_width=3;
         datelabel_position="165.0 -49.5";;
esac

# option -V: verbose (disable for prodution code)
if [ $verbose_mode -eq 1 ]; then
  verbose_flag="-V"
else
  verbose_flag=""
fi

# set general GMT options and frame setting
gmtopts="$verbose_flag $area_code $projection_code"
frame_parameter="$frame $frame_options"

# GMT settings

#${gmt_path}gmtdefaults -L > .gmtdefaults4

# ${gmt_path}gmtset PLOT_DEGREE_FORMAT = -dddF
${gmt_path}gmtset PLOT_DEGREE_FORMAT = dddF
${gmt_path}gmtset BASEMAP_TYPE       = fancy
${gmt_path}gmtset BASEMAP_FRAME_RGB  = 0/0/0
${gmt_path}gmtset ANOT_FONT_SIZE     = 16
${gmt_path}gmtset FRAME_PEN          = 5
${gmt_path}gmtset COLOR_BACKGROUND   = 255/255/255
${gmt_path}gmtset COLOR_FOREGROUND   = 0/0/0
${gmt_path}gmtset TICK_PEN           = 0.5p,0/0/0
#${gmt_path}gmtset X_ORIGIN           = 0.0c
#${gmt_path}gmtset Y_ORIGIN           = 0.0c

#${gmt_path}gmtdefaults -L >> $logfile

# check if input data is already in GMT format or if it is still (gzipped) XML
is_gzipped=`file $data_file | grep -c 'gzip' -`

# set options for pmc2gmt.py script
if [ $type_map == "m" ]; then
    pmc2gmt_opt="-p $target_probability"
else
    pmc2gmt_opt="-m $target_magnitude"
fi

# get depth layer option
if [ $target_depth != "None" ]; then
    pmc2gmt_opt="$pmc2gmt_opt -d $target_depth"
fi

if [ $is_gzipped -eq 1 ]; then

    $python_exec ./pmc2gmt.py -z $pmc2gmt_opt < $data_file > $datafile_gmt
    echo "converted gzipped XML input data, wrote to $datafile_gmt" >> $logfile
else

    # data_file is not gzipped, check if it's XML
    is_xml=`file $data_file | grep -c 'XML' -`

    if [ $is_xml -eq 1 ]; then
        $python_exec ./pmc2gmt.py $pmc2gmt_opt < $data_file > $datafile_gmt
        echo "converted XML input data, wrote to $datafile_gmt" >> $logfile
    else
        cp -f $data_file $datafile_gmt
    fi
fi

# try to extract date string from filename
date_str=`echo $data_file | grep -P -o -e '\d{4}-\d{2}-\d{2}'`
echo "extracted date string from file name $date_str" >> $logfile

# Prepare data 
sed "s/NaN/${mag_cut_max}/g;s/nan/${mag_cut_max}/g" $datafile_gmt | gawk '{printf("%g\t%g\t%s\n", $1, $2, $3)}' > $datafile_gawk

# log some stuff
echo "current directory is " `pwd` >> $logfile
echo "writing output to file " $final_file >> $logfile
echo "temporary directory is " $SKR  >> $logfile

echo "preparing input data from $datafile_gmt" >> $logfile
ls -la $datafile_gmt >> $logfile

echo "prepared input data, wrote to $datafile_gawk" >> $logfile
ls -la $datafile_gawk >> $logfile
echo "------------------------------------------" >> $logfile

# Create colormap
if [ $type_map == 'm' ]; then
  colorbar_frame_opt=$white_frame_opt
#  $python_exec ./
  pov2cpt.py --min=${mag_cut_min} --max=${mag_cut_max} --extend --integer --nan=0/71/217 --output=$colormap $pmc_pov_file
else
#  $python_exec ./
  pov2cpt.py --min=0 --max=1 --extend --integer --output=$colormap $prob_pov_file
  colorbar_frame_opt=$black_frame_opt
fi
echo " -created color table" >> $logfile

# check if GMT command produces good output (on stderr)
# ${gmt_path}psbasemap 2>> $logfile

# Create basemap
${gmt_path}psbasemap --PLOT_DEGREE_FORMAT=dddF $gmtopts $frame_parameter -P -K > $output
echo " -created basemap with GMT options $gmtopts $frame_parameter" >> $logfile

# log size of output EPS file at beginning
ls -la $output >> $logfile

# make grid from data file
# set missing nodes (outside of polygon) to maximum value (like NaNs in ASCII data file): option -N
# Plot data in background to cover off-shore areas
${gmt_path}xyz2grd $datafile_gawk -G$datafile_grid $area_code -I0.1/0.1 -N${mag_cut_max}
echo " -created data grid" >> $logfile

# resample data grid to higher resolution (because of topography)
${gmt_path}grdsample $datafile_grid -G$datafile_hires $area_code -I0.5m/0.5m $verbose_flag

# make image from data grid
${gmt_path}grdimage $datafile_hires $gmtopts -C$colormap -O -K >> $output
echo " -created data grid image" >> $logfile

# Plot topography
if [ $plot_topo -eq 1 ]; then
    # Clip image by coastline
    ${gmt_path}pscoast $gmtopts -Dh -Gc -O -K >> $output

    # Plot data on top of topography
    ${gmt_path}grdview $topo_grid_file $gmtopts -G$datafile_hires -C$colormap -I$intens_grid_file -Qi600 -O -K >> $output

    # Release clipping mask
    ${gmt_path}pscoast $gmtopts -Dh -Q -O -K >> $output

    colorbar_illumination_options="-I0.5"
    echo " -plotted topography" >> $logfile
else
    
    colorbar_illumination_options=""
    echo " -no topography requested" >> $logfile
fi

# Plot coastlines
${gmt_path}pscoast $gmtopts -Dh -W$coastline_width/0/0/0 -N1/$political_boundary_width/0/0/0 -O -K >> $output
echo " -plotted coastline" >> $logfile

# Plot border
if [ -e $boundary_file ]; then
    ${gmt_path}psxy $boundary_file $gmtopts -A -MS -: -W$political_boundary_width/0/0/0 -O -K >> $output
    echo " -plotted boundaries from file $boundary_file" >> $logfile
else
    echo " -ERROR BOUNDARIES: cannot read boundary file $boundary_file" >> $logfile
fi

# Plot various layers (events, stations, contours, authorative region)
number_layers=`expr length $plot_layers - 1`

echo " -plotting" `expr $number_layers + 1` "layers from layer code: $plot_layers" >> $logfile
for i in `seq 0 $number_layers`
do
  case ${plot_layers:$i:1} in
    c ) echo " --plotting contours, contour level $contour_level" >> $logfile;
        echo $contour_level | sed -e 's/[0-9,.]*/&\t/g' -e 's/\//\n/g' > $contours;
        ${gmt_path}grdcontour $datafile_grid $gmtopts -W2/255/255/255-At -C$contours -O -K >> $output;;
    f ) echo " --plotting faults from file $fault_file" >> $logfile;
        ${gmt_path}psxy $fault_file $gmtopts -M -W2 -O -K >> $output;;
    s ) echo " --plotting stations from file $stationfile" >> $logfile;
        cat $stationfile | gawk '{printf("%e\t%e\n", $1, $2 )}' | ${gmt_path}psxy $gmtopts -St8p -G128 -W1/0/0/0 -O -K >> $output;;
  esac
done
echo " -layers finished" >> $logfile

# date label
if [ $date_label -eq 1 ]; then
  echo "$datelabel_position $colorbar_font_options $date_str" | ${gmt_path}pstext $gmtopts -S0.5p $colorbar_frame_opt -P -O -K >> $output
  echo " -plotted date label $date_str" >> $logfile
fi

pscoast $gmtopts -Dh -W1/0/0/0 -N1/1/0/0/0 -O -K >> $output

# Colorbar
# Note: this has to be the last plot command because of '-K' option of GMT
if [ $type_map == 'm' ]; then

  # use colorbar illumination (-I) only if topography is plotted
  ${gmt_path}psscale $colorbar_position -C$colormap $colorbar_frame_m $colorbar_illumination_options $colorbar_frame_opt $colorbar_line_options -O -K >> $output
  echo "$colorbar_legend_position $colorbar_font_options M@-P@-" | ${gmt_path}pstext $gmtopts -S0.5p $colorbar_frame_opt -P -O >> $output
else
  ${gmt_path}psscale $colorbar_position -C$colormap -Bf0.25a0.5 $colorbar_illumination_options $colorbar_frame_opt $colorbar_line_options -O -K >> $output
  echo "$colorbar_legend_position $colorbar_font_options P@-${target_magnitude}@-" | ${gmt_path}pstext $gmtopts -S0.5p $colorbar_frame_opt -P -O >> $output
  echo " --probability map: target magnitude is ${target_magnitude}" >> $logfile
fi
echo " -plotted color bar" >> $logfile

# log size of output EPS file at end
ls -la $output >> $logfile

# Clean EPS
if [ $type_figure == 'e' ]; then
  ps2raster $output -Te -A
#  $eps2eps_exec $output $SKR/$final_file
#  rm -f $output
  mv $output_converted $SKR/$final_file
  echo " --converted figure to final EPS" >> $logfile;
else
  ps2raster $output -Te -A
#  $eps2eps_exec $output $output_eps
#  rm -f $output
  mv $output_converted $output_eps
  echo " --converted figure to temporary EPS" >> $logfile;
fi

case $type_figure in
  p ) echo " --converting EPS figure to preview PNG" >> $logfile;
      $convert_exec -density 144 $output_eps $output_png;
      $convert_exec -resize 600x $output_png $SKR/$final_file;;
  r ) echo " --converting EPS figure to full size PNG" >> $logfile;
      $convert_exec -density 288 $output_eps $SKR/$final_file;;
esac

if [ $leave_final_file_in_skr -eq 0 ]; then
    mv $SKR/$final_file $final_file
    echo " --moved final file out of temp dir" >> $logfile;
fi

echo "END LOGGING" >> $logfile

# garbage collection
if [ $garbage_collection -eq 1 ] && [ $new_skr_dir -eq 1 ]; then
    rm -R $SKR
fi


