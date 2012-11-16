#!/bin/bash

# download monthly g'zipped QuakeML catalog chunks from GeoNet NZ server
# filename: nz.YYYY-MM.xml.gz

# please remove in order to run program
exit

SERVER_URL="http://magma.geonet.org.nz/services/quake/quakeml/1.0.1/query?"

for startyear in `seq 1950 2009`
do
    for startmonth in `seq 1 12`
    do
        if [ $startmonth -eq 12 ]; then
            endyear=`expr $startyear + 1`
            endmonth=1;
        else
            endyear=$startyear
            endmonth=`expr $startmonth + 1`;
        fi

        # format 2-digit month with trailing zero
        printf -v filename "nz.$startyear-%02d.xml" "$startmonth"
        url="${SERVER_URL}startDate=${startyear}-${startmonth}-01&endDate=${endyear}-${endmonth}-01&includePicksAndArrivals=true"

        wget -O $filename $url

        gzip -c $filename > "$filename.gz"
        rm -f $filename
    done
done
