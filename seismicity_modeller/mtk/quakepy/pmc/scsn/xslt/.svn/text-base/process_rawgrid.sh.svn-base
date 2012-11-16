#!/bin/bash

# $1 : Input XML (not prettified)
# $2 : Output XML (prettified)
# $3 : Output Data in GMT format

xmllint --format $1 > $2
xsltproc qpgrid2gmt.xsl $2 > $3


