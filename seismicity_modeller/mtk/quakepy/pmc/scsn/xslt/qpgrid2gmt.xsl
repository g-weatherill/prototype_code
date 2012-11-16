<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet	version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:output method="text" encoding="utf-8"/>
  
  <xsl:template match="QPGrid/grid/depthLayer/cell">
        <xsl:value-of select="@lon"/><xsl:text>    </xsl:text><xsl:value-of select="@lat"/>
        <xsl:apply-templates select="PMCData" />
  </xsl:template>
  
  <xsl:template match="PMCData">
    <!-- <xsl:param name="magnitude" select="probability" /> -->
    <xsl:for-each select = "mp">
        <!-- <xsl:text>    </xsl:text><xsl:value-of select="."/> -->
        <xsl:if test="./@probability = '0.999'">
            <xsl:text>    </xsl:text><xsl:value-of select="."/>
        </xsl:if>
    </xsl:for-each>
  </xsl:template>
  
</xsl:stylesheet>

