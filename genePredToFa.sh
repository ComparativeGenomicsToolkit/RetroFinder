#
source $1
chr=$2
getRnaPred -genomeSeqs=$GENOME/$DB/$DB.2bit -cdsUpper $DB ../transMapAlnMRna.gp $chr mrnaCds.${chr}.fa
tr acgt ACGT < mrnaCds.${chr}.fa > mrna.${chr}.fa
blat $GENOME/$DB/$DB.2bit mrna.${chr}.fa ../mrna.psl/${chr}.psl -q=rna -fine -mask=lower -ooc=$GENOME/$DB/11.ooc
