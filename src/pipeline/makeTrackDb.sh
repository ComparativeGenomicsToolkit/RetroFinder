#
source $1
echo "track ucscRetroAli$VERSION" 
echo "shortLabel Retroposed Genes $VERSION.0" 
echo "longLabel Retroposed Genes, Including Pseudogenes $RUNDATE [ucscRetroAli$VERSION]" 
echo "group genes" 
echo "type psl" 
echo "priority 37.14" 
echo "color 20,0,250" 
echo "visibility pack" 
echo "nextItemButton on" 
echo "ucscRetroInfo ucscRetroInfo$VERSION" 
echo "baseColorDefault diffCodons" 
echo "baseColorUseCds table ucscRetroCds$VERSION" 
echo "baseColorUseSequence extFile ucscRetroSeq$VERSION ucscRetroExtFile$VERSION" 
echo "indelDoubleInsert on" 
echo "indelQueryInsert on" 
echo "showDiffBasesAllScales ." 
echo "showDiffBasesMaxZoom 10000.0" 
echo "showCdsAllScales ." 
echo "showCdsMaxZoom 10000.0" 
echo "" 
echo "searchName ucscRetroInfoRefSeq${VERSION}" 
echo "searchTable ucscRetroAli${VERSION}" 
echo "searchDescription Retroposed GenesV${VERSION}, Including Pseudogenes - $RUNDATE" 
echo "query select tName, tStart,tEnd, qName from %s where qName like '%s%%'" 
echo "xrefTable refLink, ucscRetroInfo${VERSION}" 
echo "dontCheckXrefQueryFormat 1" 
echo "xrefQuery select ucscRetroInfo${VERSION}.name, refLink.name from %s where refLink.name like '%s%%' and refSeq = mrnaAcc " 
echo "searchPriority 3.52" 
echo "" 
echo "searchName ucscRetroInfoMrna${VERSION}" 
echo "searchTable ucscRetroAli${VERSION}" 
echo "searchDescription Retroposed GenesV${VERSION}, Including Pseudogenes - $RUNDATE" 
echo "query select tName, tStart,tEnd, qName from %s where qName like '%s%%'" 
echo "searchPriority 3.55" 
echo "" 
echo "searchName ucscRetroUniProt${VERSION}" 
echo "searchTable ucscRetroAli${VERSION}" 
echo "searchDescription Retroposed GenesV${VERSION}, Including Pseudogenes - $RUNDATE" 
echo "query select tName, tStart,tEnd, qName from %s where qName like '%s%%'" 
echo "dontCheckXrefQueryFormat 1" 
echo "xrefTable kgXref, ucscRetroInfo${VERSION}" 
echo "xrefQuery select ucscRetroInfo${VERSION}.name, spDisplayID from %s where spDisplayID like '%s%%' and kgName = kgID " 
echo "searchPriority 3.54" 
echo "" 
echo "searchName ucscRetroKnownGene${VERSION}" 
echo "searchTable ucscRetroAli${VERSION}" 
echo "searchDescription Retroposed GenesV${VERSION}, Including Pseudogenes - $RUNDATE" 
echo "query select tName, tStart,tEnd, qName from %s where qName like '%s%%'" 
echo "dontCheckXrefQueryFormat 1" 
echo "xrefTable kgXref, ucscRetroInfo${VERSION}" 
echo "xrefQuery select ucscRetroInfo${VERSION}.name, geneSymbol from %s where geneSymbol like '%s%%' and kgName = kgID " 
echo "searchPriority 3.53" 
echo "" 