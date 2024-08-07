#install using 'BiocManager::install("RcisTarget")'
library(RcisTarget)

#using  HALLMARK_G2M_CHECKPOINT as input gene list
geneList1 <- c("ABL1",	"AMD1",	"ARID4A",	"ATF5",	"ATRX",	"AURKA",	"AURKB",	"BARD1",	"BCL3",	"BIRC5",	"BRCA2",	"BUB1",	"BUB3",	"KNL1",	"CASP8AP2",	"CBX1",	"CCNA2",	"CCNB2",	"CCND1",	"CCNF",	"CCNT1",	"CDC20",	"CDC25A",	"CDC25B",	"CDC27",	"CDC45",	"CDC6",	"CDC7",	"CDK1",	"CDK4",	"CDKN1B",	"CDKN2C",	"CDKN3",	"CENPA",	"CENPE",	"CENPF",	"CHAF1A",	"CHEK1",	"CHMP1A",	"CKS1B",	"CKS2",	"CTCF",	"CUL1",	"CUL3",	"CUL4A",	"CUL5",	"DBF4",	"DDX39A",	"DKC1",	"DMD",	"DR1",	"DTYMK",	"E2F1",	"E2F2",	"E2F3",	"E2F4",	"EFNA5",	"EGF",	"ESPL1",	"EWSR1",	"EXO1",	"EZH2",	"FANCC",	"FBXO5",	"FOXN3",	"G3BP1",	"GINS2",	"GSPT1",	"H2AZ2",	"H2AX",	"H2AZ1",	"HIF1A",	"HIRA",	"H2BC12",	"HMGA1",	"HMGB3",	"HMGN2",	"HMMR",	"JPT1",	"HNRNPD",	"HNRNPU",	"HOXC10",	"HSPA8",	"HUS1",	"ILF3",	"INCENP",	"KATNA1",	"KIF11",	"KIF15",	"KIF20B",	"KIF22",	"KIF23",	"KIF2C",	"KIF4A",	"KIF5B",	"KPNA2",	"KPNB1",	"LBR",	"LIG3",	"LMNB1",	"MAD2L1",	"MAPK14",	"MARCKS",	"MCM2",	"MCM3",	"MCM5",	"MCM6",	"MEIS1",	"MEIS2",	"MKI67",	"MNAT1",	"MT2A",	"MTF2",	"MYBL2",	"MYC",	"NASP",	"NCL",	"NDC80",	"NEK2",	"NOLC1",	"NOTCH2",	"NUMA1",	"NUP50",	"NUP98",	"NUSAP1",	"ODC1",	"ODF2",	"ORC5",	"ORC6",	"PAFAH1B1",	"TENT4A",	"PBK",	"PDS5B",	"PLK1",	"PLK4",	"PML",	"POLA2",	"POLE",	"POLQ",	"PRC1",	"PRIM2",	"PRMT5",	"PRPF4B",	"PTTG1",	"PTTG3P",	"PURA",	"RACGAP1",	"RAD21",	"RAD23B",	"RAD54L",	"RASAL2",	"RBL1",	"RBM14",	"RPA2",	"RPS6KA5",	"SAP30",	"KMT5A",	"SFPQ",	"SLC12A2",	"SLC38A1",	"SLC7A1",	"SLC7A5",	"SMAD3",	"SMARCC1",	"SMC1A",	"SMC2",	"SMC4",	"SNRPD1",	"SQLE",	"SRSF1",	"SRSF10",	"SRSF2",	"SS18",	"STAG1",	"STIL",	"STMN1",	"SUV39H1",	"SYNCRIP",	"TACC3",	"TFDP1",	"TGFB1",	"TLE3",	"TMPO",	"TNPO2",	"TOP1",	"TOP2A",	"TPX2",	"TRA2B",	"TRAIP",	"TROAP",	"TTK",	"UBE2C",	"UBE2S",	"UCK2",	"UPF1",	"NSD2",	"WRN",	"XPO1",	"YTHDC1",	"MAP3K20")
geneLists <- list(G2MCHKPT=geneList1)

#load motif database of choice
#using hg38 500bp up and 100 bp down rankings, you can use the appropriate feather file from the database on' https://resources.aertslab.org/cistarget/' 

data(motifAnnotations_hgnc)
motifRankings <- importRankings("hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather") 

#to check the annotations format
motifAnnotations[1:5,]

#easiest way to use is cistarget, it combines three steps: calcAUC, addMotifAnnotation, addSignificantGenes
#if the overlap between genelist and ranking database is less than 80%, the function won't work and error will be displayed.
motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                  motifAnnot=motifAnnotations)

resultsSubset <- motifEnrichmentTable_wGenes[1:6,]

#shows the TFs with enriched motifs in the region of 500bp up and 100 bp down of TSS in genes included in the gene list
showLogo(resultsSubset)

#addlogo
motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)

resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:10,]

#to see enriched genes for each Motif

dim(motifEnrichmentTable_wGenes)

geneSetName <- "G2MCHKPT"

#we can select top 3 using head. we can use 'sample' to randomly select 3 from the object

selectedMotifs <- c(head(motifEnrichmentTable_wGenes$motif, 3))

par(mfrow=c(2,2))

#To list the enrich genes associated with motifs selected in the previous command
#can change plotCurve to 'TRUE' to visualize the plots of enriched genes in selected motifs
getSignificantGenes(geneSet=geneLists$G2MCHKPT,, 
                    rankings=motifRankings,
                    signifRankingNames=selectedMotifs,
                    plotCurve=FALSE, maxRank=5000, genesFormat="geneList",
                    method="iCisTarget")


library(DT)

datatable(resultsSubset[,-c("enrichedGenes"), with=FALSE], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))

anotatedTfs <- lapply(split(motifEnrichmentTable_wGenes$TF_highConf,
                            motifEnrichmentTable_wGenes$geneSet),
                      function(x) {
                        genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
                        genesSplit <- unique(unlist(strsplit(genes, "; ")))
                        return(genesSplit)
                      })

anotatedTfs$G2MCHKPT

#building a network, we can choose 3 or any numberof motifs

signifMotifNames <- motifEnrichmentTable_wGenes$motif[1:3]

incidenceMatrix <- getSignificantGenes(geneLists$G2MCHKPT, 
                                       motifRankings,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000, 
                                       genesFormat="incidMatrix",
                                       method="iCisTarget")$incidMatrix

library(reshape2)
edges <- melt(incidenceMatrix)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")

library(visNetwork)
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))
visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE, 
                                        nodesIdSelection = TRUE)




