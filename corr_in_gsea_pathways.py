import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import openpyxl
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from simko_func import simko

protein_data = pd.read_csv('~/icr/simko/data/simko2_data/passport_prots.csv', index_col=0)
exp_data = pd.read_excel('~/icr/simko/data/pbrm1_experiment_data/all_pbrm1_results.xlsx')
exp_data.columns = exp_data.columns.str.replace(' ', '_').str.replace('-', '').str.replace('/', '_')
#performing simko
class_df = simko.get_classes_by_mean_abundance(ko_proteins=['PBRM1'], abundance=protein_data, n=30)
control_diffs = simko.get_control_differentials(protein_data, class_df, k=30)
diffs = simko.get_ko_differentials(abundance=protein_data, class_df=class_df)
diff_stats = simko.get_significance(diffs, control_diffs, n=30)
diff_stats['diff_simko'] = diff_stats['diff']
diff_stats.sort_values(by = 'adjusted_p', ascending=True).head(50)

#pathways/proteins of interest
#asigning proteins to the different 'keywords'
allograft_rejection = ("AARS1","ABCE1","ABI1","ACHE","ACVR2A","AKT1","APBB1","B2M","BCAT1",
                       "BCL10","BCL3","BRCA1","C2","CAPG","CARTPT","CCL11","CCL13","CCL19",
                       "CCL2","CCL22","CCL4","CCL5","CCL7","CCND2","CCND3","CCR1","CCR2",
                       "CCR5","CD1D","CD2","CD247","CD28","CD3D","CD3E","CD3G","CD4","CD40",
                       "CD40LG","CD47","CD7","CD74","CD79A","CD80","CD86","CD8A","CD8B",
                       "CD96","CDKN2A","CFP","CRTAM","CSF1","CSK","CTSS","CXCL13","CXCL9",
                       "CXCR3","DARS1","DEGS1","DYRK3","EGFR","EIF3A","EIF3D","EIF3J",
                       "EIF4G3","EIF5A","ELANE","ELF4","EREG","ETS1","F2","F2R","FAS",
                       "FASLG","FCGR2B","FGR","FLNA","FYB1","GALNT1","GBP2","GCNT1",
                       "GLMN","GPR65","GZMA","GZMB","HCLS1","HDAC9","HIF1A","HLA-A",
                       "HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DQA1","HLA-DRA",
                       "HLA-E","HLA-G","ICAM1","ICOSLG","IFNAR2","IFNG","IFNGR1",
                       "IFNGR2","IGSF6","IKBKB","IL10","IL11","IL12A","IL12B","IL12RB1",
                       "IL13","IL15","IL16","IL18","IL18RAP","IL1B","IL2","IL27RA",
                       "IL2RA","IL2RB","IL2RG","IL4","IL4R","IL6","IL7","IL9","INHBA",
                       "INHBB","IRF4","IRF7","IRF8","ITGAL","ITGB2","ITK","JAK2","KLRD1",
                       "KRT1","LCK","LCP2","LIF","LTB","LY75","LY86","LYN","MAP3K7",
                       "MAP4K1","MBL2","MMP9","MRPL3","MTIF2","NCF4","NCK1","NCR1",
                       "NLRP3","NME1","NOS2","NPM1","PF4","PRF1","PRKCB","PRKCG",
                       "PSMB10","PTPN6","PTPRC","RARS1","RIPK2","RPL39","RPL3L","RPL9",
                       "RPS19","RPS3A","RPS9","SIT1","SOCS1","SOCS5","SPI1","SRGN",
                       "ST8SIA4","STAB1","STAT1","STAT4","TAP1","TAP2","TAPBP","TGFB1",
                       "TGFB2","THY1","TIMP1","TLR1","TLR2","TLR3","TLR6","TNF","TPD52",
                       "TRAF2","TRAT1","UBE2D1","UBE2N","WARS1","WAS","ZAP70")

apoptosis = ("ADD1","AIFM3","ANKH","ANXA1","APP","ATF3","AVPR1A","BAX","BCAP31","BCL10","BCL2L1",
             "BCL2L10","BCL2L11","BCL2L2","BGN","BID","BIK","BIRC3","BMF","BMP2","BNIP3L","BRCA1",
             "BTG2","BTG3","CASP1","CASP2","CASP3","CASP4","CASP6","CASP7","CASP8","CASP9","CAV1",
             "CCNA1","CCND1","CCND2","CD14","CD2","CD38","CD44","CD69","CDC25B","CDK2","CDKN1A",
             "CDKN1B","CFLAR","CLU","CREBBP","CTH","CTNNB1","CYLD","DAP","DAP3","DCN","DDIT3",
             "DFFA","DIABLO","DNAJA1","DNAJC3","DNM1L","DPYD","EBP","EGR3","EMP1","ENO2","ERBB2",
             "ERBB3","EREG","ETF1","F2","F2R","FAS","FASLG","FDXR","FEZ1","GADD45A","GADD45B",
             "GCH1","GNA15","GPX1","GPX3","GPX4","GSN","GSR","GSTM1","GUCY2D","H1-0","HGF",
             "HMGB2","HMOX1","HSPB1","IER3","IFITM3","IFNB1","IFNGR1","IGF2R","IGFBP6","IL18",
             "IL1A","IL1B","IL6","IRF1","ISG20","JUN","KRT18","LEF1","LGALS3","LMNA","LUM",
             "MADD","MCL1","MGMT","MMP2","NEDD9","NEFH","PAK1","PDCD4","PDGFRB","PEA15","PLAT",
             "PLCB2","PLPPR4","PMAIP1","PPP2R5B","PPP3R1","PPT1","PRF1","PSEN1","PSEN2","PTK2",
             "RARA","RELA","RETSAT","RHOB","RHOT2","RNASEL","ROCK1","SAT1","SATB1","SC5D",
             "SLC20A1","SMAD7","SOD1","SOD2","SPTAN1","SQSTM1","TAP1","TGFB2","TGFBR3","TIMP1",
             "TIMP2","TIMP3","TNF","TNFRSF12A","TNFSF10","TOP2A","TSPO","TXNIP","VDAC2","WEE1","XIAP")

bile_acid_metabolism = ("ABCA1","ABCA2","ABCA3","ABCA4","ABCA5","ABCA6","ABCA8","ABCA9","ABCD1",
                        "ABCD2","ABCD3","ABCG4","ABCG8","ACSL1","ACSL5","AGXT","AKR1D1","ALDH1A1",
                        "ALDH8A1","ALDH9A1","AMACR","APOA1","AQP9","AR","ATXN1","BBOX1","BCAR3",
                        "BMP6","CAT","CH25H","CROT","CYP27A1","CYP39A1","CYP46A1","CYP7A1","CYP7B1",
                        "CYP8B1","DHCR24","DIO1","DIO2","EFHC1","EPHX2","FADS1","FADS2","FDXR","GC",
                        "GCLM","GNMT","GNPAT","GSTK1","HACL1","HAO1","HSD17B11","HSD17B4","HSD17B6",
                        "HSD3B1","HSD3B7","IDH1","IDH2","IDI1","ISOC1","KLF1","LCK","LIPE","LONP2",
                        "MLYCD","NEDD4","NPC1","NR0B2","NR1H4","NR1I2","NR3C2","NUDT12","OPTN",
                        "PAOX","PECR","PEX1","PEX11A","PEX11G","PEX12","PEX13","PEX16","PEX19",
                        "PEX26","PEX6","PEX7","PFKM","PHYH","PIPOX","PNPLA8","PRDX5","PXMP2","RBP1",
                        "RETSAT","RXRA","RXRG","SCP2","SERPINA6","SLC22A18","SLC23A1","SLC23A2",
                        "SLC27A2","SLC27A5","SLC29A1","SLC35B2","SLCO1A2","SOAT2","SOD1","SULT1B1",
                        "SULT2B1","TFCP2L1","TTR")

cholesterol_homeostasis = ("ABCA2","ACAT2","ACSS2","ACTG1","ADH4","ALCAM","ALDOC","ANTXR2","ANXA13",
                           "ANXA5","ATF3","ATF5","ATXN2","AVPR1A","CBS","CD9","CHKA","CLU","CPEB2",
                           "CTNNB1","CXCL16","CYP51A1","DHCR7","EBP","ECH1","ERRFI1","ETHE1",
                           "FABP5","FADS2","FASN","FBXO6","FDFT1","FDPS","GLDC","GNAI1","GPX8",
                           "GSTM2","GUSB","HMGCR","HMGCS1","HSD17B7","IDI1","JAG1","LDLR","LGALS3",
                           "LGMN","LPL","LSS","MAL2","MVD","MVK","NFIL3","NIBAN1","NSDHL","PCYT2",
                           "PDK3","PLAUR","PLSCR1","PMVK","PNRC1","PPARG","S100A11","SC5D","SCD",
                           "SEMA3B","SQLE","SREBF2","STARD4","STX5","TM7SF2","TMEM97","TNFRSF12A",
                           "TP53INP1","TRIB3")

coagulation = ("A2M","ACOX2","ADAM9","ANG","ANXA1","APOA1","APOC1","APOC2","APOC3","ARF4","BMP1",
               "C1QA","C1R","C1S","C2","C3","C8A","C8B","C8G","C9","CAPN2","CAPN5","CASP9","CD9",
               "CFB","CFD","CFH","CFI","CLU","COMP","CPB2","CPN1","CPQ","CRIP2","CSRP1","CTSB",
               "CTSE","CTSH","CTSK","CTSO","CTSV","DCT","DPP4","DUSP14","DUSP6","F10","F11",
               "F12","F13B","F2","F2RL2","F3","F8","F9","FBN1","FGA","FGG","FN1","FURIN","FYN",
               "GDA","GNB2","GNG12","GP1BA","GP9","GSN","HMGCS2","HNF4A","HPN","HRG","HTRA1",
               "ISCU","ITGA2","ITGB3","ITIH1","KLF7","KLK8","KLKB1","LAMP2","LEFTY2","LGMN",
               "LRP1","LTA4H","MAFF","MASP2","MBL2","MEP1A","MMP1","MMP10","MMP11","MMP14",
               "MMP15","MMP2","MMP3","MMP7","MMP8","MMP9","MSRB2","MST1","OLR1","P2RY1",
               "PDGFB","PECAM1","PEF1","PF4","PLAT","PLAU","PLEK","PLG","PREP","PROC","PROS1",
               "PROZ","PRSS23","RABIF","RAC1","RAPGEF3","RGN","S100A1","S100A13","SERPINA1",
               "SERPINB2","SERPINC1","SERPINE1","SERPING1","SH2B2","SIRT2","SPARC","TF",
               "TFPI2","THBD","THBS1","TIMP1","TIMP3","TMPRSS6","USP11","VWF","WDR1")

emt = ("ABI3BP","ACTA2","ADAM12","ANPEP","APLP1","AREG","BASP1","BDNF","BGN","BMP1","CADM1",
       "CALD1","CALU","CAP2","CAPG","CCN1","CCN2","CD44","CD59","CDH11","CDH2","CDH6",
       "COL11A1","COL12A1","COL16A1","COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1",
       "COL5A2","COL5A3","COL6A2","COL6A3","COL7A1","COL8A2","COLGALT1","COMP","COPA",
       "CRLF1","CTHRC1","CXCL1","CXCL12","CXCL6","CXCL8","DAB2","DCN","DKK1","DPYSL3",
       "DST","ECM1","ECM2","EDIL3","EFEMP2","ELN","EMP3","ENO2","FAP","FAS","FBLN1",
       "FBLN2","FBLN5","FBN1","FBN2","FERMT2","FGF2","FLNA","FMOD","FN1","FOXC2","FSTL1",
       "FSTL3","FUCA1","FZD8","GADD45A","GADD45B","GAS1","GEM","GJA1","GLIPR1","GPC1",
       "GPX7","GREM1","HTRA1","ID2","IGFBP2","IGFBP3","IGFBP4","IL15","IL32","IL6","INHBA",
       "ITGA2","ITGA5","ITGAV","ITGB1","ITGB3","ITGB5","JUN","LAMA1","LAMA2","LAMA3",
       "LAMC1","LAMC2","LGALS1","LOX","LOXL1","LOXL2","LRP1","LRRC15","LUM","MAGEE1",
       "MATN2","MATN3","MCM7","MEST","MFAP5","MGP","MMP1","MMP14","MMP2","MMP3","MSX1",
       "MXRA5","MYL9","MYLK","NID2","NNMT","NOTCH2","NT5E","NTM","OXTR","P3H1","PCOLCE",
       "PCOLCE2","PDGFRB","PDLIM4","PFN2","PLAUR","PLOD1","PLOD2","PLOD3","PMEPA1","PMP22",
       "POSTN","PPIB","PRRX1","PRSS2","PTHLH","PTX3","PVR","QSOX1","RGS4","RHOB","SAT1",
       "SCG2","SDC1","SDC4","SERPINE1","SERPINE2","SERPINH1","SFRP1","SFRP4","SGCB","SGCD",
       "SGCG","SLC6A8","SLIT2","SLIT3","SNAI2","SNTB1","SPARC","SPOCK1","SPP1","TAGLN",
       "TFPI2","TGFB1","TGFBI","TGFBR3","TGM2","THBS1","THBS2","THY1","TIMP1","TIMP3",
       "TNC","TNFAIP3","TNFRSF11B","TNFRSF12A","TPM1","TPM2","TPM4","VCAM1","VCAN","VEGFA",
       "VEGFC","VIM","WIPF1","WNT5A")

fatty_acid_metabolism = ("AADAT","ACAA1","ACAA2","ACADL","ACADM","ACADS","ACADVL","ACAT2","ACO2",
                         "ACOT2","ACOT8","ACOX1","ACSL1","ACSL4","ACSL5","ACSM3","ACSS1","ADH1C",
                         "ADH7","ADIPOR2","ADSL","ALAD","ALDH1A1","ALDH3A1","ALDH3A2","ALDH9A1",
                         "ALDOA","AOC3","APEX1","AQP7","AUH","BCKDHB","BLVRA","BMPR1B","BPHL",
                         "CA2","CA4","CA6","CBR1","CBR3","CD1D","CD36","CEL","CIDEA","CPOX",
                         "CPT1A","CPT2","CRAT","CRYZ","CYP1A1","CYP4A11","CYP4A22","D2HGDH",
                         "DECR1","DHCR24","DLD","DLST","ECH1","ECHS1","ECI1","ECI2","EHHADH",
                         "ELOVL5","ENO2","ENO3","EPHX1","ERP29","ETFDH","FABP1","FABP2","FASN",
                         "FH","FMO1","G0S2","GABARAPL1","GAD2","GAPDHS","GCDH","GLUL","GPD1",
                         "GPD2","GRHPR","GSTZ1","H2AZ1","HADH","HADHB","HAO2","HCCS","HIBCH",
                         "HMGCL","HMGCS1","HMGCS2","HPGD","HSD17B10","HSD17B11","HSD17B4",
                         "HSD17B7","HSDL2","HSP90AA1","HSPH1","IDH1","IDH3B","IDH3G","IDI1",
                         "IL4I1","INMT","KMT5A","LDHA","LGALS1","LTC4S","MAOA","MCEE","MDH1",
                         "MDH2","ME1","METAP1","MGLL","MIF","MIX23","MLYCD","NBN","NCAPH2","NSDHL",
                         "NTHL1","ODC1","OSTC","PCBD1","PDHA1","PDHB","PPARA","PRDX6","PSME1",
                         "PTPRG","PTS","RAP1GDS1","RDH11","RDH16","REEP6","RETSAT","S100A10","SDHA",
                         "SDHC","SDHD","SERINC1","SLC22A5","SMS","SUCLA2","SUCLG1","SUCLG2","TDO2",
                         "TP53INP2","UBE2L6","UGDH","UROD","UROS","VNN1","XIST","YWHAH")

heme_metabolism = ("ABCB6","ABCG2","ACKR1","ACP5","ACSL6","ADD1","ADD2","ADIPOR1","AGPAT4","AHSP","ALAD",
                   "ALAS2","ALDH1L1","ALDH6A1","ANK1","AQP3","ARHGEF12","ARL2BP","ASNS","ATG4A","ATP6V0A1",
                   "BACH1","BCAM","BLVRA","BLVRB","BMP2K","BNIP3L","BPGM","BSG","BTG2","BTRC","C3","CA1",
                   "CA2","CAST","CAT","CCDC28A","CCND3","CDC27","CDR2","CIR1","CLCN3","CLIC2","CPOX",
                   "CROCCP2","CTNS","CTSB","CTSE","DAAM1","DCAF10","DCAF11","DCUN1D1","DMTN","E2F2",
                   "EIF2AK1","ELL2","ENDOD1","EPB41","EPB42","EPOR","ERMAP","EZH1","FBXO34","FBXO7",
                   "FBXO9","FECH","FN3K","FOXJ2","FOXO3","FTCD","GAPVD1","GATA1","GCLC","GCLM","GDE1",
                   "GLRX5","GMPS","GYPA","GYPB","GYPC","GYPE","H1-0","H4C3","HAGH","HBB","HBBP1","HBD",
                   "HBQ1","HBZ","HDGF","HEBP1","HMBS","HTATIP2","HTRA2","ICAM4","IGSF3","ISCA1","KAT2B",
                   "KDM7A","KEL","KHNYN","KLF1","KLF3","LAMP2","LMO2","LPIN2","LRP10","MAP2K3","MARCHF2",
                   "MARCHF8","MARK3","MBOAT2","MFHAS1","MGST3","MINPP1","MKRN1","MOCOS","MOSPD1","MPP1",
                   "MXI1","MYL4","NARF","NCOA4","NEK7","NFE2","NFE2L1","NNT","NR3C1","NUDT4","OPTN","OSBP2",
                   "P4HA2","PC","PDZK1IP1","PGLS","PICALM","PIGQ","PPOX","PPP2R5B","PRDX2","PSMD9","RAD23A",
                   "RANBP10","RAP1GAP","RBM38","RBM5","RCL1","RHAG","RHCE","RHD","RIOK3","RNF123","RNF19A",
                   "SDCBP","SEC14L1","SELENBP1","SIDT2","SLC10A3","SLC11A2","SLC22A4","SLC25A37","SLC25A38",
                   "SLC2A1","SLC30A1","SLC30A10","SLC4A1","SLC66A2","SLC6A8","SLC6A9","SLC7A11","SMOX",
                   "SNCA","SPTA1","SPTB","SYNJ1","TAL1","TCEA1","TENT5C","TFDP2","TFRC","TMCC2","TMEM9B",
                   "TNRC6B","TNS1","TOP1","TRAK2","TRIM10","TRIM58","TSPAN5","TSPO2","TYR","UBAC1","UCP2",
                   "UROD","UROS","USP15","VEZF1","XK","XPO7","YPEL5")

il2_stat5_signalling = ("ABCB1","ADAM19","AGER","AHCY","AHNAK","AHR","ALCAM","AMACR","ANXA4","APLP1",
                        "ARL4A","BATF","BATF3","BCL2","BCL2L1","BHLHE40","BMP2","BMPR2","CA2","CAPG",
                        "CAPN3","CASP3","CCND2","CCND3","CCNE1","CCR4","CD44","CD48","CD79B","CD81","CD83",
                        "CD86","CDC42SE2","CDC6","CDCP1","CDKN1C","CISH","CKAP4","COCH","COL6A1","CSF1",
                        "CSF2","CST7","CTLA4","CTSZ","CXCL10","CYFIP1","DCPS","DENND5A","DHRS3","DRC1",
                        "ECM1","EEF1AKMT1","EMP1","ENO3","ENPP1","EOMES","ETFBKMT","ETV4","F2RL2","FAH",
                        "FGL2","FLT3LG","FURIN","GABARAPL1","GADD45B","GALM","GATA1","GBP4","GLIPR2",
                        "GPR65","GPR83","GPX4","GSTO1","GUCY1B1","HIPK2","HK2","HOPX","HUWE1","HYCC2",
                        "ICOS","IFITM3","IFNGR1","IGF1R","IGF2R","IKZF2","IKZF4","IL10","IL10RA","IL13",
                        "IL18R1","IL1R2","IL1RL1","IL2RA","IL2RB","IL3RA","IL4R","IRF4","IRF6","IRF8",
                        "ITGA6","ITGAE","ITGAV","ITIH5","KLF6","LCLAT1","LIF","LRIG1","LRRC8C","LTB",
                        "MAFF","MAP3K8","MAP6","MAPKAPK2","MUC1","MXD1","MYC","MYO1C","MYO1E","NCOA3",
                        "NCS1","NDRG1","NFIL3","NFKBIZ","NOP2","NRP1","NT5E","ODC1","P2RX4","P4HA1",
                        "PDCD2L","PENK","PHLDA1","PHTF2","PIM1","PLAGL1","PLEC","PLIN2","PLPP1","PLSCR1",
                        "PNP","POU2F1","PRAF2","PRKCH","PRNP","PTCH1","PTGER2","PTH1R","PTRH2","PUS1",
                        "RABGAP1L","RGS16","RHOB","RHOH","RNH1","RORA","RRAGD","S100A1","SCN9A","SELL",
                        "SELP","SERPINB6","SERPINC1","SH3BGRL2","SHE","SLC1A5","SLC29A2","SLC2A3","SLC39A8",
                        "SMPDL3A","SNX14","SNX9","SOCS1","SOCS2","SPP1","SPRED2","SPRY4","ST3GAL4","SWAP70",
                        "SYNGR2","SYT11","TGM2","TIAM1","TLR7","TNFRSF18","TNFRSF1B","TNFRSF21","TNFRSF4",
                        "TNFRSF8","TNFRSF9","TNFSF10","TNFSF11","TRAF1","TTC39B","TWSG1","UCK2","UMPS",
                        "WLS","XBP1")

interferon_gamma_response = ("ADAR","APOL6","ARID5B","ARL4A","AUTS2","B2M","BANK1","BATF2","BPGM","BST2",
                             "BTG1","C1R","C1S","CASP1","CASP3","CASP4","CASP7","CASP8","CCL2","CCL5",
                             "CCL7","CD274","CD38","CD40","CD69","CD74","CD86","CDKN1A","CFB","CFH","CIITA",
                             "CMKLR1","CMPK2","CMTR1","CSF2RB","CXCL10","CXCL11","CXCL9","DDX60","DHX58",
                             "EIF2AK2","EIF4E3","EPSTI1","FAS","FCGR1A","FGL2","FPR1","GBP4","GBP6","GCH1",
                             "GPR18","GZMA","HELZ2","HERC6","HIF1A","HLA-A","HLA-B","HLA-DMA","HLA-DQA1",
                             "HLA-DRB1","HLA-G","ICAM1","IDO1","IFI27","IFI30","IFI35","IFI44","IFI44L",
                             "IFIH1","IFIT1","IFIT2","IFIT3","IFITM2","IFITM3","IFNAR2","IL10RA","IL15",
                             "IL15RA","IL18BP","IL2RB","IL4R","IL6","IL7","IRF1","IRF2","IRF4","IRF5",
                             "IRF7","IRF8","IRF9","ISG15","ISG20","ISOC1","ITGB7","JAK2","KLRK1","LAP3",
                             "LATS2","LCP2","LGALS3BP","LY6E","LYSMD2","MARCHF1","MT2A","MTHFD2","MVP",
                             "MX1","MX2","MYD88","NAMPT","NCOA3","NFKB1","NFKBIA","NLRC5","NMI","NOD1",
                             "NUP93","OAS2","OAS3","OASL","OGFR","P2RY14","PARP12","PARP14","PDE4B",
                             "PELI1","PFKP","PIM1","PLA2G4A","PLSCR1","PML","PNP","PNPT1","PSMA2","PSMA3",
                             "PSMB10","PSMB2","PSMB8","PSMB9","PSME1","PSME2","PTGS2","PTPN1","PTPN2",
                             "PTPN6","RAPGEF6","RBCK1","RIGI","RIPK1","RIPK2","RNF213","RNF31","RSAD2",
                             "RTP4","SAMD9L","SAMHD1","SECTM1","SELP","SERPING1","SLAMF7","SLC25A28",
                             "SOCS1","SOCS3","SOD2","SP110","SPPL2A","SRI","SSPN","ST3GAL5","ST8SIA4",
                             "STAT1","STAT2","STAT3","STAT4","TAP1","TAPBP","TDRD7","TMT1B","TNFAIP2",
                             "TNFAIP3","TNFAIP6","TNFSF10","TOR1B","TRAFD1","TRIM14","TRIM21","TRIM25",
                             "TRIM26","TXNIP","UBE2L6","UPP1","USP18","VAMP5","VAMP8","VCAM1","WARS1",
                             "XAF1","XCL1","ZBP1","ZNFX1")

myogenesis = ("ABLIM1","ACHE","ACSL1","ACTA1","ACTC1","ACTN2","ACTN3","ADAM12","ADCY9","AEBP1","AGL","AGRN",
              "AK1","AKT2","ANKRD2","APLNR","APOD","APP","ATP2A1","ATP6AP1","BAG1","BDKRB2","BHLHE40",
              "BIN1","CACNA1H","CACNG1","CAMK2B","CASQ1","CASQ2","CAV3","CD36","CDH13","CDKN1A","CFD",
              "CHRNA1","CHRNB1","CHRNG","CKB","CKM","CKMT2","CLU","CNN3","COL15A1","COL1A1","COL3A1",
              "COL4A2","COL6A2","COL6A3","COX6A2","COX7A1","CRAT","CRYAB","CSRP3","CTF1","DAPK2","DENND2B",
              "DES","DMD","DMPK","DTNA","EFS","EIF4A2","ENO3","EPHB3","ERBB3","FABP3","FDPS","FGF2","FHL1",
              "FKBP1B","FLII","FOXO4","FST","FXYD1","GAA","GABARAPL2","GADD45B","GJA5","GNAO1","GPX3",
              "GSN","HBEGF","HDAC5","HRC","HSPB2","HSPB8","IFRD1","IGF1","IGFBP3","IGFBP7","ITGA7","ITGB1",
              "ITGB4","ITGB5","KCNH1","KCNH2","KIFC3","KLF5","LAMA2","LARGE1","LDB3","LPIN1","LSP1",
              "MAPK12","MAPRE3","MB","MEF2A","MEF2C","MEF2D","MRAS","MYBPC3","MYBPH","MYF6","MYH1","MYH11",
              "MYH2","MYH3","MYH4","MYH7","MYH8","MYH9","MYL1","MYL11","MYL2","MYL3","MYL4","MYL6B","MYL7",
              "MYLK","MYO1C","MYOG","MYOM1","MYOM2","MYOZ1","NAV2","NCAM1","NOS1","NOTCH1","NQO1","OCEL1",
              "PC","PDE4DIP","PDLIM7","PFKM","PGAM2","PICK1","PKIA","PLXNB2","PPFIA4","PPP1R3C","PRNP",
              "PSEN2","PTGIS","PTP4A3","PVALB","PYGM","RB1","REEP1","RIT1","RYR1","SCD","SCHIP1","SGCA",
              "SGCD","SGCG","SH2B1","SH3BGR","SIRT2","SLC6A8","SLN","SMTN","SOD3","SORBS1","SORBS3","SPARC",
              "SPDEF","SPEG","SPHK1","SPTAN1","SSPN","STC2","SVIL","SYNGR2","TAGLN","TCAP","TEAD4","TGFB1",
              "TNNC1","TNNC2","TNNI1","TNNI2","TNNT1","TNNT2","TNNT3","TPD52L1","TPM2","TPM3","TSC2","VIPR1",
              "WWTR1")

tnfa_via_nfkb = ("ABCA1","ACKR3","AREG","ATF3","ATP2B1","B4GALT1","B4GALT5","BCL2A1","BCL3","BCL6","BHLHE40",
                 "BIRC2","BIRC3","BMP2","BTG1","BTG2","BTG3","CCL2","CCL20","CCL4","CCL5","CCN1","CCND1",
                 "CCNL1","CCRL2","CD44","CD69","CD80","CD83","CDKN1A","CEBPB","CEBPD","CFLAR","CLCF1","CSF1",
                 "CSF2","CXCL1","CXCL10","CXCL11","CXCL2","CXCL3","CXCL6","DENND5A","DNAJB4","DRAM1","DUSP1",
                 "DUSP2","DUSP4","DUSP5","EDN1","EFNA1","EGR1","EGR2","EGR3","EHD1","EIF1","ETS2","F2RL1",
                 "F3","FJX1","FOS","FOSB","FOSL1","FOSL2","FUT4","G0S2","GADD45A","GADD45B","GCH1","GEM",
                 "GFPT2","GPR183","HBEGF","HES1","ICAM1","ICOSLG","ID2","IER2","IER3","IER5","IFIH1","IFIT2",
                 "IFNGR2","IL12B","IL15RA","IL18","IL1A","IL1B","IL23A","IL6","IL6ST","IL7R","INHBA","IRF1",
                 "IRS2","JAG1","JUN","JUNB","KDM6B","KLF10","KLF2","KLF4","KLF6","KLF9","KYNU","LAMB3","LDLR",
                 "LIF","LITAF","MAFF","MAP2K3","MAP3K8","MARCKS","MCL1","MSC","MXD1","MYC","NAMPT","NFAT5",
                 "NFE2L2","NFIL3","NFKB1","NFKB2","NFKBIA","NFKBIE","NINJ1","NR4A1","NR4A2","NR4A3","OLR1",
                 "PANX1","PDE4B","PDLIM5","PER1","PFKFB3","PHLDA1","PHLDA2","PLAU","PLAUR","PLEK","PLK2",
                 "PLPP3","PMEPA1","PNRC1","PPP1R15A","PTGER4","PTGS2","PTPRE","PTX3","RCAN1","REL","RELA",
                 "RELB","RHOB","RIGI","RIPK2","RNF19B","SAT1","SDC4","SERPINB2","SERPINB8","SERPINE1","SGK1",
                 "SIK1","SLC16A6","SLC2A3","SLC2A6","SMAD3","SNN","SOCS3","SOD2","SPHK1","SPSB1","SQSTM1",
                 "STAT5A","TANK","TAP1","TGIF1","TIPARP","TLR2","TNC","TNF","TNFAIP2","TNFAIP3","TNFAIP6",
                 "TNFAIP8","TNFRSF9","TNFSF9","TNIP1","TNIP2","TRAF1","TRIB1","TRIP10","TSC22D1","TUBB2A",
                 "VEGFA","YRDC","ZBTB10","ZC3H12A","ZFP36")

xenobiotic_metabolism = ("ABCC2","ABCC3","ABCD2","ABHD6","ACO2","ACOX1","ACOX2","ACOX3","ACP1","ACP2","ACSM1",
                         "ADH1C","ADH5","ADH7","AHCY","AKR1C2","AKR1C3","ALAS1","ALDH2","ALDH3A1","ALDH9A1",
                         "ANGPTL3","AOX1","AP4B1","APOE","AQP9","ARG1","ARG2","ARPP19","ASL","ATOH8","ATP2A2",
                         "BCAR1","BCAT1","BLVRB","BPHL","CA2","CASP6","CAT","CBR1","CCL25","CD36","CDA",
                         "CDO1","CES1","CFB","CNDP2","COMT","CROT","CRP","CSAD","CYB5A","CYFIP2","CYP17A1",
                         "CYP1A1","CYP1A2","CYP26A1","CYP27A1","CYP2C18","CYP2E1","CYP2J2","CYP2S1","CYP4F2",
                         "DCXR","DDAH2","DDC","DDT","DHPS","DHRS1","DHRS7","ECH1","ELOVL5","ENPEP","ENTPD5",
                         "EPHA2","EPHX1","ESR1","ETFDH","ETS2","F10","F11","FABP1","FAH","FAS","FBLN1",
                         "FBP1","FETUB","FMO1","FMO3","G6PC1","GABARAPL1","GAD1","GART","GCH1","GCKR","GCLC",
                         "GCNT2","GNMT","GSR","GSS","GSTA3","GSTM4","GSTO1","GSTT2","HACL1","HES6","HGFAC",
                         "HMOX1","HNF4A","HPRT1","HRG","HSD11B1","HSD17B2","ID2","IDH1","IGF1","IGFBP1",
                         "IGFBP4","IL1R1","IRF8","ITIH1","ITIH4","JUP","KARS1","KYNU","LCAT","LEAP2","LONP1",
                         "LPIN2","MAN1A1","MAOA","MARCHF6","MBL2","MCCC2","MPP2","MT2A","MTHFD1","NDRG2","NFS1",
                         "NINJ1","NMT1","NPC1","NQO1","PAPSS2","PC","PDK4","PDLIM5","PEMT","PGD","PGRMC1","PINK1",
                         "PLG","PMM1","POR","PPARD","PROS1","PSMB10","PTGDS","PTGES","PTGES3","PTGR1","PTS","PYCR1",
                         "RAP1GAP","RBP4","REG1A","RETSAT","SAR1B","SERPINA6","SERPINE1","SERTAD1","SHMT2","SLC12A4",
                         "SLC1A5","SLC22A1","SLC35B1","SLC35D1","SLC46A3","SLC6A12","SLC6A6","SMOX","SPINT2","SSR3",
                         "TAT","TDO2","TGFB2","TKFC","TMBIM6","TMEM176B","TMEM97","TNFRSF1A","TPST1","TTPA","TYR",
                         "UGDH","UPB1","UPP1","VNN1","VTN","XDH")

dna_repair = ("AAAS","ADA","ADCY6","ADRM1","AGO4","AK1","AK3","ALYREF","APRT","ARL6IP1","BCAM","BCAP31","BOLA2",
              "BRF2","CANT1","CCNO","CDA","CETN2","CLP1","CMPK2","COX17","CSTF3","DAD1","DCTN4","DDB1","DDB2",
              "DGCR8","DGUOK","DUT","EDF1","EIF1B","ELL","ELOA","ERCC1","ERCC2","ERCC3","ERCC4","ERCC5","ERCC8",
              "FEN1","GMPR2","GPX4","GSDME","GTF2A2","GTF2B","GTF2F1","GTF2H1","GTF2H3","GTF2H5","GTF3C5","GUK1",
              "HCLS1","HPRT1","IMPDH2","ITPA","LIG1","MPC2","MPG","MRPL40","NCBP2","NELFB","NELFCD","NELFE",
              "NFX1","NME1","NME3","NME4","NPR2","NT5C","NT5C3A","NUDT21","NUDT9","PCNA","PDE4B","PDE6G","PNP",
              "POLA1","POLA2","POLB","POLD1","POLD3","POLD4","POLE4","POLH","POLL","POLR1C","POLR1D","POLR1H",
              "POLR2A","POLR2C","POLR2D","POLR2E","POLR2F","POLR2G","POLR2H","POLR2I","POLR2J","POLR2K","POLR3C",
              "POLR3GL","POM121","PRIM1","RAD51","RAD52","RAE1","RALA","RBX1","REV3L","RFC2","RFC3","RFC4","RFC5",
              "RNMT","RPA2","RPA3","RRM2B","SAC3D1","SDCBP","SEC61A1","SF3A3","SMAD5","SNAPC4","SNAPC5","SRSF6",
              "SSRP1","STX3","SUPT4H1","SUPT5H","SURF1","TAF10","TAF12","TAF13","TAF1C","TAF6","TAF9","TARBP2",
              "TK2","TMED2","TP53","TSG101","TYMS","UMPS","UPF3B","USP11","VPS28","VPS37B","VPS37D","XPC","ZNF707",
              "ZWINT")

e2f_targets = ("AK2","ANP32E","ASF1A","ASF1B","ATAD2","AURKA","AURKB","BARD1","BIRC5","BRCA1","BRCA2","BRMS1L",
               "BUB1B","CBX5","CCNB2","CCNE1","CCP110","CDC20","CDC25A","CDC25B","CDCA3","CDCA8","CDK1","CDK4",
               "CDKN1A","CDKN1B","CDKN2A","CDKN2C","CDKN3","CENPE","CENPM","CHEK1","CHEK2","CIT","CKS1B","CKS2",
               "CNOT9","CSE1L","CTCF","CTPS1","DCK","DCLRE1B","DCTPP1","DDX39A","DEK","DEPDC1","DIAPH3","DLGAP5",
               "DNMT1","DONSON","DSCC1","DUT","E2F8","EED","EIF2S1","ESPL1","EXOSC8","EZH2","GINS1","GINS3",
               "GINS4","GSPT1","H2AX","H2AZ1","HELLS","HMGA1","HMGB2","HMGB3","HMMR","HNRNPD","HUS1","ILF3",
               "ING3","IPO7","JPT1","KIF18B","KIF22","KIF2C","KIF4A","KPNA2","LBR","LIG1","LMNB1","LUC7L3",
               "LYAR","MAD2L1","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MELK","MKI67","MLH1","MMS22L","MRE11",
               "MSH2","MTHFD2","MXD3","MYBL2","MYC","NAA38","NAP1L1","NASP","NBN","NCAPD2","NME1","NOLC1","NOP56",
               "NUDT21","NUP107","NUP153","NUP205","ORC2","ORC6","PA2G4","PAICS","PAN2","PCNA","PDS5B","PHF5A",
               "PLK1","PLK4","PMS2","PNN","POLA2","POLD1","POLD2","POLD3","POLE","POLE4","POP7","PPM1D","PPP1R8",
               "PRDX4","PRIM2","PRKDC","PRPS1","PSIP1","PSMC3IP","PTTG1","RACGAP1","RAD1","RAD21","RAD50",
               "RAD51AP1","RAD51C","RAN","RANBP1","RBBP7","RFC1","RFC2","RFC3","RNASEH2A","RPA1","RPA2","RPA3",
               "RRM2","SHMT1","SLBP","SMC1A","SMC3","SMC4","SMC6","SNRPB","SPAG5","SPC24","SPC25","SRSF1","SRSF2",
               "SSRP1","STAG1","STMN1","SUV39H1","SYNCRIP","TACC3","TBRG4","TCF19","TFRC","TIMELESS","TIPIN","TK1",
               "TMPO","TOP2A","TP53","TRA2B","TRIP13","TUBB","TUBG1","UBE2S","UBE2T","UBR7","UNG","USP1","WDR90",
               "WEE1","XPO1","XRCC6","ZW10")

g2m_checkpoint = ("ABL1","AMD1","ARID4A","ATF5","ATRX","AURKA","AURKB","BARD1","BCL3","BIRC5","BRCA2","BUB1",
                  "BUB3","CASP8AP2","CBX1","CCNA2","CCNB2","CCND1","CCNF","CCNT1","CDC20","CDC25A","CDC25B",
                  "CDC27","CDC45","CDC6","CDC7","CDK1","CDK4","CDKN1B","CDKN2C","CDKN3","CENPA","CENPE","CENPF",
                  "CHAF1A","CHEK1","CHMP1A","CKS1B","CKS2","CTCF","CUL1","CUL3","CUL4A","CUL5","DBF4","DDX39A",
                  "DKC1","DMD","DR1","DTYMK","E2F1","E2F2","E2F3","E2F4","EFNA5","EGF","ESPL1","EWSR1","EXO1",
                  "EZH2","FANCC","FBXO5","FOXN3","G3BP1","GINS2","GSPT1","H2AX","H2AZ1","H2AZ2","H2BC12","HIF1A",
                  "HIRA","HMGA1","HMGB3","HMGN2","HMMR","HNRNPD","HNRNPU","HOXC10","HSPA8","HUS1","ILF3","INCENP",
                  "JPT1","KATNA1","KIF11","KIF15","KIF20B","KIF22","KIF23","KIF2C","KIF4A","KIF5B","KMT5A","KNL1",
                  "KPNA2","KPNB1","LBR","LIG3","LMNB1","MAD2L1","MAP3K20","MAPK14","MARCKS","MCM2","MCM3","MCM5",
                  "MCM6","MEIS1","MEIS2","MKI67","MNAT1","MT2A","MTF2","MYBL2","MYC","NASP","NCL","NDC80","NEK2",
                  "NOLC1","NOTCH2","NSD2","NUMA1","NUP50","NUP98","NUSAP1","ODC1","ODF2","ORC5","ORC6","PAFAH1B1",
                  "PBK","PDS5B","PLK1","PLK4","PML","POLA2","POLE","POLQ","PRC1","PRIM2","PRMT5","PRPF4B","PTTG1",
                  "PTTG3P","PURA","RACGAP1","RAD21","RAD23B","RAD54L","RASAL2","RBL1","RBM14","RPA2","RPS6KA5",
                  "SAP30","SFPQ","SLC12A2","SLC38A1","SLC7A1","SLC7A5","SMAD3","SMARCC1","SMC1A","SMC2","SMC4",
                  "SNRPD1","SQLE","SRSF1","SRSF10","SRSF2","SS18","STAG1","STIL","STMN1","SUV39H1","SYNCRIP",
                  "TACC3","TENT4A","TFDP1","TGFB1","TLE3","TMPO","TNPO2","TOP1","TOP2A","TPX2","TRA2B","TRAIP",
                  "TROAP","TTK","UBE2C","UBE2S","UCK2","UPF1","WRN","XPO1","YTHDC1")

myc_targets = ("ABCE1","ACP1","AIMP2","AP3S1","APEX1","BUB3","C1QBP","CAD","CANX","CBX3","CCNA2","CCT2","CCT3",
               "CCT4","CCT5","CCT7","CDC20","CDC45","CDK2","CDK4","CLNS1A","CNBP","COPS5","COX5A","CSTF2","CTPS1",
               "CUL1","CYC1","DDX18","DDX21","DEK","DHX15","DUT","EEF1B2","EIF1AX","EIF2S1","EIF2S2","EIF3B",
               "EIF3D","EIF3J","EIF4A1","EIF4E","EIF4G2","EIF4H","EPRS1","ERH","ETF1","EXOSC7","FAM120A","FBL",
               "G3BP1","GLO1","GNL3","GOT2","GSPT1","H2AZ1","HDAC2","HDDC2","HDGF","HNRNPA1","HNRNPA2B1","HNRNPA3",
               "HNRNPC","HNRNPD","HNRNPR","HNRNPU","HPRT1","HSP90AB1","HSPD1","HSPE1","IARS1","IFRD1","ILF2",
               "IMPDH2","KARS1","KPNA2","KPNB1","LDHA","LSM2","LSM7","MAD2L1","MCM2","MCM4","MCM5","MCM6","MCM7",
               "MRPL23","MRPL9","MRPS18B","MYC","NAP1L1","NCBP1","NCBP2","NDUFAB1","NHP2","NME1","NOLC1","NOP16",
               "NOP56","NPM1","ODC1","ORC2","PA2G4","PABPC1","PABPC4","PCBP1","PCNA","PGK1","PHB1","PHB2","POLD2",
               "POLE3","PPIA","PPM1G","PRDX3","PRDX4","PRPF31","PRPS2","PSMA1","PSMA2","PSMA4","PSMA6","PSMA7"
               ,"PSMB2","PSMB3","PSMC4","PSMC6","PSMD1","PSMD14","PSMD3","PSMD7","PSMD8","PTGES3","PWP1","RACK1",
               "RAD23B","RAN","RANBP1","RFC4","RNPS1","RPL14","RPL18","RPL22","RPL34","RPL6","RPLP0","RPS10",
               "RPS2","RPS3","RPS5","RPS6","RRM1","RRP9","RSL1D1","RUVBL2","SERBP1","SET","SF3A1","SF3B3",
               "SLC25A3","SMARCC1","SNRPA","SNRPA1","SNRPB2","SNRPD1","SNRPD2","SNRPD3","SNRPG","SRM","SRPK1",
               "SRSF1","SRSF2","SRSF3","SRSF7","SSB","SSBP1","STARD7","SYNCRIP","TARDBP","TCP1","TFDP1","TOMM70",
               "TRA2B","TRIM28","TUFM","TXNL4A","TYMS","U2AF1","UBA2","UBE2E1","UBE2L3","USP1","VBP1","VDAC1",
               "VDAC3","XPO1","XPOT","XRCC6","YWHAE","YWHAQ")

#for loop for getting all results
pathways = {
    "allograft_rejection": allograft_rejection,
    "apoptosis": apoptosis,
    "bile_acid_metabolism": bile_acid_metabolism,
    "cholesterol_homeostasis": cholesterol_homeostasis,
    "coagulation": coagulation,
    "dna_repair": dna_repair,
    "e2f_targets": e2f_targets,
    "emt":emt,
    "fatty_acid_metabolism": fatty_acid_metabolism,
    "g2m_checkpoint": g2m_checkpoint,
    "heme_metabolism": heme_metabolism,
    "il2_stat5_signalling": il2_stat5_signalling,
    "interferon_gamma_response": interferon_gamma_response,
    "myc_targets": myc_targets,
    "myogenesis": myogenesis,
    "tnfa_via_nfkb": tnfa_via_nfkb,
    "xenobiotic_metabolism": xenobiotic_metabolism
}

#creating a data frame that direvtly compares like proteins - only using significant simko proteins 
prots_sf = diff_stats.loc[diff_stats['adjusted_p'] < 0.05]
prots_sf.sort_values(by='adjusted_p', ascending=True)

#only keeping relevant columns to make merge easier
prots_sf = prots_sf[['protein', 'diff_simko']]

exp_data = exp_data.rename(columns={'Gene_Names_(primary)': 'protein', 'mean_log2_KO_WT': 'diff_exp'})
exp_data = exp_data[['protein', 'diff_exp']]

simko_exp_df = pd.merge(prots_sf, exp_data, on='protein')

for pathway, protein_list in pathways.items():
    comparison_df = simko_exp_df.loc[simko_exp_df['protein'].isin(protein_list)]
    plt.figure(figsize=(8, 8))
    sns.regplot(
        x = comparison_df['diff_simko'],
        y = comparison_df['diff_exp'],
        scatter = True,
        ci = None,
        line_kws={'color': 'black', 'linewidth': 1}
    )
    plt.grid(True)
    plt.xlabel('LogFC (SimKO)')
    plt.ylabel('LogFC (Experimental)')
    plt.axhline(0, color='gray', linestyle='--', linewidth=2)
    plt.axvline(0, color='gray', linestyle='--', linewidth=2)
    plt.title(f"Pathway: {pathway}")
    plt.tight_layout()
    plt.show()

#now looking at correlation and significance for these proteins
for pathway, protein_list in pathways.items():
    comparison_df2 = simko_exp_df.loc[simko_exp_df['protein'].isin(protein_list)]
    corr, pvalue = pearsonr(comparison_df2['diff_simko'], comparison_df2['diff_exp'])
    print(f"Pathway: {pathway}, R Value: {corr}, P Value: {pvalue}")




#quickly looking at cvae results
cvae = pd.read_csv('~/icr/simko/autoencoder/continuous_cvae2_results/protein_shift_summary_with_pvalue.csv')
cvae['diff_cvae'] = cvae['diff']
cvae = cvae[['protein', 'diff_cvae', 'p_value']]
cvae = cvae.loc[cvae['p_value'] < 0.05]
cvae = cvae[['protein', 'diff_cvae']]
cvae

cvae_exp_df = pd.merge(cvae, exp_data, on='protein')

#for loop for getting correlations across pathways
for pathway, protein_list in pathways.items():
    comparison_df_cvae = cvae_exp_df.loc[cvae_exp_df['protein'].isin(protein_list)]
    plt.figure(figsize=(8, 8))
    sns.regplot(
        x = comparison_df_cvae['diff_cvae'],
        y = comparison_df_cvae['diff_exp'],
        scatter = True,
        ci = None,
        line_kws={'color': 'black', 'linewidth': 1}
    )
    plt.grid(True)
    plt.xlabel('LogFC (CVAE simulated)')
    plt.ylabel('LogFC (Experimental)')
    plt.axhline(0, color='gray', linestyle='--', linewidth=2)
    plt.axvline(0, color='gray', linestyle='--', linewidth=2)
    plt.title(f"Pathway: {pathway}")
    plt.tight_layout()
    plt.show()

#now looking at correlation and significance for these proteins
for pathway, protein_list in pathways.items():
    comparison_df_cvae2 = cvae_exp_df.loc[cvae_exp_df['protein'].isin(protein_list)]
    corr, pvalue = pearsonr(comparison_df_cvae2['diff_cvae'], comparison_df_cvae2['diff_exp'])
    print(f"Pathway: {pathway}, R Value: {corr}, P Value: {pvalue}")




#pathyway boxplot on the significant cvae results 
cvae['diff_copy'] = cvae['diff_cvae']

cvae_pathways = simko.get_pathway_assigned(cvae)

simko.boxplots_for_pathways(cvae_pathways, x='pathways', y='diff_copy',
                      title='CVAE PBRM1 KO results)',
                      figsize=(10,8),
                      xlabel=' ',
                      ylabel=' ')