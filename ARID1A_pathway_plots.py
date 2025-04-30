import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import openpyxl
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy import stats
from simko_func import simko

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

oxidative_phosphorylation = ("ABCB7","ACAA1","ACAA2","ACADM","ACADSB","ACADVL","ACAT1","ACO2","AFG3L2","AIFM1","ALAS1",
                             "ALDH6A1","ATP1B1","ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E","ATP5MC1","ATP5MC2",
                             "ATP5MC3","ATP5ME","ATP5MF","ATP5MG","ATP5PB","ATP5PD","ATP5PF","ATP5PO","ATP6AP1","ATP6V0B",
                             "ATP6V0C","ATP6V0E1","ATP6V1C1","ATP6V1D","ATP6V1E1","ATP6V1F","ATP6V1G1","ATP6V1H","BAX","BCKDHA",
                             "BDH2","CASP7","COX10","COX11","COX15","COX17","COX4I1","COX5A","COX5B","COX6A1","COX6B1","COX6C",
                             "COX7A2","COX7A2L","COX7B","COX7C","COX8A","CPT1A","CS","CYB5A","CYB5R3","CYC1","CYCS","DECR1","DLAT",
                             "DLD","DLST","ECH1","ECHS1","ECI1","ETFA","ETFB","ETFDH","FDX1","FH","FXN","GLUD1","GOT2","GPI","GPX4",
                             "GRPEL1","HADHA","HADHB","HCCS","HSD17B10","HSPA9","HTRA2","IDH1","IDH2","IDH3A","IDH3B","IDH3G","IMMT",
                             "ISCA1","ISCU","LDHA","LDHB","LRPPRC","MAOB","MDH1","MDH2","MFN2","MGST3","MPC1","MRPL11","MRPL15","MRPL34",
                             "MRPL35","MRPS11","MRPS12","MRPS15","MRPS22","MRPS30","MTRF1","MTRR","MTX2","NDUFA1","NDUFA2","NDUFA3",
                             "NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFAB1","NDUFB1","NDUFB2","NDUFB3","NDUFB4",
                             "NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFC1","NDUFC2","NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS6","NDUFS7",
                             "NDUFS8","NDUFV1","NDUFV2","NNT","NQO2","OAT","OGDH","OPA1","OXA1L","PDHA1","PDHB","PDHX","PDK4","PDP1","PHB2",
                             "PHYH","PMPCA","POLR2F","POR","PRDX3","RETSAT","RHOT1","RHOT2","SDHA","SDHB","SDHC","SDHD","SLC25A11",
                             "SLC25A12","SLC25A20","SLC25A3","SLC25A4","SLC25A5","SLC25A6","SUCLA2","SUCLG1","SUPV3L1","SURF1","TCIRG1",
                             "TIMM10","TIMM13","TIMM17A","TIMM50","TIMM8B","TIMM9","TOMM22","TOMM70","UQCR10","UQCR11","UQCRB","UQCRC1","UQCRC2",
                             "UQCRFS1","UQCRH","UQCRQ","VDAC1","VDAC2","VDAC3")

estrogen_response_early = ("ABAT","ABCA3","ABHD2","ABLIM1","ADCY1","ADCY9","ADD3","AFF1","AKAP1","ALDH3B1","AMFR","ANXA9","AQP3","AR","AREG","ARL3","ASB13","B4GALT1","BAG1","BCL11B",
                           "BCL2","BHLHE40","BLVRB","CA12","CALB2","CALCR","CANT1","CBFA2T3","CCN5","CCND1","CD44","CELSR1","CELSR2","CHPT1",
                           "CISH","CLDN7","CLIC3","CXCL12","CYP26B1","DEPTOR","DHCR7","DHRS2","DHRS3","DLC1","DYNLT3","EEIG1","EGR3","ELF1","ELF3",
                           "ELOVL2","ELOVL5","ENDOD1","ESRP2","FARP1","FASN","FCMR","FDFT1","FHL2","FKBP4","FKBP5","FLNB","FOS","FOXC1","FRK",
                           "GAB2","GFRA1","GJA1","GLA","GREB1","HES1","HR","HSPB8","IGF1R","IGFBP4","IL17RB","IL6ST","INHBB","INPP5F","ISG20L2",
                           "ITPK1","JAK2","KAZN","KCNK15","KCNK5","KDM4B","KLF10","KLF4","KLK10","KRT13","KRT15","KRT18","KRT19","KRT8","LAD1",
                           "LRIG1","MAPT","MAST4","MED13L","MED24","MICB","MINDY1","MLPH","MPPED2","MREG","MSMB","MUC1","MYB","MYBBP1A","MYBL1",
                           "MYC","MYOF","NADSYN1","NAV2","NBL1","NCOR2","NHERF1","NPY1R","NRIP1","NXT1","OLFM1","OLFML3","OPN3","OVOL2","P2RY2",
                           "PAPSS2","PDLIM3","PDZK1","PEX11A","PGR","PLAAT3","PMAIP1","PODXL","PPIF","PRSS23","PTGES","RAB17","RAB31","RAPGEFL1",
                           "RARA","RASGRP1","RBBP8","REEP1","RET","RETREG1","RHOBTB3","RHOD","RPS6KA2","RRP12","SCARB1","SCNN1A","SEC14L2","SEMA3B",
                           "SFN","SH3BP5","SIAH2","SLC16A1","SLC19A2","SLC1A1","SLC1A4","SLC22A5","SLC24A3","SLC26A2","SLC27A2","SLC2A1","SLC37A1",
                           "SLC39A6","SLC7A2","SLC7A5","SNX24","SOX3","STC2","SULT2B1","SVIL","SYBU","SYNGR1","SYT12","TBC1D30","TFAP2C","TFF1",
                           "TFF3","TGIF2","TGM2","THSD4","TIAM1","TIPARP","TJP3","TMEM164","TMPRSS3","TOB1","TPBG","TPD52L1","TSKU","TTC39A",
                           "TUBB2B","UGCG","UNC119","WFS1","WWC1","XBP1","ZNF185")

estrogen_response_late = ("ABCA3","ABHD2","ACOX2","ADD3","AFF1","AGR2","ALDH3A2","ALDH3B1","AMFR","ANXA9","AREG","ARL3","ASCL1",
                          "ASS1","ATP2B4","BAG1","BATF","BCL2","BLVRB","BTG3","CA12","CA2","CACNA2D2","CALCR","CAV1",
                          "CCN5","CCNA1","CCND1","CD44","CD9","CDC20","CDC6","CDH1","CELSR2","CHPT1","CHST8","CISH","CKB",
                          "CLIC3","COX6C","CPE","CXCL12","CXCL14","CYP26B1","CYP4F11","DCXR","DHCR7","DHRS2","DLG5",
                          "DNAJC1","DNAJC12","DUSP2","DYNLT3","EEIG1","EGR3","ELOVL5","EMP2","ETFB","FABP5","FARP1",
                          "FDFT1","FGFR3","FKBP4","FKBP5","FLNB","FOS","FOXC1","FRK","GAL","GALE","GFUS","GINS2","GJB3",
                          "GLA","GPER1","HMGCS2","HOMER2","HPRT1","HR","HSPA4L","HSPB8","ID2","IDH2","IGFBP4","IGSF1",
                          "IL17RB","IL6ST","IMPA2","ISG20","ITPK1","JAK1","JAK2","KCNK5","KIF20A","KLF4","KLK10","KLK11",
                          "KRT13","KRT19","LAMC2","LARGE1","LLGL2","LSR","LTF","MAPK13","MAPT","MDK","MEST","METTL3",
                          "MICB","MOCS2","MYB","MYOF","NAB2","NBL1","NCOR2","NHERF1","NMU","NPY1R","NRIP1","NXT1","OLFM1",
                          "OPN3","OVOL2","PAPSS2","PCP4","PDCD4","PDLIM3","PDZK1","PERP","PGR","PKP3","PLAAT3","PLAC1",
                          "PLK4","PLXNB1","PPIF","PRKAR2B","PRLR","PRSS23","PTGER3","PTGES","PTPN6","RAB31","RABEP1",
                          "RAPGEFL1","RBBP8","RET","RNASEH2A","RPS6KA2","S100A9","SCARB1","SCNN1A","SCUBE2","SEMA3B",
                          "SERPINA1","SERPINA3","SERPINA5","SFN","SGK1","SIAH2","SLC16A1","SLC1A4","SLC22A5","SLC24A3",
                          "SLC26A2","SLC27A2","SLC29A1","SLC2A8","SLC7A5","SNX10","SORD","SOX3","ST14","ST6GALNAC2",
                          "STIL","SULT2B1","TFAP2C","TFF1","TFF3","TFPI2","TH","TIAM1","TJP3","TMPRSS3","TNNC1","TOB1",
                          "TOP2A","TPBG","TPD52L1","TPSAB1","TRIM29","TSPAN13","TST","UGDH","UNC13B","WFS1","XBP1","XRCC3",
                          "ZFP36")

kras_signaling_down = ("ABCB11","ABCG4","ACTC1","ADRA2C","AKR1B10","ALOX12B","AMBN","ARHGDIG","ARPP21","ASB7","ATP4A",
                       "ATP6V1B1","BARD1","BMPR1B","BRDT","BTG2","C5","CACNA1F","CACNG1","CALCB","CALML5","CAMK1D","CAPN9",
                       "CCDC106","CCNA1","CCR8","CD207","CD40LG","CD80","CDH16","CDKAL1","CELSR2","CHRNG","CHST2","CKM","CLDN16",
                       "CLDN8","CLPS","CLSTN3","CNTFR","COL2A1","COPZ2","COQ8A","CPA2","CPB1","CPEB3","CYP11B2","CYP39A1","DCC",
                       "DLK2","DTNB","EDAR","EDN1","EDN2","EFHD1","EGF","ENTPD7","EPHA5","FGF16","FGF22","FGFR3","FGGY","FSHB",
                       "GAMT","GDNF","GP1BA","GP2","GPR19","GPR3","GPRC5C","GRID2","GTF3C5","HNF1A","HSD11B2","HTR1B","HTR1D",
                       "IDUA","IFI44L","IFNG","IGFBP2","IL12B","IL5","INSL5","IRS4","ITGB1BP2","ITIH3","KCND1","KCNE2","KCNMB1",
                       "KCNN1","KCNQ2","KLHDC8A","KLK7","KLK8","KMT2D","KRT1","KRT13","KRT15","KRT4","KRT5","LFNG","LGALS7",
                       "LYPD3","MACROH2A2","MAGIX","MAST3","MEFV","MFSD6","MSH5","MTHFR","MX1","MYH7","MYO15A","MYOT","NGB","NOS1",
                       "NPHS1","NPY4R","NR4A2","NR6A1","NRIP2","NTF3","NUDT11","OXT","P2RX6","P2RY4","PAX3","PAX4","PCDHB1","PDCD1",
                       "PDE6B","PDK2","PKP1","PLAG1","PNMT","PRKN","PRODH","PROP1","PTGFR","PTPRJ","RGS11","RIBC2","RSAD2","RYR1",
                       "RYR2","SCGB1A1","SCN10A","SELENOP","SERPINA10","SERPINB2","SGK1","SHOX2","SIDT1","SKIL","SLC12A3","SLC16A7",
                       "SLC25A23","SLC29A3","SLC30A3","SLC38A3","SLC5A5","SLC6A14","SLC6A3","SMPX","SNCB","SNN","SOX10","SPHK2",
                       "SPRR3","SPTBN2","SSTR4","STAG3","SYNPO","TAS2R4","TCF7L1","TCL1A","TENM2","TENT5C","TEX15","TFAP2B",
                       "TFCP2L1","TFF2","TG","TGFB2","TGM1","THNSL2","THRB","TLX1","TNNI3","TSHB","UGT2B17","UPK3B","VPREB1",
                       "VPS50","WNT16","YBX2","YPEL1","ZBTB16","ZC2HC1C","ZNF112")

adipogenesis = ("ABCA1","ABCB8","ACAA2","ACADL","ACADM","ACADS","ACLY","ACO2","ACOX1","ADCY6","ADIG","ADIPOQ","ADIPOR2",
                "AGPAT3","AIFM1","AK2","ALDH2","ALDOA","ANGPT1","ANGPTL4","APLP2","APOE","ARAF","ARL4A","ATL2","ATP1B3",
                "ATP5PO","BAZ2A","BCKDHA","BCL2L13","BCL6","C3","CAT","CAVIN1","CAVIN2","CCNG2","CD151","CD302","CD36",
                "CDKN2C","CHCHD10","CHUK","CIDEA","CMBL","CMPK1","COL15A1","COL4A1","COQ3","COQ5","COQ9","COX6A1","COX7B",
                "COX8A","CPT2","CRAT","CS","CYC1","CYP4B1","DBT","DDT","DECR1","DGAT1","DHCR7","DHRS7","DHRS7B","DLAT",
                "DLD","DNAJB9","DNAJC15","DRAM2","ECH1","ECHS1","ELMOD3","ELOVL6","ENPP2","EPHX2","ESRRA","ESYT1","ETFB",
                "FABP4","FAH","FZD4","G3BP2","GADD45A","GBE1","GHITM","GPAM","GPAT4","GPD2","GPHN","GPX3","GPX4","GRPEL1",
                "HADH","HIBCH","HSPB8","IDH1","IDH3A","IDH3G","IFNGR1","IMMT","ITGA7","ITIH5","ITSN1","JAGN1","LAMA4",
                "LEP","LIFR","LIPE","LPCAT3","LPL","LTC4S","MAP4K3","MCCC1","MDH2","ME1","MGLL","MGST3","MIGA2","MRAP",
                "MRPL15","MTARC2","MTCH2","MYLK","NABP1","NDUFA5","NDUFAB1","NDUFB7","NDUFS3","NKIRAS1","NMT1","OMD",
                "ORM1","PDCD4","PEMT","PEX14","PFKFB3","PFKL","PGM1","PHLDB1","PHYH","PIM3","PLIN2","POR","PPARG","PPM1B",
                "PPP1R15B","PRDX3","PREB","PTCD3","PTGER3","QDPR","RAB34","REEP5","REEP6","RETN","RETSAT","RIOK3","RMDN3",
                "RNF11","RREB1","RTN3","SAMM50","SCARB1","SCP2","SDHB","SDHC","SLC19A1","SLC1A5","SLC25A1","SLC25A10",
                "SLC27A1","SLC5A6","SLC66A3","SNCG","SOD1","SORBS1","SOWAHC","SPARCL1","SQOR","SSPN","STAT5A","STOM",
                "SUCLG1","SULT1A1","TALDO1","TANK","TKT","TOB1","TST","UBC","UBQLN1","UCK1","UCP2","UQCR10","UQCR11",
                "UQCRC1","UQCRQ","VEGFB","YWHAG")

#loading abundance data
abundance = pd.read_csv('~/icr/simko/data/simko2_data/passport_prots.csv', index_col=0)
abundance.index = abundance.index.astype(str)

#doing simko on the abundance data
#note: the p values will slightly change as the control shuffles are always different 
class_df = simko.get_classes_by_mean_abundance(ko_proteins=['ARID1A'], abundance=abundance, n=30)
control_diffs = simko.get_control_differentials(abundance, class_df, k=30)
diffs = simko.get_ko_differentials(abundance=abundance, class_df=class_df)
diffs
PBRM1_row = diffs[diffs["protein"] == "PBRM1"]
PBRM1_row
diff_stats = simko.get_significance(diffs, control_diffs, n=30)
diff_stats
#-np.log10() on pvals for volcano plot
diff_stats['log10_pval']= -np.log10(diff_stats['adjusted_p'])
diff_stats = diff_stats.sort_values(by='log10_pval', ascending = False)

#in this current data frame, the diff coloumn acts more like a function, rathwr than just a simple value
#for ease im making another coloumn (diff_copy) --> so its just a value not a 'fucntion'
diff_stats['diff_copy'] = diff_stats['diff']
diff_stats[:10]

# now i want to do a group by function?
#divide up data frame based on what pathways the proteins are a part of
diff_stats['pathways'] = diff_stats['protein'].apply(lambda protein: [
    pathway
    for pathway, proteins in {
        "apoptosis": apoptosis,
        "coagulation": coagulation,
        "emt": emt,
        "interferon gamma response": interferon_gamma_response,
        "myogenesis": myogenesis,
        "e2f targets": e2f_targets,
        "myc targets": myc_targets,
        "oxidative phosphorylation": oxidative_phosphorylation,
        "estrogen response early": estrogen_response_early,
        "estrogen response late": estrogen_response_late,
        "KRAS signaling (down)": kras_signaling_down,
        "adipogenesis": adipogenesis
    }.items()
    if protein in proteins
])

#duplicating/tripling etc protein rows that appear in multiple pathways
diff_stats = diff_stats.explode("pathways").reset_index(drop=True)
#getting rid of nan rows
diff_stats_arid1a = diff_stats.dropna(subset=["pathways"]).reset_index(drop=True)
diff_stats_arid1a

#only significant proteins for boxplot
diff_stats_sf = diff_stats_arid1a.loc[diff_stats['adjusted_p'] < 0.0005]
diff_stats_sf

#getting the colors for above of below median  - know what they are already
diff_stats_sf = diff_stats_sf.sort_values('pathways')

blue_pathways = ["coagulation", "e2f targets", "myc targets", "interferon gamma response"]
diff_stats_sf["Colour"] = diff_stats_sf["pathways"].apply(
    lambda pathway: "Mean Downregulation" if pathway in blue_pathways else "Mean Upregulation"
)

plt.figure(figsize=(15,10))
sns.boxplot(data=diff_stats_sf, x='pathways', y='diff_copy',
            palette={"Mean Downregulation": "skyblue", "Mean Upregulation": "red"},
            hue="Colour",
            dodge=False)

# Add custom legend
plt.title("Boxplot of Protein Abundance Changes by Pathway: Pan Cancer", fontsize=20)
plt.xlabel(" ", fontsize=12)
plt.ylabel("Abudnance Change (LogFC)", fontsize=16)
plt.xticks(rotation=30, ha="right", fontsize=16)
plt.axhline(0, color='black', linestyle='--', linewidth=2)
plt.tight_layout()
plt.show()





#now doing the same again but only with breast cancer cell lines

model_lists = pd.read_csv('~/icr/simko/data/simko2_data/model_list_20240110.csv', index_col=0)

#need to make sure that the cell lines in model_list matched us with cell lines in abundance data
abundance_cell_lines = set(abundance.columns)
model_list_filt = model_lists[model_lists['model_name'].isin(abundance_cell_lines)]
unique_tissues = model_list_filt['tissue'].unique()
unique_tissues
#column names to sue: model_name and tissues
# making a dictionary
grouped = model_list_filt.groupby('tissue')['model_name'].apply(list)
tissue_to_cell_lines = grouped.to_dict()
print(tissue_to_cell_lines)
#finding out how many cell lines for each tissue
for key, value in tissue_to_cell_lines.items():
    if isinstance(value, list):  # Check if the value is a list
        print(f"The length of the list for {key} is: {len(value)}")


specific_CL = tissue_to_cell_lines["Ovary"]
abundance_specific_CL = abundance[specific_CL]
abundance_specific_CL

#now we can do diff_stats on the specific cell lines
class_df_CL = simko.get_classes_by_mean_abundance(ko_proteins=['ARID1A'], abundance=abundance_specific_CL, n=15)
control_diffs_CL = simko.get_control_differentials(abundance_specific_CL, class_df_CL, k=15)
diffs_CL = simko.get_ko_differentials(abundance=abundance_specific_CL, class_df=class_df_CL)
diffs_CL
PBRM1_row_CL = diffs_CL[diffs_CL["protein"] == "PBRM1"]
PBRM1_row_CL
diff_stats_CL = simko.get_significance(diffs_CL, control_diffs_CL, n=15)
#-np.log10() on pvals for volcano plot
diff_stats_CL['log10_pval']= -np.log10(diff_stats_CL['adjusted_p'])
diff_stats_CL['diff_copy'] = diff_stats_CL['diff']
diff_stats_CL = diff_stats_CL.sort_values(by='log10_pval', ascending = False)
diff_stats_CL

# now i want to do a group by function?
#divide up data frame based on what pathways the proteins are a part of
diff_stats_CL['pathways'] = diff_stats_CL['protein'].apply(lambda protein: [
    pathway
    for pathway, proteins in {
        "apoptosis": apoptosis,
        "coagulation": coagulation,
        "emt": emt,
        "interferon gamma response": interferon_gamma_response,
        "myogenesis": myogenesis,
        "e2f targets": e2f_targets,
        "myc targets": myc_targets,
        "oxidative phosphorylation": oxidative_phosphorylation,
        "estrogen response early": estrogen_response_early,
        "estrogen response late": estrogen_response_late,
        "KRAS signaling (down)": kras_signaling_down,
        "adipogenesis": adipogenesis
    }.items()
    if protein in proteins
])

#duplicating/tripling etc protein rows that appear in multiple pathways
diff_stats_CL = diff_stats_CL.explode("pathways").reset_index(drop=True)
#getting rid of nan rows
diff_stats_arid1a_breast = diff_stats_CL.dropna(subset=["pathways"]).reset_index(drop=True)
diff_stats_arid1a_breast

#only significant proteins for boxplot
#diff_stats_sf = diff_stats_arid1a.loc[diff_stats['adjusted_p'] < 0.0005]
#diff_stats_sf

#getting the colors for above of below median  - know what they are already
diff_stats_arid1a_breast = diff_stats_arid1a_breast.sort_values('pathways')

#finding out which ones are mean up or down regulation
mean_pathway_reg = diff_stats_arid1a_breast.groupby('pathways').diff_copy.mean().reset_index()
mean_pathway_reg

blue_pathways = ["KRAS signaling (down)", "coagulation", "myogenesis"]
diff_stats_arid1a_breast["Colour"] = diff_stats_arid1a_breast["pathways"].apply(
    lambda pathway: "Mean Downregulation" if pathway in blue_pathways else "Mean Upregulation"
)

plt.figure(figsize=(15,10))
sns.boxplot(data=diff_stats_arid1a_breast, x='pathways', y='diff_copy',
            palette={"Mean Downregulation": "skyblue", "Mean Upregulation": "red"},
            hue="Colour",
            dodge=False)

# Add custom legend
plt.title("Ovary Cancer", fontsize=20)
plt.xlabel(" ", fontsize=12)
plt.ylabel("Abudnance Change (LogFC)", fontsize=16)
plt.xticks(rotation=30, ha="right", fontsize=16)
plt.axhline(0, color='black', linestyle='--', linewidth=2)
plt.tight_layout()
plt.show()
