#!/usr/bin/Rscript

###
library(BiocManager)
library(memoise)
library(sessioninfo)
library(devtools)
library(methylKit)
library(edmr)

file.list<-list("A1.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "A2.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "A3.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "A4.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "A5.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "A6.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "B1.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "B2.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "B4.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "B5.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "B6.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "C2.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "C3.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "C4.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "C6.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "D1.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "D2.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "D3.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "D4.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "D5.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "D6.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "E1.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "E2.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "E3.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "E4.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "E5.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "F_B5.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "F1.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "F2.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "F3.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "F4.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "F5.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "F6.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "G1.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "G2.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "G4.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "G5.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "G6.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "H1.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "H2.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "H3.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz",
                    "H4.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz")
                    
IDs<-list("A1","A2","A3","A4","A5","A6","B1","B2","B4","B5","B6","C2","C3","C4","C6","D1","D2","D3","D4","D5","D6","E1","E2","E3","E4","E5","F_B5" "F1","F2","F3","F4","F5","F6","G1","G2","G4","G5","G6","H1","H2","H3","H4")                    
treatment<-c(0,1,0,0,1,0,1,1,0,1,1,1,0,0,1,1,0,0,1,1,0,0,0,1,1,1,0,1,0,1,0,0,1,1,1,1,0,0,1,1,0,1)
my_obj<-methRead(file.list, sample.id = IDs, assembly = "Chick", pipeline = "bismarkCoverage", treatment= treatment, context= "CpG", mincov=5)
#t
