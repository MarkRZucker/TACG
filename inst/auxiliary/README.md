---
author: Mark Zucker
date: 14 April 2019
title: Auxiliary Material for TACG
---

# Overview
This directory contains a combination of material that performs
some of the analyses described in our manuscript.


## Analyses

1. `chlens.txt` is a simple table holding the lengths In base pairs)
   of human chromosomes. It is based on build [UNKNOWN} and was
   initialy obtained from [UNKNOWN].
2. `makeChlens.R` is a script that reads the previous file, converts
   it to a useful form, and saves it in the `sysdata.rda` file so it
   can be used by the package. In particular, it is used by functions
   that simulate `Tumor` objects to create  realistic segmentation
   data.
3. `generateSimulationSet.R` is a script that generates sets of 
   simulated chromosomes and associated data to which segmentation 
   algorithms can be applied, such as those used in our manuscript.
4. `loci.rda` consists of a set of SNP loci for the simulated 
   chromosomes used in our analysis.
