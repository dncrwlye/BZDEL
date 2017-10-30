## phylofactor bin classification script
## Daniel Becker
## daniel.becker3@montana.edu

## clear workspace
rm(list=ls()) 
graphics.off()

## packages
library(readxl) ## alternative if xlsx fails (25 July 2017)

## load in data
setwd("~/Desktop/BZDEL/Data/Olival/data")
pfactor=read_excel("Olival_w_phylobin.xlsx")

## merge in zoonotic virus trait daat
