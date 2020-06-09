#!/usr/bin/env Rscript

require(jlp)

con <- file("stdin")
data <- jlp::readpdb(con)
close(con)

data[,"Chain.Type"] <- "G"

writepdb(data)
