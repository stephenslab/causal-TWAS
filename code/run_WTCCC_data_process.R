source("input_reformat.R")

N_per_chunk = 10000 # ussed to save memory when writing chunks
gemma.prepare("/home/simingz/causalTWAS/WTCCC/bd.RData", "/home/simingz/causalTWAS/WTCCC/bd")
