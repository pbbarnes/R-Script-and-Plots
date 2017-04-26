################
## Paul Barnes##
##   CS 422   ##
## Spring 2017##
################
library(Rgraphviz)
library(gRain)
library(gRbase)
library(bnlearn)
library(caTools)
library(readr)

#loads external dataset and changes the data types for each column 
dermatology <- read_csv("~/Dropbox/SCHOOL/Spring 2017/CS 422/Bayesian Belief Networks/dermatology.csv", 
col_types = cols(A = col_number(), AGE = col_number(), BLI = col_number(), C = col_number(), CRR = col_number(), D = col_number(), DB = col_number(), DGL = col_number(), E = col_number(), EI = col_number(), 
ERR = col_number(), ES = col_number(), F = col_number(), FH = col_number(), FHP = col_number(), FH_1 = col_number(), FP = col_number(), FPD = col_number(), G = col_number(), HK = col_number(), 
I = col_number(), IMI = col_number(), KEI = col_number(), KP = col_number(), MI = col_number(), MM = col_number(), OMI = col_number(), PK = col_number(), PNL = col_number(), PP = col_number(), PPS = col_number(), S = col_number(), SI = col_number(), SP = col_number(), SS = col_number(), STAR = col_number(), TSE = col_number(), VDBL = col_number(), AI = col_number()))

# Uses Grow-Shrink Algorithm
res = gs(dermatology)
graphviz.plot(res)
res
# Using graphviz plot check what fit is best in DAG undirected nodes
choose.direction(res, dermatology, arc = c("MI", "BLI"), debug = TRUE)
choose.direction(res, dermatology, arc = c("SI", "CRR"), debug = TRUE)
choose.direction(res, dermatology, arc = c("MI", "FH_1"), debug = TRUE)
choose.direction(res, dermatology, arc = c("AI", "PNL"), debug = TRUE)
choose.direction(res, dermatology, arc = c("I", "FPD"), debug = TRUE)
choose.direction(res, dermatology, arc = c("SI", "MM"), debug = TRUE)
# Removes or sets the direction of arcs
res <- drop.edge(res, "I", "FPD")
res <- drop.arc(res, "I", "FPD")
res <- set.arc(res, "CRR", "SI")
res <- set.arc(res, "SI", "MM")
res <- set.arc(res, "BLI", "MI")
res <- drop.arc(res, "AI", "PNL")
res <- drop.edge(res, "AI", "PNL")
res <- set.arc(res, "MI", "FH_1")
res
# Plots the new graph after manipulating graph
graphviz.plot(res)
# Checks arc strength to determine strength or relation
arc.strength(res, dermatology)
check <- arc.strength(res, dermatology)
# Organizes Strength by ascending order
check[order(check$strength),]



# Uses Incremental Association Algorithm 
res2 = iamb(dermatology)
res2
# Plots second algorithm bayes net
graphviz.plot(res2)
# Checks which direction would be best for second implementation
choose.direction(res2, dermatology, arc = c("I", "KP"), debug = TRUE)
choose.direction(res2, dermatology, arc = c("FH_1", "BLI"), debug = TRUE)
choose.direction(res2, dermatology, arc = c("FPD", "PK"), debug = TRUE)
choose.direction(res2, dermatology, arc = c("DB", "PK"), debug = TRUE)
choose.direction(res2, dermatology, arc = c("KP", "FPD"), debug = TRUE)
choose.direction(res2, dermatology, arc = c("SI", "MM"), debug = TRUE)
choose.direction(res2, dermatology, arc = c("SI", "CRR"), debug = TRUE)
choose.direction(res2, dermatology, arc = c("SI", "KEI"), debug = TRUE)
choose.direction(res2, dermatology, arc = c("DB", "KP"), debug = TRUE)
choose.direction(res2, dermatology, arc = c("AI", "PNL"), debug = TRUE)
choose.direction(res2, dermatology, arc = c("OMI", "SS"), debug = TRUE)
choose.direction(res2, dermatology, arc = c("I", "FPD"), debug = TRUE)
# Modifies the iamb bayesian network to be a DAG
res2 <- drop.arc(res2, "I", "KP")
res2 <- drop.edge(res2, "I", "KP")
res2 <- set.arc(res2, "BLI", "FH_1")
res2 <- set.arc(res2, "SI", "MM")
res2 <- set.arc(res2, "CRR", "SI")
res2 <- set.arc(res2, "SI", "KEI")
res2 <- set.arc(res2, "AI", "PNL")
res2 <- drop.arc(res2, "FPD", "PK")
res2 <- drop.edge(res2, "FPD", "PK")
res2 <- drop.arc(res2, "DB", "PK")
res2 <- drop.edge(res2, "DB", "PK")
res2 <- drop.arc(res2, "KP", "FPD")
res2 <- drop.edge(res2, "KP", "FPD")
res2 <- drop.arc(res2, "DB", "KP")
res2 <- drop.edge(res2, "DB", "KP")
res2 <- drop.arc(res2, "OMI", "SS")
res2 <- drop.edge(res2, "OMI", "SS")
res2 <- drop.edge(res2, "I", "FPD")
res2
graphviz.plot(res2)

# Greedy search using Hill Climbing algorithm
res3 = hc(dermatology)
res3
graphviz.plot(res3)


# Counts the number of tests that the algorithm uses
ntests(res)
ntests(res2)
ntests(res3)

# # Allows us to watch the algorithms run
# res = gs(dermatology, debug = TRUE)
# res2 = gs(dermatology, debug = TRUE)
# res3 = hc(dermatology, debug = TRUE)

# Fit the parameters of a Bayesian network conditional on its structure.
fitted = bn.fit(hc(dermatology), dermatology) 
# Plot the new fitted bayes net
graphviz.plot(fitted)

# Queries of interest on conditional probability
cpquery(fitted, event = ((PK > 2) & (PK <= 3)), (evidence = (AI <= 5) & (AI >4)))
cpquery(fitted, event = ((PK > 2) & (PK <= 3)), (evidence = (AI <= 4) & (AI >3)))
cpquery(fitted, event = ((PK > 2) & (PK <= 3)), (evidence = (AI <= 3) & (AI >2)))
cpquery(fitted, event = ((PK > 2) & (PK <= 3)), (evidence = (AI <= 2) & (AI >1)))
cpquery(fitted, event = ((PK > 2) & (PK <= 3)), (evidence = (AI <= 1) & (AI >0)))