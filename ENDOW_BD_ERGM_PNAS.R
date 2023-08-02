library(igraph)
library(sna)
library(statnet)
library(intergraph)
library(stargazer)
library(arm)
library(latticeExtra)

# 1. setting up the data

## 1.1. loading node attribute data

node_info <- read.csv("https://raw.githubusercontent.com/JoonHwang-psu/ENDOW_BD_ERGM_PNAS/main/BD_SharingUnit_JH_R.csv", header=TRUE, fileEncoding="UTF-8-BOM")
node_info$name <- as.character(node_info$name)
node_info$logwealth <- log(node_info$WEALTH_TOTAL_VAL_USD, base=10)
node_info$head_age_sq <- node_info$head_age^2

### creating descriptive statistics of HH demography
descript_node <- node_info[,c("HH_size", "head_age", "F_age50plus", "age6under_bin", "WEALTH_TOTAL_VAL_USD", "land", "salaried", "Status")]

stargazer(descript_node, summary=TRUE, digits=2, 
          covariate.labels=c("Household size", "Household head age", "Women over age 50", "Children under age 6", "Household wealth (USD)", "Land ownership", "Salaried income", "Social status"),
          out="BD_ERGM_summary.htm")

## 1.2. adjacency matrix

### creating kinhsip network adjacency matrix
kin_edgelist<-read.csv("https://raw.githubusercontent.com/JoonHwang-psu/ENDOW_BD_ERGM_PNAS/main/BD_Edgelist_kinship_JH.csv", header=TRUE)
write.csv(as.matrix(get.adjacency(graph.data.frame(kin_edgelist, directed=FALSE))), file = "endow-BD/primary-sources/data/BD_adjmatrix_kinship_JH.csv", row.names = TRUE)
kin_adjmatrix<-read.csv("endow-BD/primary-sources/data/BD_adjmatrix_kinship_JH.csv", header=TRUE, row.names=1)
kin_adjmatrix <-as.matrix(kin_adjmatrix)
diag(kin_adjmatrix)<-0

### creating financial support network adjacency matrix
cash_edgelist<-read.csv("https://raw.githubusercontent.com/JoonHwang-psu/ENDOW_BD_ERGM_PNAS/main/BD_Edgelist_Q1Q2_no_ext_JH.csv", header=TRUE)
cash_graph_proto <- graph.data.frame(cash_edgelist, directed=TRUE)
cash_graph_proto <- cash_graph_proto + vertex("BDSU044") #adding vertex with 0 edge (not represented in adjacency matrix)
write.csv(as.matrix(get.adjacency(cash_graph_proto)), file = "endow-BD/primary-sources/data/BD_adjmatrix_Q1Q2_no_ext_JH.csv", row.names = TRUE)
cash_adjmatrix<-read.csv("endow-BD/primary-sources/data/BD_adjmatrix_Q1Q2_no_ext_JH.csv", header=TRUE, row.names=1)
cash_adjmatrix <-as.matrix(cash_adjmatrix)
diag(cash_adjmatrix)<-0 #deleting self-loops

### creating material support network adjacency matrix
material_edgelist<-read.csv("https://raw.githubusercontent.com/JoonHwang-psu/ENDOW_BD_ERGM_PNAS/main/BD_Edgelist_Q3Q4_no_ext_JH.csv", header=TRUE)
material_graph_proto <- graph.data.frame(material_edgelist, directed=TRUE)
material_graph_proto <- material_graph_proto + vertex("BDSU044")
write.csv(as.matrix(get.adjacency(material_graph_proto)), file = "endow-BD/primary-sources/data/BD_adjmatrix_Q3Q4_no_ext_JH.csv", row.names = TRUE)
material_adjmatrix<-read.csv("endow-BD/primary-sources/data/BD_adjmatrix_Q3Q4_no_ext_JH.csv", header=TRUE, row.names=1)
material_adjmatrix <-as.matrix(material_adjmatrix)
diag(material_adjmatrix)<-0

### creating labor support network adjacency matrix (inverse)
inv_help_edgelist<-read.csv("https://raw.githubusercontent.com/JoonHwang-psu/ENDOW_BD_ERGM_PNAS/main/BD_Edgelist_Q5Q6_no_ext_JH.csv", header=TRUE)
write.csv(as.matrix(get.adjacency(graph.data.frame(inv_help_edgelist, directed=TRUE))), file = "endow-BD/primary-sources/data/BD_adjmatrix_Q5Q6_no_ext_JH(inv).csv", row.names = TRUE)
inv_help_adjmatrix<-read.csv("endow-BD/primary-sources/data/BD_adjmatrix_Q5Q6_no_ext_JH(inv).csv", header=TRUE, row.names=1)
inv_help_adjmatrix <-as.matrix(inv_help_adjmatrix)
diag(inv_help_adjmatrix)<-0

### creating inverse matrix (to see the effects of between-network reciprocity)
inv_cash_adjmatrix <- t(cash_adjmatrix)
inv_material_adjmatrix <- t(material_adjmatrix)
help_adjmatrix <- t(inv_help_adjmatrix) #because the original help graph is inverse, re-inversing makes original graph where the direction of ties matches the flow of services

## 1.3. creating graphs from adj. matrices

### creating graphs from matrices
cash_graph<-graph_from_adjacency_matrix(cash_adjmatrix, mode="directed", weighted=TRUE, diag=FALSE)
material_graph<-graph_from_adjacency_matrix(material_adjmatrix, mode="directed", weighted=TRUE, diag=FALSE)
inv_help_graph<-graph_from_adjacency_matrix(inv_help_adjmatrix, mode="directed", weighted=TRUE, diag=FALSE)
kin_graph<-graph_from_adjacency_matrix(kin_adjmatrix, mode="undirected", weighted=TRUE, diag=FALSE)

### creating graphs with reversed ties (to see the effects of between-network reciprocity)
inv_cash_graph<-graph_from_adjacency_matrix(inv_cash_adjmatrix, mode="directed", weighted=TRUE, diag=FALSE)
inv_material_graph<-graph_from_adjacency_matrix(inv_material_adjmatrix, mode="directed", weighted=TRUE, diag=FALSE)
help_graph<-graph_from_adjacency_matrix(help_adjmatrix, mode="directed", weighted=TRUE, diag=FALSE)

## 1.4. assigning node attributes to graphs

### financial support network graph
V(cash_graph)$name
V(cash_graph)$head_age <- sapply(V(cash_graph)$name, function(x) node_info$head_age[node_info$name == x])
V(cash_graph)$head_age_sq <- sapply(V(cash_graph)$name, function(x) node_info$head_age_sq[node_info$name == x])
V(cash_graph)$land <- sapply(V(cash_graph)$name, function(x) node_info$land[node_info$name == x])
V(cash_graph)$salaried <- sapply(V(cash_graph)$name, function(x) node_info$salaried[node_info$name == x])
V(cash_graph)$highest_edu <- sapply(V(cash_graph)$name, function(x) node_info$highest_edu[node_info$name == x])
V(cash_graph)$avg_edu <- sapply(V(cash_graph)$name, function(x) node_info$avg_edu[node_info$name == x])
V(cash_graph)$HH_size <- sapply(V(cash_graph)$name, function(x) node_info$HH_size[node_info$name == x])
V(cash_graph)$wealth <- sapply(V(cash_graph)$name, function(x) node_info$WEALTH_TOTAL_VAL_USD[node_info$name == x])
V(cash_graph)$wealth1000 <- sapply(V(cash_graph)$name, function(x) node_info$wealth1000[node_info$name == x])
V(cash_graph)$logwealth <- sapply(V(cash_graph)$name, function(x) node_info$logwealth[node_info$name == x])
V(cash_graph)$labormigrant <- sapply(V(cash_graph)$name, function(x) node_info$LaborMigrant[node_info$name == x])
V(cash_graph)$producer <- sapply(V(cash_graph)$name, function(x) node_info$prod[node_info$name == x])
V(cash_graph)$male_producer <- sapply(V(cash_graph)$name, function(x) node_info$M_prod[node_info$name == x])
V(cash_graph)$female_producer <- sapply(V(cash_graph)$name, function(x) node_info$F_prod[node_info$name == x])
V(cash_graph)$female_age50plus <- sapply(V(cash_graph)$name, function(x) node_info$F_age50plus[node_info$name == x])
V(cash_graph)$CPratio <- sapply(V(cash_graph)$name, function(x) node_info$CPratio[node_info$name == x])
V(cash_graph)$status <- sapply(V(cash_graph)$name, function(x) node_info$Status[node_info$name == x])
V(cash_graph)$labormigrant <- sapply(V(cash_graph)$name, function(x) node_info$LaborMigrant[node_info$name == x])
V(cash_graph)$under19 <- sapply(V(cash_graph)$name, function(x) node_info$under19[node_info$name == x])
V(cash_graph)$age6under_bin <- sapply(V(cash_graph)$name, function(x) node_info$age6under_bin[node_info$name == x])

### material support network graph
V(material_graph)$name
V(material_graph)$head_age <- sapply(V(material_graph)$name, function(x) node_info$head_age[node_info$name == x])
V(material_graph)$head_age_sq <- sapply(V(material_graph)$name, function(x) node_info$head_age_sq[node_info$name == x])
V(material_graph)$land <- sapply(V(material_graph)$name, function(x) node_info$land[node_info$name == x])
V(material_graph)$salaried <- sapply(V(material_graph)$name, function(x) node_info$salaried[node_info$name == x])
V(material_graph)$highest_edu <- sapply(V(material_graph)$name, function(x) node_info$highest_edu[node_info$name == x])
V(material_graph)$avg_edu <- sapply(V(material_graph)$name, function(x) node_info$avg_edu[node_info$name == x])
V(material_graph)$HH_size <- sapply(V(material_graph)$name, function(x) node_info$HH_size[node_info$name == x])
V(material_graph)$wealth <- sapply(V(material_graph)$name, function(x) node_info$WEALTH_TOTAL_VAL_USD[node_info$name == x])
V(material_graph)$wealth1000 <- sapply(V(material_graph)$name, function(x) node_info$wealth1000[node_info$name == x])
V(material_graph)$logwealth <- sapply(V(material_graph)$name, function(x) node_info$logwealth[node_info$name == x])
V(material_graph)$under19 <- sapply(V(material_graph)$name, function(x) node_info$under19[node_info$name == x])
V(material_graph)$status <- sapply(V(material_graph)$name, function(x) node_info$Status[node_info$name == x])
V(material_graph)$age6under_bin <- sapply(V(material_graph)$name, function(x) node_info$age6under_bin[node_info$name == x])
V(material_graph)$female_age50plus <- sapply(V(material_graph)$name, function(x) node_info$F_age50plus[node_info$name == x])
V(material_graph)$labormigrant <- sapply(V(material_graph)$name, function(x) node_info$LaborMigrant[node_info$name == x])

### labor support network graph
V(help_graph)$name
V(help_graph)$head_age <- sapply(V(help_graph)$name, function(x) node_info$head_age[node_info$name == x])
V(help_graph)$head_age_sq <- sapply(V(help_graph)$name, function(x) node_info$head_age_sq[node_info$name == x])
V(help_graph)$land <- sapply(V(help_graph)$name, function(x) node_info$land[node_info$name == x])
V(help_graph)$salaried <- sapply(V(help_graph)$name, function(x) node_info$salaried[node_info$name == x])
V(help_graph)$highest_edu <- sapply(V(help_graph)$name, function(x) node_info$highest_edu[node_info$name == x])
V(help_graph)$avg_edu <- sapply(V(help_graph)$name, function(x) node_info$avg_edu[node_info$name == x])
V(help_graph)$HH_size <- sapply(V(help_graph)$name, function(x) node_info$HH_size[node_info$name == x])
V(help_graph)$wealth <- sapply(V(help_graph)$name, function(x) node_info$WEALTH_TOTAL_VAL_USD[node_info$name == x])
V(help_graph)$wealth1000 <- sapply(V(help_graph)$name, function(x) node_info$wealth1000[node_info$name == x])
V(help_graph)$logwealth <- sapply(V(help_graph)$name, function(x) node_info$logwealth[node_info$name == x])
V(help_graph)$labormigrant <- sapply(V(help_graph)$name, function(x) node_info$LaborMigrant[node_info$name == x])
V(help_graph)$producer <- sapply(V(help_graph)$name, function(x) node_info$prod[node_info$name == x])
V(help_graph)$male_producer <- sapply(V(help_graph)$name, function(x) node_info$M_prod[node_info$name == x])
V(help_graph)$female_producer <- sapply(V(help_graph)$name, function(x) node_info$F_prod[node_info$name == x])
V(help_graph)$female_age50plus <- sapply(V(help_graph)$name, function(x) node_info$F_age50plus[node_info$name == x])
V(help_graph)$age6under <- sapply(V(help_graph)$name, function(x) node_info$age6under[node_info$name == x])
V(help_graph)$age6under_bin <- sapply(V(help_graph)$name, function(x) node_info$age6under_bin[node_info$name == x])
V(help_graph)$CPratio<- sapply(V(help_graph)$name, function(x) node_info$CPratio[node_info$name == x])
V(help_graph)$status<- sapply(V(help_graph)$name, function(x) node_info$Status[node_info$name == x])
V(help_graph)$labormigrant<- sapply(V(help_graph)$name, function(x) node_info$LaborMigrant[node_info$name == x])
V(help_graph)$under19 <- sapply(V(help_graph)$name, function(x) node_info$under19[node_info$name == x])
V(help_graph)$status <- sapply(V(help_graph)$name, function(x) node_info$Status[node_info$name == x])

## 1.5. creating networks from graphs

cash_network <- asNetwork(cash_graph)
material_network <- asNetwork(material_graph)
help_network <- asNetwork(help_graph)
kin_network <- asNetwork(kin_graph)

### creating reversed networks
inv_cash_network <-asNetwork(inv_cash_graph)
inv_material_network <-asNetwork(inv_material_graph)
inv_help_network <-asNetwork(inv_help_graph)

### checking networks

#### financial support network
summary(cash_network)
get.node.attr(cash_network, 'salaried')
get.node.attr(cash_network, 'avg_edu')
get.node.attr(cash_network, 'wealth')
get.node.attr(cash_network, 'HH_size')
get.node.attr(cash_network, 'female_producer')

#### material network
summary(material_network)
get.node.attr(material_network, 'salaried')
get.node.attr(material_network, 'avg_edu')
get.node.attr(material_network, 'wealth')
get.node.attr(material_network, 'HH_size')
get.node.attr(material_network, 'female_producer')

#### labor support network
summary(help_network)
get.node.attr(help_network, 'salaried')
get.node.attr(help_network, 'avg_edu')
get.node.attr(help_network, 'wealth')
get.node.attr(help_network, 'producer')
get.node.attr(help_network, 'female_producer')

# 2. ERGM

## 2.1. financial support network ERGM

cash.ERGM <- ergm(cash_network~edges+nodecov("HH_size")+nodecov("head_age")+nodecov("head_age_sq")
                 +nodeofactor("female_age50plus")+nodeifactor("age6under_bin")
                 +nodeifactor("land")+nodeofactor("land")+nodematch("land")
                 +nodeifactor("salaried")+nodeofactor("salaried")+nodematch("salaried")
                 +nodeifactor("status")+nodeofactor("status")+nodematch("status")
                 +diff("wealth1000", pow=1, dir="t-h", sign.action="identity")
                 +edgecov(inv_help_network)+edgecov(inv_material_network)+edgecov(kin_network)
                 +mutual+dgwesp(0.5, T)+dgwdsp(0.5, T)
                 ,control=control.ergm(MCMC.burnin=20000, MCMC.samplesize=10000)
                 )
summary(cash.ERGM)

### MCMC diagnostics and GOF

mcmc.diagnostics(cash.ERGM)
GOF.cash.ERGM<-gof(cash.ERGM)
plot(GOF.cash.ERGM)

coefplot(ERGM2.18.1)

## 2.2. material support network ERGM

material.ERGM <- ergm(material_network~edges+nodecov("HH_size")+nodecov("head_age")+nodecov("head_age_sq")
                      +nodeofactor("female_age50plus")+nodeifactor("age6under_bin")
                      +nodeofactor("salaried")+nodeifactor("salaried")+nodematch("salaried")
                      +nodeofactor("land")+nodeifactor("land")+nodematch("land")
                      +nodeofactor("status")+nodeifactor("status")+nodematch("status")
                      +diff("wealth1000", pow=1, dir="t-h", sign.action="identity")
                      +edgecov(inv_cash_network)+edgecov(inv_help_network)+edgecov(kin_network)
                      +mutual+dgwesp(0.5, T)+dgwdsp(0.5,T)
                      ,control=control.ergm(MCMC.burnin=20000, MCMC.samplesize=10000)
)
summary(material.ERGM)

### MCMC diagnostics and GOF

mcmc.diagnostics(material.ERGM)
GOF.material.ERGM<-gof(material.ERGM)
plot(GOF.material.ERGM)

## 2.3. labor support network ERGM

help.ERGM <- ergm(help_network~edges+nodecov("HH_size")+nodecov("head_age")+nodecov("head_age_sq")
                   +nodeofactor("female_age50plus")+nodeifactor("age6under_bin")
                   +nodeofactor("salaried")+nodeifactor("salaried")+nodematch("salaried")
                   +nodeofactor("land")+nodeifactor("land")+nodematch("land")
                   +nodeofactor("status")+nodeifactor("status")+nodematch("status")
                   +diff("wealth1000", pow=1, dir="t-h", sign.action="identity")
                   +edgecov(inv_cash_network)+edgecov(inv_material_network)+edgecov(kin_network)
                   +mutual+dgwesp(0.5,T)+dgwdsp(0.5,T)
                   ,control=control.ergm(MCMC.burnin=20000, MCMC.samplesize=10000)
                  )
summary(help.ERGM)

### MCMC diagnostics and GOF

mcmc.diagnostics(help.ERGM)
GOF.help.ERGM<-gof(help.ERGM)
plot(GOF.help.ERGM)

dev.off()

### creating ERGM result table

par(mar=c(1,1,1,1))

stargazer(cash.ERGM, material.ERGM, help.ERGM, star.cutoffs=c(.05, .01, .001), out="BDSU_ERGM_table.htm")


