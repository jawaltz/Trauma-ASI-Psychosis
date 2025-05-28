library("RColorBrewer")
library("readxl")
library(igraph)

#####################################################################################
#####################################################################################
################################## preprocessing ####################################
#####################################################################################
#####################################################################################

var <- readxl::read_xlsx("Master_CTQ_ASI_BPRS.xlsx", sheet = 1)
dim(var)
var = var[-which(var[,3]=="FDR"),]
names <- names(var)[-c(1:4)]

var1 = data.frame(var[,5:81])

cormat = matrix(1,dim(var1)[2],dim(var1)[2])

for(i in 1:(dim(var1)[2]-1)){
  for(j in (i+1):dim(var1)[2]){
    cormat[i,j] = cor(na.omit(data.frame(as.numeric(var1[,i]),as.numeric(var1[,j]) )))[1,2]
    cormat[j,i]=cormat[i,j]
  }
}

threshold= 0.34
W = cormat

lambda =seq(1,2,by=0.1)

result_00 = lapply(lambda, function(x) dense(abs(W),threshold=threshold,lambda=x))


W00=matrix(0,dim(cormat)[1],dim(cormat)[2])
W00[W>threshold] = 1
diag(W00) = 0
kld_values = unlist(lapply(result_00, function(x) kld(W00,x$Clist,x$CID,length(x$CID))))

max_index = which.max(kld_values)
result = result_00[[max_index]]
Clist = result$Clist
CID0 = result$CID


########################################################
######################## SZ ########################
########################################################

var_SZ = data.frame(var[var$group=="SZ",5:81])

cormat_SZ = matrix(1,dim(var_SZ)[2],dim(var_SZ)[2])

for(i in 1:(dim(var_SZ)[2]-1)){
  for(j in (i+1):dim(var_SZ)[2]){
    cormat_SZ[i,j] = cor(na.omit(data.frame(as.numeric(var_SZ[,i]),as.numeric(var_SZ[,j]) )))[1,2]
    cormat_SZ[j,i]=cormat_SZ[i,j]
  }
}

W_SZ = cormat_SZ
results_SZ=dense(abs(W_SZ),threshold=0.34,lambda=1.25,Advanced=TRUE)
CID0_SZ= results_SZ$CID
Clist_SZ=results_SZ$Clist
CIDsum_SZ = cumsum(results_SZ$CID) 

########################################################
######################## SZ- high hal ########################
########################################################

var_SZ = data.frame(var[var$group=="SZ",5:81])

var_sz_h_ind = var_SZ$X12_Hallucination>=3
var_sz_l_ind = var_SZ$X12_Hallucination<=2

var0 <- readxl::read_xlsx("/Users/Hongju.Park/Dropbox/Sonia_network/CTQ_BPRS_ASI domains scores for SEM.xlsx", sheet = 1,na="ms")
#var0_SZ= var0[var0$group=="SZ",]
#var0_SZ_hp_ind = var0_SZ$`Positive-BPRS` > median(var0_SZ$`Positive-BPRS`[!is.na(var0_SZ$`Positive-BPRS`)])
#var0_SZ_lp_ind = var0_SZ$`Positive-BPRS` <= median(var0_SZ$`Positive-BPRS`[!is.na(var0_SZ$`Positive-BPRS`)])
#table(var_sz_h_ind,var0_SZ_hp_ind)


sum((as.numeric(var0_SZ$`Positive-BPRS` == median(var0_SZ$`Positive-BPRS`[!is.na(var0_SZ$`Positive-BPRS`)]))))

var_SZ_h = var_SZ[var_SZ$X12_Hallucination>=3,]


cormat_SZ_h = matrix(1,dim(var_SZ_h)[2],dim(var_SZ_h)[2])

for(i in 1:(dim(var_SZ_h)[2]-1)){
  for(j in (i+1):dim(var_SZ_h)[2]){
    cormat_SZ_h[i,j] = cor(na.omit(data.frame(as.numeric(var_SZ_h[,i]),as.numeric(var_SZ_h[,j]) )))[1,2]
    cormat_SZ_h[j,i]=cormat_SZ_h[i,j]
  }
}

W_SZ_h = cormat_SZ_h

results_SZ_h=dense(abs(W_SZ_h),threshold=0.34,lambda=1.2,Advanced=TRUE)
CID0_SZ_h= results_SZ_h$CID
Clist_SZ_h=results_SZ_h$Clist
CIDsum_SZ_h = cumsum(results_SZ_h$CID) 



########################################################
######################## SZ- low hal ########################
########################################################

var_SZ = data.frame(var[var$group=="SZ",5:81])
var_SZ_l = var_SZ[var_SZ$X12_Hallucination<=2,]


cormat_SZ_l = matrix(1,dim(var_SZ_l)[2],dim(var_SZ_l)[2])

for(i in 1:(dim(var_SZ_l)[2]-1)){
  for(j in (i+1):dim(var_SZ_l)[2]){
    cormat_SZ_l[i,j] = cor(na.omit(data.frame(as.numeric(var_SZ_l[,i]),as.numeric(var_SZ_l[,j]) )))[1,2]
    cormat_SZ_l[j,i]=cormat_SZ_l[i,j]
  }
}

W_SZ_l = cormat_SZ_l

results_SZ_l=dense(abs(W_SZ_l),threshold=0.34,lambda=1.2,Advanced=TRUE)
CID0_SZ_l= results_SZ_l$CID
Clist_SZ_l=results_SZ_l$Clist
CIDsum_SZ_l = cumsum(results_SZ_l$CID) 


########################################################
######################## HC ########################
########################################################

var_HC = data.frame(var[var$group=="HC",5:81])
var_HC[,42]
var_HC[,51]
cormat_HC = diag(dim(var_HC)[2])

for(i in 1:(dim(var_HC)[2]-1)){
  for(j in (i+1):dim(var_HC)[2]){
    cormat_HC[i,j] = cor(na.omit(data.frame(as.numeric(var_HC[,i]),as.numeric(var_HC[,j]) )))[1,2]
    cormat_HC[j,i]=cormat_HC[i,j]
  }
}
W_HC = cormat_HC
W_HC[W_HC==NA] = 0
W_HC[, c(35,46)] = 0
W_HC[c(35,46),] = 0
W_HC[49:77,c(40,41)] = 0
W_HC[c(40,41),49:77] = 0



results_HC=dense(abs(W_HC),threshold=0.34,lambda=1.2)
CID0_HC= results_HC$CID
CID0_HC = c(CID0_HC[1:7],sum(CID0_HC[-(1:7)]))
Clist_HC=results_HC$Clist
CIDsum_HC = cumsum(results_HC$CID) 



########################################################
######################## NCVH ########################
########################################################

var_NCVH = data.frame(var[var$group=="NCVH",5:81])
cormat_NCVH = diag(dim(var_NCVH)[2])


for(i in 1:(dim(var_NCVH)[2]-1)){
  for(j in (i+1):dim(var_NCVH)[2]){
    cormat_NCVH[i,j] = cor(na.omit(data.frame(as.numeric(var_NCVH[,i]),as.numeric(var_NCVH[,j]) )))[1,2]
    cormat_NCVH[j,i] = cormat_NCVH[i,j]
  }
}


cormat_NCVH[,c(46,47)] = 0
cormat_NCVH[c(46,47),] = 0

cormat_NCVH[c(30,34,37,38),c(30,34,37,38)]
cormat_NCVH[c(45,48),c(45,48)]

W_NCVH = cormat_NCVH

results_NCVH=dense(abs(W_NCVH),threshold=0.34,lambda=1.25)
CID0_NCVH = results_NCVH$CID
CID0_NCVH = c(CID0_NCVH[1:7],sum(CID0_NCVH[-(1:7)]))
Clist_NCVH=results_NCVH$Clist
CIDsum_NCVH = cumsum(results_NCVH$CID) 

###################################################################################################
###################################################################################################
###################################################################################################
########################################## igraph ################################################
###################################################################################################
###################################################################################################
###################################################################################################




connection = matrix(0,77,77)
W0 = matrix(0,77,77)
W0[abs(W)>0.34] = 1
sum(W0)/77^2
ind=1


for(i in 1:length(CID0)){
  connection[Clist[ind:(ind+CID0[i]-1)],Clist[ind:(ind+CID0[i]-1)]]=1
  ind = ind + CID0[i]
}


con2 = c()
ind=1
for(i in 1:76){
  for(j in (i+1):77){
    if(connection[i,j] == 1){
      con2 = c(con2,i,j)
    }
  } 
}

con3 = c()
ind=1
strength = c()
for(i in 1:76){
  for(j in (i+1):77){
    if(W0[i,j] == 1){
      con3 = c(con3,i,j)
      strength = c(strength,abs(W[i,j])) 
    }
  } 
}

colors0 = c()
for(i in 1:76){
  for(j in (i+1):77){
    if(W0[i,j] == 1){
      colors0 = c(colors0,sign(W[i,j]) ) 
    }
  } 
}



edges = matrix(con3,ncol=2,byrow=T)
edges = cbind(edges,strength)


colnames(connection) = names



names2 = c()

CID = rep(1:8,times=c(CID0[1],CID0[2],CID0[3],CID0[4],CID0[5],CID0[6],CID0[7],sum(CID0[-(1:7)]) ))

for(i in 1:77){
  if (Clist[i]<10){
    names2[i] = paste("  ",Clist[i],": ",names[Clist[i]],sep="")
  } else {
    names2[i] = paste(Clist[i],": ",names[Clist[i]],sep="")  
  }
}


g <- make_empty_graph()
g <- make_graph(edges = con3, n = 77, directed = FALSE)
col00 = c(brewer.pal(n = 11, name = "RdYlGn")[3],brewer.pal(n = 11, name = "RdYlGn")[4],brewer.pal(n = 11, name = "RdYlGn")[5],brewer.pal(n = 8, name = "Set2")[5],
          brewer.pal(n = 11, name = "Spectral")[8],brewer.pal(n = 11, name = "BrBG")[8],brewer.pal(n = 8, name = "Set2")[3],"#B3B3B3")
V(g)$color[Clist] <- col00[CID]
E(g)$width <- strength/mean(strength)
E(g)$color <- colors0
E(g)$color[E(g)$color==1] = brewer.pal(n = 9, name = "Reds")[4]
E(g)$color[E(g)$color==-1] = brewer.pal(n = 9, name = "Blues")[4]

#minC <- rep(-Inf, vcount(g))
#maxC <- rep(Inf, vcount(g))
#minC[1] <- maxC[1] <- 0
#location <- layout_with_fr(g, minx=minC, maxx=maxC,miny=minC, maxy=maxC)
#location=ig_location(location,c(13),0.1,0.1)
#setwd("/Users/hongju/Downloads/Dropbox/igraph")
#write.csv(location,"location_total.csv")

location0=read.csv("/Users/Hongju.Park/Dropbox/Sonia_network/igraph/location_total.csv")[,2:3]
location0=as.matrix(location0)

location =location0

#location=ig_location(location,c(13),0,0.1)
#location0 = location
plot(g, layout=location0, vertex.size = 7,vertex.label.dist = 0, vertex.label.cex=0.8,edge.width=E(g)$width^3,main="Total")


#edge.width = ifelse(E(g)$is_formal

legend("topleft",legend=names2[1:38],col=col00[CID][1:38],pch=16,cex=0.7)
legend("topright",legend=names2[39:77],col=col00[CID][39:77],pch=16,cex=0.7)




########################################################################################
########################################################################################
############################### SZ #######################################################
########################################################################################
########################################################################################


connection_SZ = matrix(0,77,77)


W0_SZ= matrix(0,77,77)
W0_SZ[abs(W_SZ)>0.34] = 1
sum(W0_SZ)/77^2
ind=1


for(i in 1:length(CID0_SZ)){
  connection[Clist_SZ[ind:(ind+CID0_SZ[i]-1)],Clist_SZ[ind:(ind+CID0_SZ[i]-1)]]=1
  ind = ind + CID0[i]
}


con2_SZ = c()
ind=1
for(i in 1:76){
  for(j in (i+1):77){
    if(connection_SZ[i,j] == 1){
      con2_SZ = c(con2_SZ,i,j)
    }
  } 
}

con3_SZ = c()
ind=1
strength_SZ = c()
for(i in 1:76){
  for(j in (i+1):77){
    if(W0_SZ[i,j] == 1){
      con3_SZ = c(con3_SZ,i,j)
      strength_SZ = c(strength_SZ,abs(W_SZ[i,j])) 
    }
  } 
}


colors0 = c()
for(i in 1:76){
  for(j in (i+1):77){
    if(W0_SZ[i,j] == 1){
      colors0 = c(colors0,sign(W_SZ[i,j]) ) 
    }
  } 
}

sum(W_SZ< -0.34)

#edges = matrix(con3,ncol=2,byrow=T)
#edges = cbind(edges,strength)

colnames(connection_SZ) = names

dim(connection)


names2_SZ = c()

CID_SZ = rep(c(1:8),times=c(CID0_SZ[1],CID0_SZ[2],CID0_SZ[3],CID0_SZ[4],CID0_SZ[5],CID0_SZ[6],CID0_SZ[7],sum(CID0_SZ[-(1:7)]) ))

for(i in 1:77){
  if (Clist_SZ[i]<10){
    names2_SZ[i] = paste("  ",Clist_SZ[i],": ",names[Clist_SZ[i]],sep="")
  } else {
    names2_SZ[i] = paste(Clist_SZ[i],": ",names[Clist_SZ[i]],sep="")  
  }
}

g_SZ <- make_empty_graph()
g_SZ <- make_graph(edges = con3_SZ, n = 77, directed = FALSE)
V(g_SZ)$color[Clist_SZ] <- col00[CID_SZ]
E(g_SZ)$width <- strength_SZ/mean(strength_SZ)
E(g_SZ)$color <- colors0
E(g_SZ)$color[E(g_SZ)$color==1] = brewer.pal(n = 9, name = "Reds")[4]
E(g_SZ)$color[E(g_SZ)$color==-1] = brewer.pal(n = 9, name = "Blues")[4]



minC_SZ <- rep(-Inf, vcount(g_SZ))
maxC_SZ <- rep(Inf, vcount(g_SZ))
minC_SZ[1] <- maxC_SZ[1] <- 0
location_SZ <- layout_with_fr(g_SZ, minx=minC_SZ, maxx=maxC_SZ,miny=minC_SZ, maxy=maxC_SZ)

plot(g_SZ, layout=location, vertex.size = 7,vertex.label.dist = 0, vertex.label.cex=0.8,edge.width=E(g_SZ)$width^3,main="SZ")
legend("topleft",legend=names2_SZ[1:40],col=col00[CID_SZ][1:40],pch=16,cex=0.7)
legend("topright",legend=names2_SZ[41:77],col=col00[CID_SZ][41:77],pch=16,cex=0.7)


plot(g_SZ,  vertex.size = 7,vertex.label.dist = 0, vertex.label.cex=0.8,edge.width=E(g_SZ)$width^3,main="SZ")
legend("topleft",legend=names2_SZ[1:40],col=col00[CID_SZ][1:40],pch=16,cex=0.7)
legend("topright",legend=names2_SZ[41:77],col=col00[CID_SZ][41:77],pch=16,cex=0.7)



########################################################################################
########################################################################################
############################### SZ high #######################################################
########################################################################################
########################################################################################


connection_SZ_h = matrix(0,77,77)


W0_SZ_h= matrix(0,77,77)
W0_SZ_h[abs(W_SZ_h)>0.34] = 1
sum(W0_SZ_h)/77^2
ind=1


for(i in 1:length(CID0_SZ_h)){
  connection[Clist_SZ_h[ind:(ind+CID0_SZ_h[i]-1)],Clist_SZ[ind:(ind+CID0_SZ_h[i]-1)]]=1
  ind = ind + CID0[i]
}


con2_SZ_h = c()
ind=1
for(i in 1:76){
  for(j in (i+1):77){
    if(connection_SZ_h[i,j] == 1){
      con2_SZ_h = c(con2_SZ_h,i,j)
    }
  } 
}

con3_SZ_h = c()
ind=1
strength_SZ_h = c()
for(i in 1:76){
  for(j in (i+1):77){
    if(W0_SZ_h[i,j] == 1){
      con3_SZ_h = c(con3_SZ_h,i,j)
      strength_SZ_h = c(strength_SZ_h,abs(W_SZ_h[i,j])) 
    }
  } 
}


colors0 = c()
for(i in 1:76){
  for(j in (i+1):77){
    if(W0_SZ_h[i,j] == 1){
      colors0 = c(colors0,sign(W_SZ_h[i,j]) ) 
    }
  } 
}

#edges = matrix(con3,ncol=2,byrow=T)
#edges = cbind(edges,strength)

colnames(connection_SZ_h) = names
names2_SZ_h = c()

CID_SZ_h = rep(c(1:8),times=c(CID0_SZ_h[1],CID0_SZ_h[2],CID0_SZ_h[3],CID0_SZ_h[4],CID0_SZ_h[5],CID0_SZ_h[6],CID0_SZ_h[7],sum(CID0_SZ_h[-(1:7)]) ))

for(i in 1:77){
  if (Clist_SZ_h[i]<10){
    names2_SZ_h[i] = paste("  ",Clist_SZ_h[i],": ",names[Clist_SZ_h[i]],sep="")
  } else {
    names2_SZ_h[i] = paste(Clist_SZ_h[i],": ",names[Clist_SZ_h[i]],sep="")  
  }
}

g_SZ_h <- make_empty_graph()
g_SZ_h <- make_graph(edges = con3_SZ_h, n = 77, directed = FALSE)
V(g_SZ_h)$color[Clist_SZ_h] <- col00[CID_SZ_h]
E(g_SZ_h)$width <- strength_SZ_h/mean(strength_SZ_h)
E(g_SZ_h)$color <- colors0
E(g_SZ_h)$color[E(g_SZ_h)$color==1] = brewer.pal(n = 9, name = "Reds")[4]
E(g_SZ_h)$color[E(g_SZ_h)$color==-1] = brewer.pal(n = 9, name = "Blues")[4]



minC_SZ_h <- rep(-Inf, vcount(g_SZ_h))
maxC_SZ_h <- rep(Inf, vcount(g_SZ_h))
minC_SZ_h[1] <- maxC_SZ_h[1] <- 0
location_SZ_h <- layout_with_fr(g_SZ_h, minx=minC_SZ_h, maxx=maxC_SZ_h,miny=minC_SZ_h, maxy=maxC_SZ_h)

plot(g_SZ_h, layout=location0, vertex.size = 7,vertex.label.dist = 0, vertex.label.cex=0.8,edge.width=E(g_SZ_h)$width^3,main="SZ high-Hallucination")
legend("topleft",legend=names2_SZ_h[1:40],col=col00[CID_SZ_h][1:40],pch=16,cex=0.7)
legend("topright",legend=names2_SZ_h[41:77],col=col00[CID_SZ_h][41:77],pch=16,cex=0.7)


plot(g_SZ_h,  vertex.size = 7,vertex.label.dist = 0, vertex.label.cex=0.8,edge.width=E(g_SZ_h)$width^3,main="SZ high-Hallucination")
legend("topleft",legend=names2_SZ_h[1:40],col=col00[CID_SZ_h][1:40],pch=16,cex=0.7)
legend("topright",legend=names2_SZ_h[41:77],col=col00[CID_SZ_h][41:77],pch=16,cex=0.7)


########################################################################################
########################################################################################
############################### SZ low #######################################################
########################################################################################
########################################################################################


connection_SZ_l = matrix(0,77,77)


W0_SZ_l= matrix(0,77,77)
W0_SZ_l[abs(W_SZ_l)>0.34] = 1
sum(W0_SZ_l)/77^2
ind=1


for(i in 1:length(CID0_SZ_l)){
  connection[Clist_SZ_l[ind:(ind+CID0_SZ_l[i]-1)],Clist_SZ[ind:(ind+CID0_SZ_l[i]-1)]]=1
  ind = ind + CID0[i]
}


con2_SZ_l = c()
ind=1
for(i in 1:76){
  for(j in (i+1):77){
    if(connection_SZ_l[i,j] == 1){
      con2_SZ_l = c(con2_SZ_l,i,j)
    }
  } 
}

con3_SZ_l = c()
ind=1
strength_SZ_l = c()
for(i in 1:76){
  for(j in (i+1):77){
    if(W0_SZ_l[i,j] == 1){
      con3_SZ_l = c(con3_SZ_l,i,j)
      strength_SZ_l = c(strength_SZ_l,abs(W_SZ_l[i,j])) 
    }
  } 
}


colors0 = c()
for(i in 1:76){
  for(j in (i+1):77){
    if(W0_SZ_l[i,j] == 1){
      colors0 = c(colors0,sign(W_SZ_l[i,j]) ) 
    }
  } 
}

#edges = matrix(con3,ncol=2,byrow=T)
#edges = cbind(edges,strength)

colnames(connection_SZ_l) = names
names2_SZ_l = c()

CID_SZ_l = rep(c(1:8),times=c(CID0_SZ_l[1],CID0_SZ_l[2],CID0_SZ_l[3],CID0_SZ_l[4],CID0_SZ_l[5],CID0_SZ_l[6],CID0_SZ_l[7],sum(CID0_SZ_l[-(1:7)]) ))

for(i in 1:77){
  if (Clist_SZ_l[i]<10){
    names2_SZ_l[i] = paste("  ",Clist_SZ_l[i],": ",names[Clist_SZ_l[i]],sep="")
  } else {
    names2_SZ_l[i] = paste(Clist_SZ_l[i],": ",names[Clist_SZ_l[i]],sep="")  
  }
}

g_SZ_l <- make_empty_graph()
g_SZ_l <- make_graph(edges = con3_SZ_l, n = 77, directed = FALSE)
V(g_SZ_l)$color[Clist_SZ_l] <- col00[CID_SZ_l]
E(g_SZ_l)$width <- strength_SZ_l/mean(strength_SZ_l)
E(g_SZ_l)$color <- colors0
E(g_SZ_l)$color[E(g_SZ_l)$color==1] = brewer.pal(n = 9, name = "Reds")[4]
E(g_SZ_l)$color[E(g_SZ_l)$color==-1] = brewer.pal(n = 9, name = "Blues")[4]



minC_SZ_l <- rep(-Inf, vcount(g_SZ_l))
maxC_SZ_l <- rep(Inf, vcount(g_SZ_l))
minC_SZ_l[1] <- maxC_SZ_l[1] <- 0
location_SZ_l <- layout_with_fr(g_SZ_l, minx=minC_SZ_l, maxx=maxC_SZ_l,miny=minC_SZ_l, maxy=maxC_SZ_l)

plot(g_SZ_l, layout=location0, vertex.size = 7,vertex.label.dist = 0, vertex.label.cex=0.8,edge.width=E(g_SZ_l)$width^3,main="SZ Low-Hallucination")
legend("topleft",legend=names2_SZ_l[1:40],col=col00[CID_SZ_l][1:40],pch=16,cex=0.7)
legend("topright",legend=names2_SZ_l[41:77],col=col00[CID_SZ_l][41:77],pch=16,cex=0.7)


plot(g_SZ_l,  vertex.size = 7,vertex.label.dist = 0, vertex.label.cex=0.8,edge.width=E(g_SZ_l)$width^3,main="SZ Low-Hallucination")
legend("topleft",legend=names2_SZ_l[1:40],col=col00[CID_SZ_l][1:40],pch=16,cex=0.7)
legend("topright",legend=names2_SZ_l[41:77],col=col00[CID_SZ_l][41:77],pch=16,cex=0.7)





########################################################################################
########################################################################################
############################### HC #######################################################
########################################################################################
########################################################################################


connection_HC = matrix(0,77,77)


W0_HC= matrix(0,77,77)
W0_HC[abs(W_HC)>0.34] = 1
sum(W0_HC)/77^2
ind=1


for(i in 1:length(CID0_HC)){
  connection[Clist_HC[ind:(ind+CID0_HC[i]-1)],Clist_HC[ind:(ind+CID0_HC[i]-1)]]=1
  ind = ind + CID0_HC[i]
}


con2_HC = c()
ind=1
for(i in 1:76){
  for(j in (i+1):77){
    if(connection_HC[i,j] == 1){
      con2_HC = c(con2_HC,i,j)
    }
  } 
}

con3_HC = c()
ind=1
strength_HC = c()
for(i in 1:76){
  for(j in (i+1):77){
    if(W0_HC[i,j] == 1){
      con3_HC = c(con3_HC,i,j)
      strength_HC = c(strength_HC,abs(W_HC[i,j]) ) 
    }
  } 
}

colors0 = c()
for(i in 1:76){
  for(j in (i+1):77){
    if(W0_HC[i,j] == 1){
      colors0 = c(colors0,sign(W_HC[i,j]) ) 
    }
  } 
}
#edges = matrix(con3,ncol=2,byrow=T)
#edges = cbind(edges,strength)

colnames(connection_HC) = names

dim(connection)


names2_HC = c()

CID_HC = rep(1:8,times=c(CID0_HC[1],CID0_HC[2],CID0_HC[3],CID0_HC[4],CID0_HC[5],CID0_HC[6],CID0_HC[7],sum(CID0_HC[-(1:7)]) ))

for(i in 1:77){
  if (Clist_HC[i]<10){
    names2_HC[i] = paste("  ",Clist_HC[i],": ",names[Clist_HC[i]],sep="")
  } else {
    names2_HC[i] = paste(Clist_HC[i],": ",names[Clist_HC[i]],sep="")  
  }
}

g_HC <- make_empty_graph()
g_HC <- make_graph(edges = con3_HC, n = 77, directed = FALSE)
V(g_HC)$color[Clist_HC] <- col00[CID_HC]
E(g_HC)$width <- strength_HC/mean(strength_HC)
E(g_HC)$color <- colors0
E(g_HC)$color[E(g_HC)$color==1] = brewer.pal(n = 9, name = "Reds")[4]
E(g_HC)$color[E(g_HC)$color==-1] = brewer.pal(n = 9, name = "Blues")[4]



minC_HC <- rep(-Inf, vcount(g_HC))
maxC_HC <- rep(Inf, vcount(g_HC))
minC_HC[1] <- maxC_HC[1] <- 0
location_HC <- layout_with_fr(g_HC, minx=minC_HC, maxx=maxC_HC,miny=minC_HC, maxy=maxC_HC)

plot(g_HC, layout=location, vertex.size = 7,vertex.label.dist = 0, vertex.label.cex=0.8,edge.width=E(g_HC)$width^3,main="HC")
legend("topleft",legend=names2_HC[1:37],col=col00[CID_HC][1:37],pch=16,cex=0.7)
legend("topright",legend=names2_HC[38:77],col=col00[CID_HC][38:77],pch=16,cex=0.7)

plot(g_HC,  vertex.size = 7,vertex.label.dist = 0, vertex.label.cex=0.8,edge.width=E(g)$width^3,main="HC")
legend("topleft",legend=names2_HC[1:37],col=col00[CID_HC][1:37],pch=16,cex=0.7)
legend("topright",legend=names2_HC[38:77],col=col00[CID_HC][38:77],pch=16,cex=0.7)


########################################################################################
########################################################################################
############################### NCVH #######################################################
########################################################################################
########################################################################################


connection_NCVH = matrix(0,77,77)


W0_NCVH= matrix(0,77,77)
W0_NCVH[abs(W_NCVH)>0.34] = 1
sum(W0_NCVH)/77^2
ind=1


for(i in 1:length(CID0_NCVH)){
  connection[Clist_NCVH[ind:(ind+CID0_NCVH[i]-1)],Clist_NCVH[ind:(ind+CID0_NCVH[i]-1)]]=1
  ind = ind + CID0_NCVH[i]
}


con2_NCVH = c()
ind=1
for(i in 1:76){
  for(j in (i+1):77){
    if(connection_NCVH[i,j] == 1){
      con2_NCVH = c(con2_NCVH,i,j)
    }
  } 
}

con3_NCVH = c()
ind=1
strength_NCVH = c()
for(i in 1:76){
  for(j in (i+1):77){
    if(W0_NCVH[i,j] == 1){
      con3_NCVH = c(con3_NCVH,i,j)
      strength_NCVH = c(strength_NCVH,abs(W_NCVH[i,j]) ) 
    }
  } 
}



colors0 = c()
for(i in 1:76){
  for(j in (i+1):77){
    if(W0_NCVH[i,j] == 1){
      colors0 = c(colors0,sign(W_NCVH[i,j]) ) 
    }
  } 
}


edges_NCVH = matrix(con3_NCVH,ncol=2,byrow=T)
edges_NCVH = cbind(edges_NCVH,strength_NCVH)

#edges = matrix(con3,ncol=2,byrow=T)
#edges = cbind(edges,strength)

colnames(connection_NCVH) = names

names2_NCVH = c()

CID_NCVH = rep(1:8,times=c(CID0_NCVH[1],CID0_NCVH[2],CID0_NCVH[3],CID0_NCVH[4],CID0_NCVH[5],CID0_NCVH[6],CID0_NCVH[7],sum(CID0_NCVH[-(1:7)]) ))

for(i in 1:77){
  if (Clist_NCVH[i]<10){
    names2_NCVH[i] = paste("  ",Clist_NCVH[i],": ",names[Clist_NCVH[i]],sep="")
  } else {
    names2_NCVH[i] = paste(Clist_NCVH[i],": ",names[Clist_NCVH[i]],sep="")  
  }
}

g_NCVH <- make_empty_graph()
g_NCVH <- make_graph(edges = con3_NCVH, n = 77, directed = FALSE)
V(g_NCVH)$color[Clist_NCVH] <- col00[CID_NCVH]
E(g_NCVH)$width <- strength_NCVH/mean(strength_NCVH)
E(g_NCVH)$color <- colors0
E(g_NCVH)$color[E(g_NCVH)$color==1] = brewer.pal(n = 9, name = "Reds")[4]
E(g_NCVH)$color[E(g_NCVH)$color==-1] = brewer.pal(n = 9, name = "Blues")[4]



minC_NCVH <- rep(-Inf, vcount(g_NCVH))
maxC_NCVH <- rep(Inf, vcount(g_NCVH))
minC_NCVH[1] <- maxC_NCVH[1] <- 0
location_NCVH <- layout_with_fr(g_NCVH, minx=minC_NCVH, maxx=maxC_NCVH,miny=minC_NCVH, maxy=maxC_NCVH)

plot(g_NCVH, layout=location, vertex.size = 7,vertex.label.dist = 0, vertex.label.cex=0.8,edge.width=E(g_NCVH)$width^3,main="NCVH")
legend("topleft",legend=names2_NCVH[1:40],col=col00[CID_NCVH][1:40],pch=16,cex=0.7)
legend("topright",legend=names2_NCVH[41:77],col=col00[CID_NCVH][41:77],pch=16,cex=0.7)

plot(g_NCVH,  vertex.size = 7,vertex.label.dist = 0, vertex.label.cex=0.8,edge.width=E(g_NCVH)$width^3,main="NCVH")
legend("topleft",legend=names2_NCVH[1:40],col=col00[CID_NCVH][1:40],pch=16,cex=0.7)
legend("topright",legend=names2_NCVH[41:77],col=col00[CID_NCVH][41:77],pch=16,cex=0.7)




########################################################################################
########################################################################################
############################### High/pos #######################################################
########################################################################################
########################################################################################

var <- readxl::read_xlsx("/Users/Hongju.Park/Dropbox/Sonia_network/CTQ_BPRS_ASI domains scores for SEM.xlsx", sheet = 1,na="ms")
var_SZ0 = var[var$group=="SZ",]
var_SZ= var_SZ[!is.na(var_SZ0$`Positive-BPRS`),]
var_SZ0 = var_SZ0[!is.na(var_SZ0$`Positive-BPRS`),] 
md_pos = median(var_SZ0$`Positive-BPRS`)

index_high_pos_sz = which(var_SZ0$`Positive-BPRS` > md_pos)
index_low_pos_sz = which(var_SZ0$`Positive-BPRS` <= md_pos)

var_high_pos_sz = var_SZ[index_high_pos_sz,]
var_low_pos_sz = var_SZ[index_low_pos_sz,]


W_Hpos = diag(1,77)

for(i in 1:76){
  for(j in (i+1):77){
    W_Hpos[i,j] = cor(na.omit(data.frame(as.numeric(var_high_pos_sz[,i]),as.numeric(var_high_pos_sz[,j]) )))[1,2]
    W_Hpos[j,i] = W_Hpos[i,j] 
  }
}

results_Hpos=dense(abs(W_Hpos),threshold=0.34,lambda=1.15,Advanced=TRUE)
CID0_Hpos= results_Hpos$CID
CID0_Hpos = c(CID0_Hpos[1:7],sum(CID0_Hpos[-(1:7)]))
Clist_Hpos=results_Hpos$Clist
CIDsum_Hpos = cumsum(results_Hpos$CID) 

W0_Hpos= matrix(0,77,77)
W0_Hpos[abs(W_Hpos)>0.34] = 1
sum(W0_Hpos)/77^2

connection_Hpos = matrix(0,77,77)

for(i in 1:length(CID0_Hpos)){
  connection_Hpos[Clist_Hpos[ind:(ind+CID0_Hpos[i]-1)],Clist_Hpos[ind:(ind+CID0_Hpos[i]-1)]]=1
  ind = ind + CID0_Hpos[i]
}


con2_Hpos = c()
ind=1
for(i in 1:76){
  for(j in (i+1):77){
    if(connection_Hpos[i,j] == 1){
      con2_Hpos = c(con2_Hpos,i,j)
    }
  } 
}

colors0 = c()
for(i in 1:76){
  for(j in (i+1):77){
    if(W0_Hpos[i,j] == 1){
      colors0 = c(colors0,sign(W0_Hpos[i,j]) ) 
    }
  } 
}


con3_Hpos = c()
ind=1
strength_Hpos = c()
for(i in 1:76){
  for(j in (i+1):77){
    if(W0_Hpos[i,j] == 1){
      con3_Hpos = c(con3_Hpos,i,j)
      strength_Hpos = c(strength_Hpos,abs(W_Hpos[i,j])) 
    }
  } 
}


edges_Hpos = matrix(con3_Hpos,ncol=2,byrow=T)
edges_Hpos = cbind(edges_Hpos,strength_Hpos)

#edges = matrix(con3,ncol=2,byrow=T)
#edges = cbind(edges,strength)

colnames(connection_Hpos) = names



names2_Hpos = c()

CID_Hpos = rep(1:8,times=c(CID0_Hpos[1],CID0_Hpos[2],CID0_Hpos[3],CID0_Hpos[4],CID0_Hpos[5],CID0_Hpos[6],CID0_Hpos[7],sum(CID0_Hpos[-(1:7)]) ))

for(i in 1:77){
  if (Clist_Hpos[i]<10){
    names2_Hpos[i] = paste("  ",Clist_Hpos[i],": ",names[Clist_Hpos[i]],sep="")
  } else {
    names2_Hpos[i] = paste(Clist_Hpos[i],": ",names[Clist_Hpos[i]],sep="")  
  }
}

g_Hpos <- make_empty_graph()
g_Hpos <- make_graph(edges = con3_Hpos, n = 77, directed = FALSE)
V(g_Hpos)$color[Clist_Hpos] <- col00[CID_Hpos]
E(g_Hpos)$width <- strength_Hpos/mean(strength_Hpos)

E(g_Hpos)$color <- colors0
E(g_Hpos)$color[E(g_Hpos)$color==1] = brewer.pal(n = 9, name = "Reds")[4]
E(g_Hpos)$color[E(g_Hpos)$color==-1] = brewer.pal(n = 9, name = "Blues")[4]



minC_Hpos <- rep(-Inf, vcount(g_Hpos))
maxC_Hpos <- rep(Inf, vcount(g_Hpos))
minC_Hpos[1] <- maxC_Hpos[1] <- 0
location_Hpos <- layout_with_fr(g_Hpos, minx=minC_Hpos, maxx=maxC_Hpos,miny=minC_Hpos, maxy=maxC_Hpos)

plot(g_Hpos, layout=location, vertex.size = 7,vertex.label.dist = 0, vertex.label.cex=0.8,edge.width=E(g_Hpos)$width^3,main="High-positive symptom SZ")
legend("topleft",legend=names2_Hpos[1:38],col=col00[CID_Hpos][1:38],pch=16,cex=0.7)
legend("topright",legend=names2_Hpos[39:77],col=col00[CID_Hpos][39:77],pch=16,cex=0.7)

plot(g_Hpos,  vertex.size = 7,vertex.label.dist = 0, vertex.label.cex=0.8,edge.width=E(g_Hpos)$width^3,main="High-positive symptom SZ")
legend("topleft",legend=names2_Hpos[1:38],col=col00[CID_Hpos][1:38],pch=16,cex=0.7)
legend("topright",legend=names2_Hpos[39:77],col=col00[CID_Hpos][39:77],pch=16,cex=0.7)




########################################################################################
########################################################################################
############################### Low/pos #######################################################
########################################################################################
########################################################################################



W_Lpos = diag(1,77)

for(i in 1:76){
  for(j in (i+1):77){
    W_Lpos[i,j] = cor(na.omit(data.frame(as.numeric(var_low_pos_sz[,i]),as.numeric(var_low_pos_sz[,j]) )))[1,2]
    W_Lpos[j,i] = W_Lpos[i,j] 
  }
}


results_Lpos=dense(abs(W_Lpos),threshold=0.34,lambda=1.15,Advanced=TRUE)
CID0_Lpos= results_Lpos$CID
CID0_Lpos = c(CID0_Lpos[1:7],sum(CID0_Lpos[-(1:7)]))
Clist_Lpos=results_Lpos$Clist
CIDsum_Lpos = cumsum(results_Lpos$CID) 

W0_Lpos= matrix(0,77,77)
W0_Lpos[abs(W_Lpos)>0.34] = 1
sum(W0_Lpos)/77^2


connection_Lpos= matrix(0,77,77)
for(i in 1:length(CID0_Lpos)){
  connection_Lpos[Clist_Lpos[ind:(ind+CID0_Lpos[i]-1)],Clist_Lpos[ind:(ind+CID0_Lpos[i]-1)]]=1
  ind = ind + CID0_Lpos[i]
}


con2_Lpos = c()
ind=1
for(i in 1:76){
  for(j in (i+1):77){
    if(connection_Lpos[i,j] == 1){
      con2_Lpos = c(con2_Lpos,i,j)
    }
  } 
}

con3_Lpos = c()
ind=1
strength_Lpos = c()
for(i in 1:76){
  for(j in (i+1):77){
    if(W0_Lpos[i,j] == 1){
      con3_Lpos = c(con3_Lpos,i,j)
      strength_Lpos = c(strength_Lpos,abs(W_Lpos[i,j])) 
    }
  } 
}

colors0 = c()
for(i in 1:76){
  for(j in (i+1):77){
    if(W0_Lpos[i,j] == 1){
      colors0 = c(colors0,sign(W_Lpos[i,j]) ) 
    }
  } 
}


edges_Lpos = matrix(con3_Lpos,ncol=2,byrow=T)
edges_Lpos = cbind(edges_Lpos,strength_Lpos)

#edges = matrix(con3,ncol=2,byrow=T)
#edges = cbind(edges,strength)

colnames(connection_Lpos) = names



names2_Lpos = c()

CID_Lpos = rep(1:8,times=c(CID0_Lpos[1],CID0_Lpos[2],CID0_Lpos[3],CID0_Lpos[4],CID0_Lpos[5],CID0_Lpos[6],CID0_Lpos[7],sum(CID0_Lpos[-(1:7)]) ))

for(i in 1:77){
  if (Clist_Lpos[i]<10){
    names2_Lpos[i] = paste("  ",Clist_Lpos[i],": ",names[Clist_Lpos[i]],sep="")
  } else {
    names2_Lpos[i] = paste(Clist_Lpos[i],": ",names[Clist_Lpos[i]],sep="")  
  }
}

g_Lpos <- make_empty_graph()
g_Lpos <- make_graph(edges = con3_Lpos, n = 77, directed = FALSE)
V(g_Lpos)$color[Clist_Lpos] <- col00[CID_Lpos]
E(g_Lpos)$width <- strength_Lpos/mean(strength_Lpos)
E(g_Lpos)$color <- colors0
E(g_Lpos)$color[E(g_Lpos)$color==1] = brewer.pal(n = 9, name = "Reds")[4]
E(g_Lpos)$color[E(g_Lpos)$color==-1] = brewer.pal(n = 9, name = "Blues")[4]

str(g_Lpos)
E(g_Lpos)$width

minC_Lpos <- rep(-Inf, vcount(g_Lpos))
maxC_Lpos <- rep(Inf, vcount(g_Lpos))
minC_Lpos[1] <- maxC_Lpos[1] <- 0
location_Lpos <- layout_with_fr(g_Lpos, minx=minC_Lpos, maxx=maxC_Lpos,miny=minC_Lpos, maxy=maxC_Lpos)

plot(g_Lpos, layout=location, vertex.size = 7,vertex.label.dist = 0, vertex.label.cex=0.8,edge.width=E(g_Lpos)$width^3,main="Low-positive symptom SZ")
legend("topleft",legend=names2_Lpos[1:38],col=col00[CID_Lpos][1:38],pch=16,cex=0.7)
legend("topright",legend=names2_Lpos[39:77],col=col00[CID_Lpos][39:77],pch=16,cex=0.7)

plot(g_Lpos,  vertex.size = 7,vertex.label.dist = 0, vertex.label.cex=0.8,edge.width=E(g_Lpos)$width^3,main="Low-positive symptom SZ")
legend("topleft",legend=names2_Lpos[1:38],col=col00[CID_Lpos][1:38],pch=16,cex=0.7)
legend("topright",legend=names2_Lpos[39:77],col=col00[CID_Lpos][39:77],pch=16,cex=0.7)



