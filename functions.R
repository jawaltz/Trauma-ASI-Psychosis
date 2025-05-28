result_process = function(result00,W){
  
  
  xdim = dim(W)[1]
  ydim = dim(W)[2]
  
  lambda_dim=dim(result00)[2]-7
  dim1 = dim(result00)[1]
  result = list()
  for(i in 1:lambda_dim){
    max_ind = which.max(result00[,7+i])
    result[[i]] = list()
    x1 = unique(result00[1:max_ind,6])
    y1 = unique(result00[1:max_ind,7])
    result[[i]]$x = (1:xdim)[-x1]
    result[[i]]$y = (1:ydim)[-y1]
  }
  return(result)
}


kld_cross8 = function(W,index){
  
  if(length(index)==0) {
    result = -Inf
  } else {
    
    
    
    
    p = mean(W)
    W0 = W[index$x,index$y]
    p1 = mean(W0)
    W = as.matrix(W)
    p0 = (sum(W) - sum(W0)) / (dim(W)[1]*dim(W)[2] - length(unlist(W0)) )
    
    
    #n11 = sum(W0)
    #n12 = sum(W0 == 0)
    #n1 = sum(W)
    #n2 = sum(W==0)
    
    const = -length(W0)*(p1*log(p) + (1-p1)*log(1-p)) - (length(W)-length(W0))*(p0*log(p) + (1-p0)*log(1-p))
    
    
    result = 0
    
    if (length(index$x)==0 || length(index$y)==0) {
      result = - Inf
    } else  if (p1 == 1){
      result = length(W0)*(p1*log(p1/p) ) + (length(W)-length(W0))*(p0*log(p0/p) +  (1-p0)*log((1-p0)/(1-p))) 
    } else if (length(index$x)+length(index$y) < dim(W)[1]+dim(W)[2] ){
      result =  length(W0)*(p1*log(p1/p) +  (1-p1)*log((1-p1)/(1-p))) + (length(W)-length(W0))*(p0*log(p0/p) +  (1-p0)*log((1-p0)/(1-p))) 
    } else {
      result = -Inf
    }
    
  }
  
  result/const
}

kld_cross7 = function(W,index){
  
  if(length(index)==0) {
    result = -Inf
  } else {
    
    
    
    
    p = mean(W)
    W0 = W[index$x,index$y]
    p1 = mean(W0)
    W = as.matrix(W)
    p0 = (sum(W) - sum(W0)) / (dim(W)[1]*dim(W)[2] - length(unlist(W0)) )
    
    
    #n11 = sum(W0)
    #n12 = sum(W0 == 0)
    #n1 = sum(W)
    #n2 = sum(W==0)
    
    const = -length(W0)*(p1*log(p1) + (1-p1)*log(1-p1)) - (length(W)-length(W0))*(p0*log(p0) + (1-p0)*log(1-p0))
    
    
    result = 0
    
    if (length(index$x)==0 || length(index$y)==0) {
      result = - Inf
    } else  if (p1 == 1){
      result = length(W0)*(p1*log(p1/p) ) + (length(W)-length(W0))*(p0*log(p0/p) +  (1-p0)*log((1-p0)/(1-p))) 
    } else if (length(index$x)+length(index$y) < dim(W)[1]+dim(W)[2] ){
      result =  length(W0)*(p1*log(p1/p) +  (1-p1)*log((1-p1)/(1-p))) + (length(W)-length(W0))*(p0*log(p0/p) +  (1-p0)*log((1-p0)/(1-p))) 
    } else {
      result = -Inf
    }
    
  }
  
  result/const
}


kld_cross6 = function(W,index){
  
  if(length(index)==0) {
    result = -Inf
  } else {
    
    
    
    
    p = mean(W)
    W0 = W[index$x,index$y]
    p1 = mean(W0)
    W = as.matrix(W)
    p0 = (sum(W) - sum(W0)) / (dim(W)[1]*dim(W)[2] - length(unlist(W0)) )
    
    
    #n11 = sum(W0)
    #n12 = sum(W0 == 0)
    #n1 = sum(W)
    #n2 = sum(W==0)
    
    const = -length(W)*(p*log(p) + (1-p)*log(1-p))
    
    
    result = 0
    
    if (length(index$x)==0 || length(index$y)==0) {
      result = - Inf
    } else  if (p1 == 1){
      result = Inf
    } else if (length(index$x)+length(index$y) < dim(W)[1]+dim(W)[2] ){
      result =  length(W0)*(p*log(p/p1) +  (1-p)*log((1-p)/(1-p1))) + (length(W)-length(W0))*(p*log(p/p0) +  (1-p)*log((1-p)/(1-p0))) 
    } else {
      result = -Inf
    }
    
  }
  
  result/const
}

kld_cross5 = function(W,index){
  
  if(length(index)==0) {
    result = -Inf
  } else {
    
    p = mean(W)
    W0 = W[index$x,index$y]
    p1 = mean(W0)
    W = as.matrix(W)
    p0 = (sum(W) - sum(W0)) / (dim(W)[1]*dim(W)[2] - length(unlist(W0)) )
    
    
    #n11 = sum(W0)
    #n12 = sum(W0 == 0)
    #n1 = sum(W)
    #n2 = sum(W==0)
    
    const = -length(W0)*(p1*log(p1) + (1-p1)*log(1-p1)) - (length(W)-length(W0))*(p0*log(p0) + (1-p0)*log(1-p0))
    
    
    result = 0
    
    if (length(index$x)==0 || length(index$y)==0) {
      result = - Inf
    } else  if (p1 == 1){
      result = length(W0)*(p1*log(p1/p) ) + (length(W)-length(W0))*(p0*log(p0/p) +  (1-p0)*log((1-p0)/(1-p))) 
    } else if (length(index$x)+length(index$y) < dim(W)[1]+dim(W)[2] ){
      result =  length(W0)*(p1*log(p1/p) +  (1-p1)*log((1-p1)/(1-p))) + (length(W)-length(W0))*(p0*log(p0/p) +  (1-p0)*log((1-p0)/(1-p))) 
    } else {
      result = -Inf
    }
    
  }
  
  result
}

kld_cross3 = function(W,index){
  
  if(length(index)==0) {
    result = -Inf
  } else {
    #index = result_00[[1]]
    #Clist=result1$Clist
    #CID=result1$CID
    num_c=length(index)
    
    p = mean(W>0)
    
    
    p1 = c()
    
    W0 = W[index$x,index$y]
    
    
    p1 = mean(unlist(W0))
    W = as.matrix(W)
    p0 = (sum(W>0) - sum(unlist(W0)>0)) / (dim(W)[1]*dim(W)[2] - length(unlist(W0)))
    n11 = sum(W0)
    n12 = sum(W0 == 0)
    n1 = sum(W)
    n2 = sum(W==0)
    
    const = -(n1*p*log(p) + n2*(1-p)*log(1-p))
    
    
    result = 0
    
    if (length(index$x)==0 || length(index$y)==0) {
      result = - Inf
    } else  if (p1 == 1){
      result = n11*log(p1/p)  +   (n1-n11)*p0*log(p0/p)   + (n2-n12)*(1-p0)*log((1-p0)/(1-p))
    } else if (length(index$x)+length(index$y) < dim(W)[1]+dim(W)[2] ){
      result = n11*p1*log(p1/p) +  n12*(1-p1)*log((1-p1)/(1-p) ) + (n1-n11)*p0*log(p0/p) +  (n2-n12)*(1-p0)*log((1-p0)/(1-p))   
    } else {
      result = -Inf
    }
    
  }
  
  result/const
}

kld_cross4 = function(W,index){
  
  if(length(index)==0) {
    result = -Inf
  } else {
    #index = result_00[[1]]
    #Clist=result1$Clist
    #CID=result1$CID
    num_c=length(index)
    
    p = mean(W>0)
    
    
    p1 = c()
    
    W0 = W[index$x,index$y]
    
    
    p1 = mean(unlist(W0))
    W = as.matrix(W)
    p0 = (sum(W>0) - sum(unlist(W0)>0)) / (dim(W)[1]*dim(W)[2] - length(unlist(W0)))
    
    
    n11 = sum(W0)
    n12 = sum(W0 == 0)
    
    n1 = sum(W)
    n2 = sum(W==0)
    
    const = -(n11*p1*log(p1) + n12*(1-p1)*log(1-p1) + + (n1-n11)*p0*log(p0) +  (n2-n12)*(1-p0)*log((1-p0)))
    
    
    result = 0
    
    if (length(index$x)==0 || length(index$y)==0) {
      result = - Inf
    } else  if (p1 == 1){
      result = n11*log(p1/p)  +   (n1-n11)*p0*log(p0/p)   + (n2-n12)*(1-p0)*log((1-p0)/(1-p))
    } else if (length(index$x)+length(index$y) < dim(W)[1]+dim(W)[2] ){
      result = n11*p1*log(p1/p) +  n12*(1-p1)*log((1-p1)/(1-p) ) + (n1-n11)*p0*log(p0/p) +  (n2-n12)*(1-p0)*log((1-p0)/(1-p))   
    } else {
      result = -Inf
    }
    
  }
  
  result/const
}

kld_cross = function(W,index){
  
  if(length(index)==0) {
    result = -Inf
  } else {
    #index = result_00[[1]]
    #Clist=result1$Clist
    #CID=result1$CID
    num_c=length(index)
    
    p = sum(W>0)/(dim(W)[1]*dim(W)[2])
    
    const = -(sum(W>0)*p*log(p) + sum(W==0)*(1-p)*log(1-p))
    
    
    p1 = c()

    W0 = W[index$x,index$y]

    p1 = mean(unlist(W0))
    W = as.matrix(W)
    p0 = (sum(W>0) - sum(unlist(W0)>0)) / (dim(W)[1]*dim(W)[2] - length(unlist(W0)))
    
    result = 0
    
    if (length(index$x)==0 || length(index$y)==0) {
      result = - Inf
    } else  if (p1 == 1){
      result = sum(unlist(W0)>0)*log(p1/p)  +   (sum(W)-sum(unlist(W0)))*p0*log(p0/p)   +  ( length(which(W==0)) - length(which(unlist(W0)==0)) )*(1-p0)*log((1-p0)/p)
    } else if (length(index$x)+length(index$y) < dim(W)[1]+dim(W)[2] ){
      result = sum(unlist(W0)>0)*p1*log(p1/p) +  (length(unlist(W0))-sum(unlist(W0)) )*(1-p1)*log((1-p1)/p) + (sum(W)-sum(unlist(W0)))*p0*log(p0/p) +  ( length(which(W==0)) - length(which(unlist(W0)==0)) )*(1-p0)*log((1-p0)/p)   
    } else {
      result = -Inf
    }
    
  }
  
  result/const
}

kld_cross2 = function(W,index){
  
  if(length(index)==0) {
    result = -Inf
  } else {
    #index = result_00[[1]]
    #Clist=result1$Clist
    #CID=result1$CID
    num_c=length(index)
    
    p = sum(W>0)/(dim(W)[1]*dim(W)[2])
    
    W0 = list()
    p1 = c()
    for(i in 1:num_c){
      W0[[i]] = W[index[[i]]$x,index[[i]]$y]
    }
    
    p1 = c()
    
    for(i in 1:num_c){
      p1[i] = mean(W0[[i]])
    }
    
    W = as.matrix(W)
    p0 = (sum(W>0) - sum(unlist(W0)>0)) / (dim(W)[1]*dim(W)[2] - length(unlist(W0)))
    
    result = (sum(W)-sum(unlist(W0)))*p0*log(p0/p) +  ( length(which(W==0)) - length(which(unlist(W0)==0)) )*(1-p0)*log((1-p0)/p) 
    
    
    for(i in 1:num_c){
      if (length(index[[1]]$x)==0 || length(index[[1]]$y)==0) {
        result = - Inf
      } else  if (p1[i] == 1){
        result = result + sum(unlist(W0[[i]])>0)*log(p1[i]/p)       
      } else if (length(index[[i]]$x)+length(index[[i]]$y) < dim(W)[1]+dim(W)[2] ){
        result = result + sum(unlist(W0[[i]])>0)*p1[i]*log(p1[i]/p) +  (length(unlist(W0[[i]]))-sum(unlist(W0[[i]])) )*(1-p1[i])*log((1-p1[i])/p) 
      } else {
        result = -Inf
      }
    }
    
  }
  
  result
}


kld = function(W,Clist,CID,num_c){
  
  #Clist=results[[11]]$Clist
  #CID=results[[11]]$CID
  #num_c=length(results[[11]]$CID)
  
  
  p = sum(W>0)/(dim(W)[1]*dim(W)[2])
  index = list()
  cum_CID = cumsum(CID)
  initial = 1
  
  count_cluster = 0
  
  
  for (i in 1:num_c){
    index[[i]] = Clist[initial:cum_CID[i]]
    initial =  cum_CID[i] + 1
    count_cluster = count_cluster+sum(W[index[[i]],index[[i]]])
  }
  
  p1 = count_cluster/(sum(CID^2)-sum(CID))
  p0 = (sum(W>0) -count_cluster) / (dim(W)[1]*(dim(W)[1]-1) -sum(CID^2) + sum(CID))
  
  if (p1 == 1){
    result = count_cluster*log(p1/p)  + (sum(W>0)-count_cluster)*p0*log(p0/p) + (dim(W)[1]*(dim(W)[1]-1)-sum(CID^2)-(sum(W>0)-count_cluster))*(1-p0)*log((1-p0)/p)
  } else if (CID[1] < length(Clist)){
    result = count_cluster*p1*log(p1/p) + (sum(CID^2)-count_cluster)*(1-p1)*log((1-p1)/p) + (sum(W>0)-count_cluster)*p0*log(p0/p) + (dim(W)[1]*(dim(W)[1]-1)-(sum(CID^2)-sum(CID)) -(sum(W>0)-count_cluster))*(1-p0)*log((1-p0)/p)   
  } else {
    result = -Inf
  }
  result
}

greedy_bipar3 <- function(W, lambda,  index) {
  
  
  W1 <- W
  #W <- W0
  #dST_vec <- numeric(length(c_vec))
  
  #W <- W00
  
  deleteS_list <- list()
  deleteT_list <- list()
  density_list <- matrix(0, nrow = (nrow(W) + ncol(W) - 1), ncol = 8)
  
  #c = 1
  nr = nrow(W)
  nc = ncol(W)
  row <- nr
  col <- nc
  R <- rowSums(W)
  C <- colSums(W)
  r_which_0 = integer(0)
  c_which_0 = integer(0)
  
  for (i in 1:(nr + nc - 2)) {
    
    if (length(r_which_0) > nr-2 ||length(c_which_0) > nc -2) break
    
    row <- nr - length(r_which_0)
    col <- nc - length(c_which_0)
    #k1 <- max(1,floor(row / k))
    #k2 <- max(1,floor(col / k))
    k1=1
    k2=1
    IndS <- order(R)[(length(r_which_0) + 1):(length(r_which_0) + k1)]
    resultR <- R[IndS]
    dS <- mean(resultR)
    IndT <- order(C)[(length(c_which_0) + 1):(length(c_which_0) + k2)]
    resultC <- C[IndT]
    dT <- mean(resultC)
    
    c = col/row
    
    if (c * dS <= dT ) {
      if(length(IndS)>1){
        C <- C - colSums(W[IndS, ])
      } else {
        C <- C - W[IndS, ]
      }
      r_which_0 = c(r_which_0,IndS)
      row <- nr - length(r_which_0)
      W[IndS, ] <- 0
      R[IndS] <- 0
      density_list[i, 1] <- sum(R != 0)
      density_list[i, 2] <- sum(C != 0)
      density_list[i, 3] <- 1
      
    } else {
      if(length(IndT)>1){
        R <- R - rowSums(W[, IndT])
      } else {
        R <- R - W[, IndT]
      }
      c_which_0 = c(c_which_0, IndT)
      col <- nc - length(c_which_0)
      W[, IndT] <- 0
      C[IndT] <- 0
      density_list[i, 1] <- row
      density_list[i, 2] <- col
      density_list[i, 3] <- 2
      
    }
    
    deleteS_list[[i]] <- IndS
    deleteT_list[[i]] <- IndT
    density_list[i, 8] <- (sum(R)) / (density_list[i, 1] * density_list[i, 2])^lambda
    
    #cat("i: ",i,"~time0:", time1-time0,"~time1:", time2-time1,"~time2:", time2-time1, "\n")
    #cat("i: ",i, "\n")
  }
  
  dST <- max(density_list[, 8])
  indST <- which.max(density_list[, 8])
  remove_s <- integer(0)
  remove_t <- integer(0)
  
  
  
  if (indST>0){
    for (j in 1:indST) {
      if (density_list[j, 3] == 1) {
        remove_s <- union(remove_s, deleteS_list[[j]])
      } else {
        remove_t <- union(remove_t, deleteT_list[[j]])   
      }
    }
  }
  
  
  s_out0 <- remove_s
  t_out0 <- remove_t
  #t_out <- sort(which(colSums(W) == 0))
  #s_out <- which(rowSums(W) == 0)
  
  if (sum(W1) / (dim(W1)[1] * dim(W1)[2])^lambda >  dST){
    return(list(s_out = c(), t_out = c()))
  } else{
    return(list(s_out = index$x[s_out0], t_out = index$y[t_out0]))
  }  
}


subg1 <- function(W, lambda) {
  
  #W1 <- W00
  #W = W00
  #W = abs(W00)
  #W1 <- W
  #W1 = abs(W1)

  W <- abs(W00)
  
  deleteS_list <- list()
  deleteT_list <- list()
  density_list <- matrix(0, nrow = (nrow(W) + ncol(W) - 1), ncol = 7 + length(lambda))
  nr = nrow(W)
  nc = ncol(W)
  row <- nr
  col <- nc
  R <- rowSums(W)
  C <- colSums(W)
  r_which_0 = integer(0)
  c_which_0 = integer(0)
  row <- nr - length(r_which_0)
  col <- nc - length(c_which_0)
  
  for (i in 1:(nr + nc - 2)) {
    
    if (length(r_which_0) > nr-2 ||length(c_which_0) > nc -2) break
    
    IndS <- order(R)[length(r_which_0) + 1]
    resultR <- R[IndS]
    dS <- mean(resultR)
    IndT <- order(C)[length(c_which_0) + 1]
    resultC <- C[IndT]
    dT <- mean(resultC)
    
    c = col/row
    
    if (c * dS <= dT ) {
      if(length(IndS)>1){
        C <- C - colSums(W[IndS, ])
      } else {
        C <- C - W[IndS, ]
      }
      r_which_0 = c(r_which_0,IndS)
      row <- nr - length(r_which_0)
      W[IndS, ] <- 0
      R[IndS] <- 0
      density_list[i, 1] <- sum(R != 0)
      density_list[i, 2] <- sum(C != 0)
      density_list[i, 3] <- 1
      
    } else {
      if(length(IndT)>1){
        R <- R - rowSums(W[, IndT])
      } else {
        R <- R - W[, IndT]
      }
      c_which_0 = c(c_which_0, IndT)
      col <- nc - length(c_which_0)
      W[, IndT] <- 0
      C[IndT] <- 0
      density_list[i, 1] <- row
      density_list[i, 2] <- col
      density_list[i, 3] <- 2
      
    }
    
    density_list[i, 6] <- IndS
    density_list[i, 7] <- IndT
    density_list[i, 8:(7+length(lambda))] <- (sum(R)) / (density_list[i, 1] * density_list[i, 2])^lambda
    
    
    #cat("i: ",i,"~time0:", time1-time0,"~time1:", time2-time1,"~time2:", time2-time1, "\n")
    #cat("i: ",i, "\n")
  }
  
  dST <- apply(density_list[, 8:(7+length(lambda))],2,max)
  indST <- apply(density_list[, 8:(7+length(lambda))],2,which.max)
  
  nodes_out <- list()
  
  for(i in 1:length(lambda)){
    nodes_out[[i]] = list()
    if(indST[i]>1){
      nodes_out[[i]]$x = sort(unique(density_list[indST[i]:dim(density_list)[1], 6]))[-1]
    } else {
      nodes_out[[i]]$x = c()
    }
    
    if(indST[i]>1){
      nodes_out[[i]]$y = sort(unique(density_list[indST[i]:dim(density_list)[1], 7]))[-1]
    } else {
      nodes_out[[i]]$y = c()
    }
    
  }
  
  
  return(nodes_out)
  
}




subg2 <- function(W, lambda) {
  
  #W1 <- W00
  #W = W00
  #W = abs(W00)
  #W1 <- W
  #W1 = abs(W1)
  
  #W <- abs(W00)
  
  deleteS_list <- list()
  deleteT_list <- list()
  density_list <- matrix(0, nrow = (nrow(W) + ncol(W) - 1), ncol = 7 + length(lambda))
  nr = nrow(W)
  nc = ncol(W)
  row <- nr
  col <- nc
  R <- rowSums(W)
  C <- colSums(W)
  r_which_0 = integer(0)
  c_which_0 = integer(0)
  
  for (i in 1:(nr + nc - 2)) {
    
    if (length(r_which_0) > nr-2 ||length(c_which_0) > nc -2) break

    row <- nr - length(r_which_0)
    col <- nc - length(c_which_0)
    
    if (length(r_which_0)>0){
      IndS <- ((1:nr)[-r_which_0])[which.min(R[(1:nr)[-r_which_0]])] 
    } else {
      IndS <- which.min(R)
    }
    
    resultR <- R[IndS]
    dS <- mean(resultR)
    
    if (length(c_which_0)>0){
      IndT <- ((1:nc)[-c_which_0])[which.min(C[(1:nc)[-c_which_0]])] 
    } else {
      IndT <- which.min(C)
    }
    
    resultC <- C[IndT]
    dT <- mean(resultC)
    
    c = col/row
    
    if (c * dS <= dT ) {
      if(length(IndS)>1){
        C <- C - colSums(W[IndS, ])
      } else {
        C <- C - W[IndS, ]
      }
      r_which_0 = c(r_which_0,IndS)
      row <- nr - length(r_which_0)
      W[IndS, ] <- 0
      R[IndS] <- 0
      density_list[i, 1] <- row
      density_list[i, 2] <- col
      density_list[i, 3] <- 1
      density_list[i, 6] <- IndS
    } else {
      if(length(IndT)>1){
        R <- R - rowSums(W[, IndT])
      } else {
        R <- R - W[, IndT]
      }
      c_which_0 = c(c_which_0, IndT)
      col <- nc - length(c_which_0)
      W[, IndT] <- 0
      C[IndT] <- 0
      density_list[i, 1] <- row
      density_list[i, 2] <- col
      density_list[i, 3] <- 2
      density_list[i, 7] <- IndT
    }


    
    density_list[i, 8:(7+length(lambda))] <- (sum(R)) / (density_list[i, 1] * density_list[i, 2])^lambda
    
    density_list[(i-1):(i+1),]
    #cat("i: ",i,"~time0:", time1-time0,"~time1:", time2-time1,"~time2:", time2-time1, "\n")
    #cat("i: ",i, "\n")
  }
  
  
  dST <- apply(density_list[, 8:(7+length(lambda))],2,max)
  indST <- apply(density_list[, 8:(7+length(lambda))],2,which.max)
  
  return(density_list)
  
}


subg = function(W, lambda=0.5, num_subg){
  start.time <- Sys.time()
  W0=W
  result = list()
  i = 1
  subg0 = list()
  index = list()
  index[[i]] = list(x=1:dim(W0)[1],y=1:dim(W0)[2])
  result[[i]] = greedy_bipar3(W0, lambda=lambda, index=index[[i]])
  
  while (i <= num_subg && length(result[[i]]$s_out)>1 && length(result[[i]]$t_out)>1 && length(result[[i]]$s_out) < length(index[[i]]$x) && length(result[[i]]$t_out) < length(index[[i]]$y)){
    subg0[[i]] = list(x=index[[i]]$x[ !index[[i]]$x %in% result[[i]]$s_out],y=index[[i]]$y[!index[[i]]$y %in% result[[i]]$t_out])
    index[[i+1]] = list(x=result[[i]]$s_out,y=result[[i]]$t_out)
    W11=W0[index[[i+1]]$x,index[[i+1]]$y]
    if(sum(W11)==0){
      break
    }
    result[[i+1]] = greedy_bipar3(W=W11, lambda=lambda, index=index[[i+1]])
    if ( length(index[[i+1]]$x[ !index[[i+1]]$x %in% result[[i+1]]$s_out]) <= 2 && length(index[[i]]$y[!index[[i]]$y %in% result[[i]]$t_out])<=2 ) {
      break
    }
    
    cat("lambda: ",lambda, " sub g i",i,"length x: ",length(result[[i+1]]$s_out),"length y: ",length(result[[i+1]]$t_out) ,"\n")
    i = i+1
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat("time: ", time.taken, "\n")
  results = subg0
  return(results)
}

greedy_peeling <- function(W_greedy, lambda, Node_Seq, removing_node, Advanced){
  # Number of nodes in the graph
  
  
  #Node_Seq = 1:3
  #removing_node= c()
  #W_greedy = W
  
  
  N <- length(Node_Seq)
  # N <- sqrt(length(W_greedy))
  Recording_Matrix <- matrix(0, nrow = N, ncol = 2)
  Recording_Clist <- 1:N
  W_temp <- W_greedy[1:N,1:N]
  
  C <- colSums(W_temp)
  remove_idx <- integer(0)
  ite <- 0
  
  for (i in seq(N, 1, -1)) {
    
    # Calculate sum ignoring removed indices
    
    C[remove_idx] <- Inf
    idx_min_temp <- which.min(C)
    # Update remove_idx
    remove_idx <- c(remove_idx, idx_min_temp)
    # Update sum
    C <- C - W_temp[idx_min_temp, ]
    C[remove_idx] <- 0
    sum_W_temp <- sum(C)
    
    # Calculate score_temp
    ite <- ite + 1
    score_temp <- sum_W_temp / ((N - ite) ^ lambda)
    Recording_Matrix[N - i + 1, ] <- c(Recording_Clist[idx_min_temp], score_temp)
    
  }
  
  
  # Find the index of the maximum score
  max_idx <- which.max(Recording_Matrix[, 2])
  
  if (Advanced == TRUE) {
    
    # Extract the removing nodes and node sequence
    Node_Seq_temp = Recording_Matrix[N:(max_idx + 1), 1]
    removing_node_temp = Recording_Matrix[1:max_idx, 1]
    temp_index = c(Node_Seq_temp,removing_node_temp)
    
    split=sum(W_greedy[Node_Seq_temp,Node_Seq_temp])/length(Node_Seq_temp)^lambda + sum(W_greedy[removing_node_temp,removing_node_temp])/length(removing_node_temp)^lambda
    
    total=sum(W_greedy)/N^lambda
    
    if (split > total){
      # Concatenate node sequence and removing nodes
      removing_node <- c(removing_node,Node_Seq[Recording_Matrix[1:max_idx,1]])
      Node_Seq <- Node_Seq[Recording_Matrix[N:(max_idx + 1), 1]]
      
      Clist <- c(Node_Seq, removing_node)
      # Extract subgraph based on Clist
      N2 = nrow(W_greedy)
      if (N2 == length(temp_index) ){
        W_interim <- W_greedy[temp_index, temp_index]
      } else {
        W_interim <- W_greedy[c(temp_index,(length(temp_index)+1):N2), c(temp_index,(length(temp_index)+1):N2)]
      }
      greedy_peeling(W_interim,lambda,Node_Seq,removing_node,Advanced)
    } else {
      return(list(W_interim = W_greedy,
                  Clist = c(Node_Seq,removing_node),
                  Node_Seq = Node_Seq,
                  removing_node = removing_node))  
    }
    
  } else {
    # Extract the removing nodes and node sequence
    removing_node <- Recording_Matrix[1:max_idx, 1]
    if (length(removing_node) != 1) {
      Node_Seq <- Recording_Matrix[N:(max_idx + 1), 1]
    } else {
      Node_Seq <- Recording_Matrix[, 1]
    }
    # Concatenate node sequence and removing nodes
    Clist <- c(Node_Seq, removing_node)
    # Extract subgraph based on Clist
    W_interim <- W_greedy[Clist, Clist]
    
    return(list(W_interim = W_interim,
                Clist = Clist,
                Node_Seq = Node_Seq,
                removing_node = removing_node))
  }
  
}



subg3 = function(W, lambda=0.5,lambda1,lambda2, xn, yn, num_subg){
  
  start.time <- Sys.time()
  W0=W
  result = list()
  i = 1
  subg0 = list()
  index = list()
  index[[i]] = list(x=1:dim(W0)[1],y=1:dim(W0)[2])
  result[[i]] = greedy_bipar5(W0, lambda=lambda, lambda1,lambda2, index=index[[i]],xn,yn)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat("time: ", time.taken, "\n")
  
  results= list()
  results[[1]] = list(x=index[[i]]$x[ !index[[i]]$x %in% result[[i]]$s_out],y=index[[i]]$y[!index[[i]]$y %in% result[[i]]$t_out])
  return(results)
}



changeind=function(x,n){
  vec1 = unlist(lapply(1:n, function(i) i*(n-(i+1)/2) ))
  i = which.max(x <= vec1)
  if(i==1){
    j = x + i
  } else{
    j = x-vec1[i-1]+i
  } 
  return(c(i,j))
}





greedy_bipar5 <- function(W, lambda,lambda1,lambda2,  index, xn, yn) {
  
  #lambda=0.8
  #lambda1=2
  #lambda2=2
  #xn=20
  #yn=25
  #W <- W00
  
  W1 <- W
  deleteS_list <- list()
  deleteT_list <- list()
  density_list <- matrix(0, nrow = (nrow(W) + ncol(W) - 1), ncol = 11)
  
  #c = 1
  nr = nrow(W)
  nc = ncol(W)
  row <- nr
  col <- nc
  R <- rowSums(W)
  C <- colSums(W)
  r_which_0 = c()
  c_which_0 = c()
  
  sum0 = sum(R)
  c1 = (length(R)*length(C))^lambda/sum0
  c2 = xn^lambda1/length(R)
  c3 = yn^lambda2/length(C)
  
  
  x_ind = 1:(xn*(xn-1)/2)
  y_ind = 1:(yn*(yn-1)/2)
  
  tx = table(unlist(lapply(x_ind,function(x) changeind(x,xn))))
  ty = table(unlist(lapply(y_ind,function(x) changeind(x,yn))))
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  for (i in 1:(nr + nc - 2)) {
    #i=1
    #i=i+1
    if (length(r_which_0) > nr-2 ||length(c_which_0) > nc -2) break
    
    row <- nr - length(r_which_0)
    col <- nc - length(c_which_0)
    
    if(length(r_which_0)==0){
      IndS <- (1:length(R))[which.min(R[(1:length(R))])]
    } else{
      IndS <- (1:length(R))[-r_which_0][which.min(R[(1:length(R))[-r_which_0]])]
      resultR <- R[IndS]
      
    }
    if (length(c_which_0)==0){
      IndT <- (1:length(C))[which.min(C[(1:length(C))])]  
    } else {
      IndT <- (1:length(C))[-c_which_0][which.min(C[(1:length(C))[-c_which_0]])]
    }
    resultC <- C[IndT]  
    resultR <- R[IndS]
    dS <- mean(resultR)
    dT <- mean(resultC)
    
    c = col/row
    
    if (c * dS <= dT ) {
      C <- C - W[IndS, ]
      r_which_0 = c(r_which_0,IndS)
      row <- nr - length(r_which_0)
      W[IndS, ] <- 0
      sum0 = sum0-R[IndS] 
      R[IndS] <- 0
      density_list[i, 1] <- row
      density_list[i, 2] <- col
      density_list[i, 3] <- 1
      density_list[i, 4] <- IndS
      ind_s = changeind(IndS,xn)
      tx[as.numeric(names(tx)) %in% ind_s ] = tx[as.numeric(names(tx)) %in% ind_s ]-1
      if (prod(tx[as.numeric(names(tx)) %in% ind_s ])==0) tx = tx[tx=!0]
    } else {
      R <- R - W[, IndT]
      c_which_0 = c(c_which_0, IndT)
      col <- nc - length(c_which_0)
      W[, IndT] <- 0
      sum0 = sum0-C[IndT] 
      C[IndT] <- 0
      density_list[i, 1] <- row
      density_list[i, 2] <- col
      density_list[i, 3] <- 2
      density_list[i, 5] <- IndT
      ind_t = changeind(IndT,yn)
      ty[as.numeric(names(ty)) %in% ind_t ] = ty[as.numeric(names(ty)) %in% ind_t ]-1
      if( prod(ty[as.numeric(names(ty)) %in% ind_t ]) == 0) ty = ty[ty=!0]
    }
    
    
    deleteS_list[[i]] <- IndS
    deleteT_list[[i]] <- IndT
    
    density_list[i, 8] <- c1*(sum0) / (density_list[i, 1] * density_list[i, 2])^lambda + c2*(sum(tx)/2)/length(tx)^lambda1 + c3*(sum(ty)/2)/length(ty)^lambda2
    density_list[i, 9] <- c1*(sum0) / (density_list[i, 1] * density_list[i, 2])^lambda 
    density_list[i, 10] <- c2*(sum(tx)/2)/length(tx)^lambda1
    density_list[i, 11] <- c3*(sum(ty)/2)/length(ty)^lambda2
    
    
    #cat("i: ", i, "sum: ",sum0, "\n")
    #cat("i: ",i, "\n")
  }
  
  #density_list[401:489,]
  #dim(density_list)
  
  dST <- max(density_list[, 8])
  indST <- which.max(density_list[, 8])
  remove_s <- integer(0)
  remove_t <- integer(0)
  
  if (indST>0){
    for (j in 1:indST) {
      if (density_list[j, 3] == 1) {
        remove_s <- union(remove_s, deleteS_list[[j]])
      } else {
        remove_t <- union(remove_t, deleteT_list[[j]])   
      }
    }
  }
  
  
  s_out0 <- remove_s
  t_out0 <- remove_t
  #t_out <- sort(which(colSums(W) == 0))
  #s_out <- which(rowSums(W) == 0)
  
  #(1:(xn*(xn-1)/2))[-sort(s_out0)]
  #(1:(yn*(yn-1)/2))[-sort(t_out0)]
  
  return(list(s_out = index$x[s_out0], t_out = index$y[t_out0]))
  
}



#' Adaptive dense subgraph extraction
#'
#' Extract dense subgraphs from the entire graph using the greedy algorithm with \eqn{l_0}{l_0} shrinkage.
#'
#' @param W_original Input adjacency matrix
#' @param threshold Threshold value for filtering edges, default to 0.5
#' @param lambda Tuning parameter \eqn{\lambda}{lambda} for greedy peeling algorithm, default to 0.5
#' @return A list containing:\tabular{ll}{
#'    \code{W_dense} \tab Reordered adjacency matrix \cr
#'    \tab \cr
#'    \code{Clist} \tab Node index of the reordered adjacency matrix with detected dense subgraphs\cr
#'    \tab \cr
#'    \code{CID} \tab Number of nodes in each dense subgraph \cr
#' }
#' @references
#' Chen, S., Zhang, Y., Wu, Q., Bi, C., Kochunov, P., & Hong, L. E. (2024).
#' Identifying covariate-related subnetworks for whole-brain connectome analysis.
#' Biostatistics (Oxford, England), 25(2), 541–558. https://doi.org/10.1093/biostatistics/kxad007
#'
#' Wu, Q., Huang, X., Culbreth, A. J., Waltz, J. A., Hong, L. E., & Chen, S. (2022). Extracting brain disease‐related connectome subgraphs by adaptive dense subgraph discovery.
#' Biometrics, 78(4), 1566–1578. https://doi.org/10.1111/biom.13537
#' @export
#' @examples
#' data(sim)
#' matrix <- cor(sim)
#' dense(matrix, 0.6, 0.5)

dense <- function(W_original, threshold = 0.3, lambda = 1, Advanced = FALSE){
  
  #W_original= W
  #lambda=1
  
  if (!all(diag(W_original) == 0)) {
    W_original <- W_original - diag(diag(W_original))
  }
  #initialize
  W_greedy <- W_original
  W_greedy[W_greedy < threshold] <- 0
  
  #store all node lists from each cluster
  Clist <- numeric()
  CID <- numeric()
  orig_Node <- 1:length(W_original)
  
  iter=0
  
  index = 1:nrow(W_greedy)
  removing_node = c()
  
  
  while (length(Clist) < ncol(W_original) - 1) {
    # perform greedy peeling to detect a dense subgraph, remove it, perform on the remaining nodes
    
    index = 1:nrow(W_greedy)
    removing_node = c()
    result <- greedy_peeling(W_greedy, lambda, index, removing_node, Advanced)
    
    index = result$Node_Seq
    removing_node = result$removing_node
    
    Clist_temp <- result$Clist
    Node_Seq <- result$Node_Seq
    remaining_node <- result$removing_node
    
    # record nodes
    Clist <- c(Clist, orig_Node[Clist_temp[1:length(Node_Seq)]])
    CID <- c(CID, length(Node_Seq))
    
    if (length(remaining_node) == 1) break
    # updates
    orig_Node <- orig_Node[remaining_node]
    W_greedy <- W_greedy[remaining_node, remaining_node]
    iter = iter+1
  }
  
  W_dense <- W_original[Clist, Clist]
  
  return(list(W_dense = W_dense, Clist = Clist, CID = CID))
}


cov1 = function(result,sample_cov){
  covmat = matrix(0,dim(sample_cov)[1],dim(sample_cov)[2])
  
  for(i in 1:length(result)){
    covmat[result[[i]][[1]],result[[i]][[2]] ]  = mean(sample_cov[result[[i]][[1]],result[[i]][[2]] ] )
  }
  covmat
}


#install.packages("sRDA")
#library(sRDA)
#help(get_cross_validated_penalty_parameters)

datagen = function (nr_LVs = 1, n = 50, nr_correlated_Xs = c(5), nr_uncorrelated_Xs = 250, cor=c(0.5), 
                    nr_correlated_Ys = c(5), nr_uncorrelated_Ys = 350) 
{
  
  etasd= sqrt(cor)
  sd1 = sqrt(1-cor)
  
  number_of_ksi <- nr_LVs
  number_of_patients <- n
  number_of_Xs_associated_with_ksis <- nr_correlated_Xs
  number_of_not_associated_Xs <- nr_uncorrelated_Xs
  number_of_Ys_associated_with_ksis <- nr_correlated_Ys
  number_of_not_associated_Ys <- nr_uncorrelated_Ys
  
  sig_mat = matrix(0,number_of_ksi,number_of_ksi)
  
  for (i in 1:dim(sig_mat)[1]){
    sig_mat[i,i] = etasd[i]^2
  }
  
  ksi = rmvnorm(number_of_patients, mean = rep(0, number_of_ksi), 
                sigma =sig_mat )
  
  x = matrix(NA, nrow = number_of_patients, ncol = (sum(number_of_Xs_associated_with_ksis) + 
                                                      number_of_not_associated_Xs))
  columncount = 0
  for (i in 1:length(number_of_Xs_associated_with_ksis)) {
    for (j in 1:number_of_Xs_associated_with_ksis[i]) {
      columncount = columncount + 1
      x[, columncount] = rnorm(number_of_patients, mean=ksi[, i], sd=sd1[i])
    }
  }
  for (j in 1:number_of_not_associated_Xs) {
    columncount = columncount + 1
    x[, columncount] = rnorm(number_of_patients, mean=0, sd=1)
  }
  y = matrix(NA, nrow = number_of_patients, ncol = (sum(number_of_Ys_associated_with_ksis) + 
                                                      number_of_not_associated_Ys))
  columncount = 0
  for (i in 1:length(number_of_Ys_associated_with_ksis)) {
    for (j in 1:number_of_Ys_associated_with_ksis[i]) {
      columncount = columncount + 1
      y[, columncount] = rnorm(number_of_patients, mean = ksi[, i], sd = sd1[i])
    }
  }
  for (j in 1:number_of_not_associated_Ys) {
    columncount = columncount + 1
    y[, columncount] = rnorm(number_of_patients, mean=0, sd=1)
  }
  X <- x
  Y <- y
  data_info <- data.frame(number_of_ksi, number_of_patients, 
                          number_of_Xs_associated_with_ksis, number_of_not_associated_Xs, 
                          number_of_Ys_associated_with_ksis, number_of_not_associated_Ys)
  result <- list(X, Y, data_info)
  names(result) <- c("X", "Y", "data_info")
  result
  
}




linedf = function(cluster){
  x = c()
  y = c()
  xend = c()
  yend = c()
  for(i in 1:(length(cluster)-1)){
    x[4*(i-1)+1:4] = c(cluster[i] - 0.5, cluster[i+1] - 0.5, cluster[i+1] - 0.5, cluster[i] - 0.5)
    y[4*(i-1)+1:4] = c(cluster[i] - 0.5, cluster[i] - 0.5, cluster[i+1] - 0.5, cluster[i+1] - 0.5)
    xend[4*(i-1)+1:4] = c(cluster[i+1] - 0.5,cluster[i+1] - 0.5,cluster[i] - 0.5,cluster[i] - 0.5)
    yend[4*(i-1)+1:4] = c(cluster[i] - 0.5,cluster[i+1] - 0.5,cluster[i+1] - 0.5,cluster[i] - 0.5)
  }
  data.frame(x,y,xend,yend)
}




linedf_cross = function(result,number){
  clusterx = 0
  clustery = 0
  for(i in 1:number){
    clusterx = c(clusterx,clusterx[length(clusterx)] + length(result[[i]]$x) )
    clustery = c(clustery,clustery[length(clustery)] + length(result[[i]]$y) )
  }
  x = c()
  y = c()
  xend = c()
  yend = c()
  for(i in 1:number){
    x[4*(i-1)+1:4] = c(clusterx[i] + 0.5, clusterx[i+1] + 0.5, clusterx[i+1] + 0.5, clusterx[i] + 0.5)
    y[4*(i-1)+1:4] = c(clustery[i] + 0.5, clustery[i] + 0.5, clustery[i+1] + 0.5, clustery[i+1] + 0.5)
    xend[4*(i-1)+1:4] = c(clusterx[i+1] + 0.5,clusterx[i+1] + 0.5,clusterx[i] + 0.5,clusterx[i] + 0.5)
    yend[4*(i-1)+1:4] = c(clustery[i] + 0.5,clustery[i+1] + 0.5,clustery[i+1] + 0.5,clustery[i] + 0.5)
  }
  data.frame(x,y,xend,yend)
}

reorder = function(result,mat){
  for(i in 1:length(result)){
    if(length(result[[i]]$y)>1 && length(result[[i]]$x)>1 ) {
      mat_temp = mat[result[[i]]$x,result[[i]]$y]
      result[[i]]$x = result[[i]]$x[order(rowSums(mat_temp),decreasing=TRUE)]
      result[[i]]$y = result[[i]]$y[order(colSums(mat_temp),decreasing=TRUE)]
    }
  }
  result
}


can_cor = function(xlist,ylist,X,Y){
  Sx = cor(X)[unlist(xlist),unlist(xlist)]
  Sy = cor(Y)[unlist(ylist),unlist(ylist)]
  Sxy = cor(X,Y)[unlist(xlist),unlist(ylist)]
  svd0 =svd(Sxy)
  a3 = svd0$u[,1]
  b3 =  svd0$v[,1]
  (t(a3) %*%  Sxy %*% b3 ) / sqrt(as.numeric(t(a3) %*% Sx %*% a3) * as.numeric(t(b3)  %*% Sy %*% b3))
}




