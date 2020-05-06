rm(list=ls())

to_one <- function(Ori_matrix){
  Ori_max <- max(Ori_matrix)
  Ori_min <- min(Ori_matrix)
  for(a in 1:dim(Ori_matrix)[1]){
    for(b in 1: dim(Ori_matrix)[2]){
      Ori_matrix[a,b] <- (Ori_matrix[a,b]-Ori_min)/(Ori_max-Ori_min)
    }
  }
  return(Ori_matrix)
}


to_laplacian_one <- function(Ori_matr){
  D <- diag(rowSums(Ori_matr))
  L <- Ori_matr
  for(a in 1:dim(L)[1]){
    for(b in 1:dim(L)[2]){                
      # L[a,b] <- L[a,b]/((D[a,a]*D[b,b])^0.5)
      L[a,b] <- (D[a,a]^-0.5)*L[a,b]*(D[b,b]^-0.5)
    }
  }
  return(L)
}

KNN <- function(adj_matrix,sim1,sim2){
  kk <- 1     
  Ym <- matrix(nrow = dim(adj_matrix)[1],ncol = dim(adj_matrix)[2])
  Yd <- matrix(nrow = dim(adj_matrix)[1],ncol = dim(adj_matrix)[2])
  for(a in 1:dim(adj_matrix)[1]){
    tmp_index <- order(sim1[a,],decreasing = TRUE)[1:kk]
    tmp_value <- sim1[a,][order(sim1[a,],decreasing = TRUE)[1:kk]]
    tmpSum <- vector(mode = "numeric",length = dim(adj_matrix)[2])
    tmpSim <- 0
    for(aa in 1:kk){
      tmpSum <- tmpSum + (0.5^(aa-1))*tmp_value[aa]*adj_matrix[tmp_index[aa],]
    }
    for(a2 in 1:kk){
      tmpSim <- tmpSim + sim1[a,tmp_index[a2]]
    }
    Ym[a,] <- tmpSum/tmpSim
  }
  for(b in 1:dim(adj_matrix)[2]){
    tmp_index <- order(sim2[b,],decreasing = TRUE)[1:kk]
    tmp_value <- sim2[b,][order(sim2[b,],decreasing = TRUE)[1:kk]]
    tmpSum <- vector(mode = "numeric",length = dim(adj_matrix)[1])
    tmpSim <- 0
    for(bb in 1:kk){
      tmpSum <- tmpSum + (0.5^(bb-1))*tmp_value[bb]*adj_matrix[,tmp_index[bb]]
    }
    for(b2 in 1:kk){
      tmpSim <- tmpSim + sim2[b,tmp_index[b2]]
    }
    Yd[,b] <- tmpSum/tmpSim
  }
  Ymd <- (Ym+Yd)/2
  Y_new <- pmax(adj_matrix,Ymd)
  return(Y_new)
}

model <- function(Y_new){
  
  Ysd_svd <- svd(Ysd)
  Ymd_svd <- svd(Ymd)
  Ysm_svd <- svd(Y_new)

  diag_A <- matrix(0,nrow = dim(Ysd)[1], ncol = dim(Ysd)[1])
  diag_B <- matrix(0,nrow = dim(Ymd)[2], ncol = dim(Ymd)[2])
  diag_C <- matrix(0,nrow = dim(Ysm)[1], ncol = dim(Ysm)[1])

  for(i in 1:dim(diag_A)[1]){
    diag_A[i,i] <- Ysd_svd$d[i]
  }
  for(j in 1:dim(diag_B)[1]){
    diag_B[j,j] <- Ymd_svd$d[j]
  }
  for(k in 1:dim(diag_C)[1]){
    diag_C[k,k] <- Ysm_svd$d[k]
  }
  
  A <- to_one(Ysd_svd$u)%*%diag_A%*%to_one(t(Ysd_svd$v))
  B <- to_one(Ymd_svd$u)%*%diag_B%*%to_one(t(Ymd_svd$v))
  C <- to_one(Ysm_svd$u)%*%diag_C%*%to_one(t(Ysm_svd$v))

  OA <- matrix(data = 1,nrow=dim(A)[1],ncol=dim(A)[2])
  OB <- matrix(data = 1,nrow=dim(B)[1],ncol=dim(B)[2])
  OC <- matrix(data = 1,nrow=dim(C)[1],ncol=dim(C)[2])
  
  Ysd_new <- KNN(Ysd,Ss,Sd)
  Ymd_new <- KNN(Ymd,Sm,Sd)
  Y_new_new <- KNN(Y_new,Ss,Sm)

  Ss_new <- to_laplacian_one(Ss)
  Sm_new <- to_laplacian_one(Sm)
  Sd_new <- to_laplacian_one(Sd)
  
  
  iterator <- 0
  while(iterator<50){
    C_numerator <- 4*Ss_new%*%C + 4*gamma_C*C%*%Sm_new + 2*alpha_C*Y_new_new*Y_new_new + 2*beta*A%*%t(B)
    C_denominator <- 4*C%*%t(C)%*%C + 4*gamma_C*C%*%t(C)%*%C + 2*alpha_C*Y_new_new*C + 2*beta*C + delta_C*OC
    C_tem <- C_numerator/C_denominator
    C  <- C_tem*C

    B_numerator <- 4*gamma_B*Sm_new%*%B +4*gamma_B*B%*%Sd_new + 2*alpha_B*Ymd_new*Ymd_new + 2*beta*t(C)%*%A
    B_denominator <- 8*gamma_B*B%*%t(B)%*%B + 2*alpha_B*Ymd_new*B + 2*beta*B%*%t(A)%*%A + delta_B*OB
    B_tem <- B_numerator/B_denominator
    B <- B_tem*B

    A_numerator <- 4*gamma_A*Ss_new%*%A + 4*gamma_A*A%*%Sd_new + 2*alpha_A*Ysd_new*Ysd_new + 2*beta*C%*%B
    A_denominator <- 8*gamma_A*A%*%t(A)%*%A + 2*alpha_A*Ysd_new*A + 2*beta*A%*%t(B)%*%B + delta_A*OA
    A_tem <- A_numerator/A_denominator
    A <- A_tem*A
    
    iterator <- iterator+1
  }
  return(C)
}

setwd("E:/")

f_drug_drug1 <- read.table("SS-Similarity_Matrix_Drug_ATC.TXT")
f_drug_drug2 <- read.table("SS-Similarity_Matrix_Drug_Smiles.txt")
f_drug_drug <- (f_drug_drug1+f_drug_drug2)/2   

f_miRNA_miRNA1 <- read.table("SS-Similarity_Matrix_miRNA_BMA.txt")  #BMA算出来的相似性
f_miRNA_miRNA2 <- read.table("SS-Similarity_Matrix_miRNA_seq.txt")
f_miRNA_miRNA <- (f_miRNA_miRNA1+f_miRNA_miRNA2)/2

f_disease_disease <- read.table("SS-Similarity_Matrix_disease.txt")

f_drug_miRNA <- read.table("SS-Matrix_Drug_miRNA.txt")
f_drug_disease <- read.table("SS-Matrix_drug_disease.txt")
f_miRNA_disease <- read.table("SS-Matrix_miRNA_disease.txt")

Ss <- as.matrix(f_drug_drug)
Sm <- as.matrix(f_miRNA_miRNA)
Sd <- as.matrix(f_disease_disease)

Ysm <- as.matrix(f_drug_miRNA)
Ysd <- as.matrix(f_drug_disease)
Ymd <- as.matrix(f_miRNA_disease)


alpha_A <- alpha_B <- alpha_C <- 0.01
gamma_A <- gamma_B <- gamma_C <- 5
delta_A <- delta_B <- delta_C <- 100
beta <- 5

position <- matrix(nrow = sum(Ysm),ncol = 3)
p <- 1

for(i in 1:dim(Ysm)[1]){
  for(j in 1:dim(Ysm)[2]){
    if(Ysm[i,j]==1) {
      position[p,1] <- i
      position[p,2] <- j
      position[p,3] <- Ysm[i,j]
      p <- p+1
    }
  }
}

tep_pos_set <- sample(dim(position)[1],dim(position)[1])  
num_tep <- floor(dim(position)[1]*0.2)
q <- 5
TPR_n_all <- matrix(0,dim(Ysm)[1]*dim(Ysm)[2]-(sum(Ysm)-num_tep)+1,q)
FPR_n_all <- matrix(0,dim(Ysm)[1]*dim(Ysm)[2]-(sum(Ysm)-num_tep)+1,q)
PRE_n_all <- matrix(0,dim(Ysm)[1]*dim(Ysm)[2]-(sum(Ysm)-num_tep)+1,q)
GM_n_all <- matrix(0,dim(Ysm)[1]*dim(Ysm)[2]-(sum(Ysm)-num_tep)+1,q)
ACC_n_all <- matrix(0,dim(Ysm)[1]*dim(Ysm)[2]-(sum(Ysm)-num_tep)+1,q)
SPE_n_all <- matrix(0,dim(Ysm)[1]*dim(Ysm)[2]-(sum(Ysm)-num_tep)+1,q)
for (x in 1:q) {
  t_p <- 1
  Y_new <- Ysm
  for (j in ((x-1)*num_tep+1):(x*num_tep)) {
    Y_new[position[tep_pos_set[j],1],position[tep_pos_set[j],2]] <- 0
  }
  Y_newAB<- model(Y_new)
  #next, we calculate the AUC
}
