

treatments = factor(1:6)
subBlocks = factor(rep(1:4,each=6))
Blocks = factor(rep(1:3,each=8))
blocks = data.frame(Blocks,subBlocks)
design(treatments,blocks)




treatments = factor(1:6)
MainCols = factor(rep(rep(1:4,each=3),4))
MainRows = factor(rep(1:4,each=12))
Columns = factor(rep(1:12,4))
blocks = data.frame(MainCols,MainRows,Columns)
design(treatments,blocks,searches=100,weighting=0)




treatments=factor(1:6)
reps=factor(rep(1:4,each=6))
sub1=factor(rep(1:3,each=8))
blocks=data.frame(reps,sub1)
design(treatments,blocks,searches=1)



blocks(6,4,3)
blocks(6,4,c(3,4))


treatments=factor(1:6)
replicates=factor(rep(1:4,each=6))
rows=factor(rep(rep(1:3,each=4),2))
cols=factor(rep(rep(1:8,each=3),1))
blocks=data.frame(replicates,rows,cols)
design(treatments,blocks,searches=1)






repeat{
  v1 <- sample.int(6, 6, replace = F)
  V1 <- c(v1, sample.int(6, 2, replace = F))
  V2 <- c(rev(v1), sample.int(6, 2, replace = F))
  # V3 <- c(sample.int(6, 2, replace = F), sample.int(6, 6, replace = F))
  
  # V1 <- c(sample.int(6, 6, replace = F), sample.int(6, 2, replace = F))
  # V2 <- c(sample.int(6, 6, replace = F), sample.int(6, 2, replace = F))
  # V3 <- c(sample.int(6, 6, replace = F), sample.int(6, 2, replace = F))
  # V4 <- c(sample.int(6, 6, replace = F), sample.int(6, 2, replace = F))
  # mat1 <- matrix(sample.int(6, 24, replace = T), ncol=6)
  # mat2 <- matrix(sample.int(6, 8, replace = F), ncol=2)
  mat <- rbind(V1, V2)
  colUq <- sum(apply(mat, 2, function(x) sum(duplicated(x)) ))
  # colUq <- sum(apply(mat, 1, function(x) sum(duplicated(x)) ) > 2 )
  valUq <- uniqueN(mat[apply(mat, 2, duplicated)])
  
  if(colUq == 0 & valUq == 4){
    break
  }
}


repeat{
  
  v1 <- c(6,3,1,2,6,5,2,4)
  v2 <- c(3,4,6,3,2,4,1,6)
    vt <- c(1,1,2,3,4,5,5,6)
    id <- sample.int(8, 8, replace = F)
  v3 <- vt[id]
  mat <- rbind(v1,v2,v3)
  colUq <- sum(apply(mat, 2, function(x) sum(duplicated(x)) ))
  
  if(colUq == 0){
    break  print(mat)
  }
 
  }

  
  
  
  
  
