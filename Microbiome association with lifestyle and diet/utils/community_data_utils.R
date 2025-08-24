# custom functions
convert_to_CPM <- function(comm.data){
  comm.df <- comm.data # abundance matrix with columns as sampleid
  comm.df <- comm.df*10^6 
}

prevalence_filter <- function(x, n) {
  x <- reads # matrix of abundance with sampleid as columns
  y <- x
  y[y>0] <- 1
  y <- y[rowSums(y)>n,]
  x <- x[rownames(x) %in% rownames(y),]
  x
}

headx <- function(data.frame,n=NULL){
  if(is.null(n)){
    data.frame[1:5,1:5]
  }else{
    data.frame[1:n,1:n]
  }
}
