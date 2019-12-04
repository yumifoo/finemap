merge_region <- function(pos_matrix, window_mb =NULL,data = NULL){
  # pos_matrix is a two column matrix for the region position
  ## 1st column is the left bound of the region
  ## 2nd column is the second bound of the region
  ## position should be in base pairs(bp)
  # plot_res, logical option, wheather to plot the manhattan plot with merged area being marked
  # window_mb, window size used for picking peaks, required if plot_res = TRUE
  # data, summary statistics of GWAS result, this is required if you want plot the result
  ## it should at least contain SNP, CHR, BP, PVAL
  mergeBP <- pos_matrix
  if(is.null(nrow(mergeBP))){ # that is no significant P vals
    return(mergeBP)
  }
  if(!is.null(nrow(mergeBP))){
    if(nrow(mergeBP) == 1|nrow(mergeBP) == 0){
      return(mergeBP)
    }
    if(nrow(mergeBP) > 1){
      mergeBP <- mergeBP[order(mergeBP[,1]),]
      for (i in 1:nrow(mergeBP)){
        # pick row i and compare it with the other rows below i (i.e. compare with row j)
        # if a merge is needed, update i row and marked merged j row as (-1,-1)
        for (j in 1:nrow(mergeBP)){
          if(mergeBP[i,1] < mergeBP[i,2] & mergeBP[j,1] < mergeBP[j,2]){
            if(mergeBP[i,1] < mergeBP[j,1] & mergeBP[j,1] < mergeBP[i,2]){
              mergeBP[i,2] <- mergeBP[j,2]
              mergeBP[j,] <- c(-1,-1)
            } 
            if(mergeBP[i,1] < mergeBP[j,2] & mergeBP[j,2] < mergeBP[i,2]){
              mergeBP[i,1] <- mergeBP[j,1]
              mergeBP[j,] <- c(-1,-1)
            }
          }
          if(mergeBP[i,1] == mergeBP[i,2] & mergeBP[j,1] < mergeBP[j,2]){
            if(mergeBP[j,1]<mergeBP[i,1] & mergeBP[i,1]<mergeBP[j,2]){
              mergeBP[i,1] <- mergeBP[j,1]
              mergeBP[i,2] <- mergeBP[j,2]
              mergeBP[j,] <- c(-1,-1)
            }
          }
          if(mergeBP[i,1] < mergeBP[i,2] & mergeBP[j,1]== mergeBP[j,2]){
            if(mergeBP[i,1] < mergeBP[j,1] & mergeBP[j,1] < mergeBP[i,2]){
              mergeBP[j,] <- c(-1,-1)
            }
          }
        }
      }
      # remove the rows that had been merged to other rows, i.e. remove the rows with (-1,-1) values
      mergeBP <- mergeBP[mergeBP[,1] != -1,]
      if(is.null(nrow(mergeBP))){ # that is only one region after merge
        return(mergeBP)
      }
      
      if(!is.null(nrow(mergeBP))){
        mergeBP <- mergeBP[order(mergeBP[,1]),]
        return(mergeBP)
      }
    } 
  }
}