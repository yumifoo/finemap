peak_picker <- function(data, window_mb, plot_step =FALSE){
  library(ggplot2)
  library(gridExtra)
  #data is the summary stats, should at least contain SNP, CHR, BP, PVAL
  #window_mb is the window size around the SNP with the smallest p value for each cycle
  #plot_step is the logical option whether to generate a plot for each region pick step
  window_bp <- (window_mb*1e6)/2 # window/2 is the length that should be taken from each side of the top SNP
  ####################### Step 1: Pick the regions to be merged
  position_matrix <- vector(mode = "numeric", length = 0)
  it <- c()
  for(i in 1:nrow(data)){
    if(min(data$PVAL, na.rm = TRUE) < 5e-8){
      sig <- data[which.min(data$PVAL),]
      bp <- sig$BP
      it <- rep(i, length(bp))
      pos <- cbind(it, bp-window_bp, bp, bp+window_bp)
      position_matrix <- rbind(position_matrix, pos)
      for (j in 1:length(bp)){
        if(plot_step == TRUE){
          # The region to be removed
          region <- subset(data, BP > (bp-window_bp)[j]& BP < (bp+window_bp)[j])
          # manhattan plot for the region that will be removed
          plot1 <- ggplot(region, aes(x = BP, y = -log10(PVAL))) +
            geom_point(alpha = 0.8, size = 1.3, colour = "#999999") +
            geom_hline(yintercept=-log10(5e-8),linetype="dashed", color = "red")+
            geom_point(data = sig, stat = "identity",
                       position = "identity", color = "#56B4E9", size = 2) +
            ggtitle(paste("region", (bp-window_bp)[j], "to", (bp+window_bp)[j], ", SNP", sig$SNP, bp[j], "BP"))
          # manhattan plot show the whole x chromosome
          plot2 <- ggplot(data,aes(x = BP, y = -log10(PVAL)))+
            geom_point(alpha = 0.8, size = 1.3, colour = "#999999")+
            geom_hline(yintercept=-log10(5e-8),linetype="dashed", color = "red")+
            geom_vline(xintercept = (bp-window_bp)[j], linetype="dashed", color = "green") +
            geom_vline(xintercept = (bp+window_bp)[j], linetype="dashed", color = "green")+
            geom_point(data = sig, stat = "identity",
                       position = "identity", color = "#56B4E9", size = 2)+
            ggtitle(paste("Before removing region",  (bp-window_bp)[j], "to", (bp+window_bp)[j]))
          grid.arrange(plot1, plot2, ncol=2)
        }
        # remove the delected region
        data <- subset(data, BP < (bp-window_bp)[j]| BP > (bp+window_bp)[j])
      }
    }
    else
      break
  }
  if(is.null(nrow(position_matrix))){
    return(position_matrix)
  }
  if(!is.null(nrow(position_matrix))){
    colnames(position_matrix) <- c("iteration","left_bound","top_SNP","right_bound")
    return(position_matrix)
  }
} 