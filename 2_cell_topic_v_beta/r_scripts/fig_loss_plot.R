library(ggplot2)
library(gridExtra)
library(reshape)
library(yaml)
library(RColorBrewer)
source("Util.R")
options(warn=-1)


plot_cell_topic_loss <- function(loss_file,f) {
    print(loss_file)
    df = read.table(loss_file, sep = ",", header=TRUE)
    print(head(df))
    colnames(df) = c("Log-likelihood","KL loss","KLB loss")
    df$epoch <- 1:nrow(df)
    dfm = melt(df,id="epoch")
    
    p1 <-
      .gg.plot(dfm[dfm$variable=="Log-likelihood",], aes(x=epoch, y=value)) +  
        geom_point(stroke = 0, color="gray", size=1) +
        geom_smooth(color="red", se=FALSE, size=1) +
        labs(x = "epoch", title = "", y = "Log-likelihood")
      

    p2 <-
      .gg.plot(dfm[dfm$variable=="KL loss",], aes(x=epoch, y=value)) +   geom_point(stroke = 0, color="gray", size=1) +
      geom_smooth(color="red", se=FALSE, size=1) +
      labs(x = "epoch", title = "", y = "KL loss")

    p3 <-
      .gg.plot(dfm[dfm$variable=="KLB loss",], aes(x=epoch, y=value)) +   geom_point(stroke = 0, color="gray", size=1) +
      geom_smooth(color="red", se=FALSE, size=1) +
      labs(x = "epoch", title = "", y = "KLB loss")


    plotlist = list()
    plotlist[[1]] = p1 
    plotlist[[2]] = p2 
    plotlist[[3]] = p3

    stplt <- grid.arrange(grobs=plotlist,ncol=3,
    heights = c(1/3, 1/3, 1/3))
    ggsave(f,stplt)
}


plot_interaction_topic_loss <- function(loss_file,f){

    df = read.table(loss_file, sep = ",", header=FALSE)

    df = df[,3:dim(df)[2]]
    colnames(df) = c("Log-likelihood","KL loss ligands","KL loss receptors","KL loss beta1","KL loss beta2")
    df$epoch <- 1:nrow(df)
    dfm = melt(df,id="epoch")
    
    p1 <-
      .gg.plot(dfm[dfm$variable=="Log-likelihood",], aes(x=epoch, y=value)) +  
        geom_point(stroke = 0, color="gray", size=1) +
        geom_smooth(color="red", se=FALSE, size=1) +
        labs(x = "Epoch", title = "", y = "Log-likelihood")
      

    p2 <-
      .gg.plot(dfm[dfm$variable=="KL loss ligands",], aes(x=epoch, y=value)) +   geom_point(stroke = 0, color="gray", size=1) +
      geom_smooth(color="red", se=FALSE, size=1) +
      labs(x = "Epoch", title = "", y = "KL loss ligands")

    p3 <-
      .gg.plot(dfm[dfm$variable=="KL loss receptors",], aes(x=epoch, y=value)) +   geom_point(stroke = 0, color="gray", size=1) +
      geom_smooth(color="red", se=FALSE, size=1) +
      labs(x = "Epoch", title = "", y = "KL loss receptors")

    p4 <-
      .gg.plot(dfm[dfm$variable=="KL loss beta1",], aes(x=epoch, y=value)) +   geom_point(stroke = 0, color="gray", size=1) +
      geom_smooth(color="red", se=FALSE, size=1) +
      labs(x = "Epoch", title = "", y = "KL loss beta1")

    p5 <-
      .gg.plot(dfm[dfm$variable=="KL loss beta2",], aes(x=epoch, y=value)) +   geom_point(stroke = 0, color="gray", size=1) +
      geom_smooth(color="red", se=FALSE, size=1) +
      labs(x = "Epoch", title = "", y = "KL loss beta2")


    plotlist = list()
    plotlist[[1]] = p1 
    plotlist[[2]] = p2 
    plotlist[[3]] = p3 
    plotlist[[4]] = p4 
    plotlist[[5]] = p5 

    stplt <- grid.arrange(grobs=plotlist,ncol=5,heights = c(1/5, 1/5, 1/5,1/5,1/5))

    ggsave(f,stplt)
}