library(sleuth)
library(tximport)
library(edgeR)
library(ggplot2)
setwd("PATH/TO/WORKIING/DIR/deg")

s2c <- data.frame(path=list.dirs("../quant")[2:7])
s2c$sample <- sub("../quant/","",s2c$path)
s2c$condition <- "case"
s2c$condition[grep("cont",s2c$sample)] <- "ctrl"
s2c$path <- as.character(s2c$path)
so <- sleuth_prep(s2c,extra_bootstrap_summary=T)
kallisto.df <- kallisto_table(so)
write.table(kallisto.df,"kallisto_res.txt",row.names=F,sep="\t",quote=F)

# edgeR
kallisto.files <- paste0(s2c$path,"/abundance.tsv")
names(kallisto.files) <- s2c$sample
tx.exp <- tximport(kallisto.files, type = "kallisto", txOut = TRUE)
counts_raw <- as.data.frame(tx.exp$counts)
counts_raw$transcript <- rownames(counts_raw)
counts <- round(tx.exp$counts)
nTranscript <- nrow(counts)
group <- factor(s2c$condition)
design <- model.matrix(~ group)
d <- DGEList(counts = counts, group = group)
d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
fit <- glmFit(d, design)
nGroup <- length(unique(group))
# two-group comparison
lrt <- glmLRT(fit,coef=2)
res <- as.data.frame(topTags(lrt,n=nTranscript))
res$transcript <- rownames(res)
res <- merge(res,counts_raw,by="transcript")
edgeR.res <- res <- res[order(res$FDR),]
write.table(res,"LRT.edgeR.res",row.names=F,quote=F,sep="\t")

# volcano plot
df <- read.table("LRT.edgeR.res",header=T,stringsAsFactors=F)
df2 <- df[,c("transcript","logFC","FDR")]
df2$sig <- "F"
df2$sig[df2$FDR <= 0.05 & abs(df2$logFC)>=3] <- "T"
df3=df2
df3$sig[df3$FDR <= 0.05 & df3$logFC>=3] <- "U"
df3$sig[df3$FDR <= 0.05 & df3$logFC<=-3] <- "D"
df3$shape <- "normal"
df3$shape[df3$FDR < 1e-5] <- "triangle"
df3$FDR2 <- df3$FDR
df3$FDR2[df3$FDR < 1e-05] <- 1e-05
g <- ggplot()+geom_hline(yintercept=-log10(0.05),linetype="dotted",colour="gray50")+geom_vline(xintercept=c(-3,3),linetype="dotted",colour="gray50")+geom_point(data=df3[df3$sig=="F",],aes(x=logFC,y=-log10(FDR2),shape=shape),colour="#c1c1c1",size=1)+geom_point(data=df3[df3$sig=="U" | df3$sig=="D" ,],aes(x=logFC,y=-log10(FDR2),shape=shape,colour=sig),size=1)+scale_colour_manual(values=c("U"="#d11141","D"="#00aedb"))+theme_bw()+theme(legend.position="none")+xlab(expression(paste(log[2](FC))))+ylab(expression(paste(-log[10](italic(q)))))+theme_classic()+theme(strip.background=element_rect(fill=NA, color=NA),legend.position="NONE")
pdf("volcanoe.pdf", useDingbats=FALSE,width=3,height=3)
g
dev.off()
