library(tidyverse)
library(reshape2)
library(igraph)
library(gprofiler2)

# reshape and export dataset
prot=read.table("proteinID_geneID_transcriptID_geneName.tsv",header=F,stringsAsFactors=F,sep="\t")
blastp_up=read.table("transcript_fdr0.05edgeR.upregulated.transdecoder.blastp.Eval1e-5.tophit.res",header=F,sep="\t",stringsAsFactors=F)
blastp_down=read.table("transcript_fdr0.05edgeR.downregulated.transdecoder.blastp.Eval1e-5.tophit.res",header=F,sep="\t",stringsAsFactors=F)
blast=rbind(blastp_up[,1:2],blastp_down[,1:2])
colnames(blast)=c("contigID","proteinID_blastp")
blast2=merge(blast,prot,by.x="proteinID_blastp",by.y="V1",all.x=T,all.y=F)
colnames(blast2)[3:5]=c("geneID","transcriptID","geneSymbol")
blast2$contigID2=sub("\\.p[0-9]*","",blast2$contigID)
df2 <- df[,c("transcript","logFC","FDR")]
df2$sig <- "F"
df2$sig[df2$FDR <= 0.05 & abs(df2$logFC)>=3] <- "T"
df3=df2[df2$sig=="T",c("transcript","logFC","FDR")]
colnames(df3)[1]="contigID"
df4=merge(df3,blast2,by.x="contigID",by.y="contigID2",all.x=T,all.y=F)
df4$type="upregulated"
df4[df4$logFC<0,"type"]="downregulated"
dfup=df4[df4$type=="upregulated",]
dfdown=df4[df4$type=="downregulated",]
write.table(rbind(dfup[order(dfup$FDR),],dfdown[order(dfdown$FDR),]),"resultTable.tsv",row.names=F,quote=F,sep="\t")

# load STRING results
df1=read.table("string_interactions.tsv",header=T,sep="\t",stringsAsFactors=F)
deg=read.table("resultTable.tsv",sep="\t",header=T,stringsAsFactors=F)
annot=read.table("string_protein_annotations2.tsv",header=T,stringsAsFactors=F)
df1$node1=toupper(df1$node1)
df1$node2=toupper(df1$node2)
prot=unique(c(df1$node1,df1$node2))
df2=data.frame(prot=prot,stringsAsFactors=F)
for(p1 in prot){
 count=NULL
 for(p2 in prot){
  n=nrow(df1[df1$node1==p1 & df1$node2==p2,])
  count=c(count,n)
 }
 df2=cbind(df2,count)
 colnames(df2)[ncol(df2)]=p1
}
rownames(df2)=df2$prot
df3=df2[,-1]
dfout=df2
df4=graph.adjacency(as.matrix(df3), mode = "undirected")
df4=simplify(df4)
# add up/down info
genelist=toupper(V(df4)$name)
deg2=na.omit(deg[,c("proteinID_blastp","geneSymbol","type")])
deg2$id2=sub("\\.[0-9]","",deg2$proteinID_blastp)
# merge by geneSymbol
tmp1=unique(deg2[toupper(deg2$geneSymbol)%in%genelist,c("geneSymbol","type")])
tmp1$geneSymbol2=toupper(tmp1$geneSymbol)
res1=tmp1[tmp1$geneSymbol2%in%genelist,c("geneSymbol2","type")]
# merge by ens id
annot2=unique(annot[toupper(annot$node)%in%genelist,])
tmp2=unique(merge(annot2,deg2,by.x="identifier",by.y="id2",all=T))
tmp2$node2=toupper(tmp2$node)
res2=tmp2[tmp2$node2%in%genelist,c("node2","type")]
colnames(res2)[1]="geneSymbol2"
# summary
res=na.omit(unique(rbind(res1,res2)))
rownames(res)=res$geneSymbol2
# Three proteins have conflicting symbol between databases...
setdiff(genelist,res$geneSymbol2)
# [1] "PHKA1"   "RAPGEF2" "UNC13D" 
# Here, I add info whether these genes were upregulated or downregulated
res=rbind(res,c("PHKA1","upregulated"),c("RAPGEF2","upregulated"),c("RAPGEF2","downregulated"),c("UNC13D","downregulated"))
# There are genes which appeared in both upregureted and downregulated DEG list
down0=res[grep("downregulated",res$type),"geneSymbol2"]
up0=res[grep("upregulated",res$type),"geneSymbol2"]
down=setdiff(down0,up0)
up=setdiff(up0,down0)
both=intersect(up0,down0)
res2=data.frame(gene=c(up,down,both),type=c(rep("up",length(up)),rep("down",length(down)),rep("both",length(both))),stringsAsFactors=F)
rownames(res2)=res2$gene

# degree centrality
degree=as.data.frame(degree(df4),stringsAsFactors=F)
colnames(degree)="degree"
dfout=cbind(dfout,degree)
V(df4)$degree=degree(df4) # assignment
V(df4)$type=res2[genelist,"type"]
V(df4)$frame.color <- "white"
pdf("ppi_degree1.pdf",width=10,height=10)
plot(df4, vertex.label.cex = .6, vertex.label.color = "black", vertex.size = V(df4)$degree, vertex.label.cex = .1) # sized by degree
dev.off()
pdf("ppi_degree2.pdf",width=6,height=6)
plot(df4, 
     vertex.label.cex = .6, 
     vertex.label.color = "black", 
     vertex.size = V(df4)$degree*4, # use a scalar to increase the value of the degree but maintain the ratio.
     vertex.color = c( "gray", "skyblue", "pink")[as.numeric(as.factor(V(df4)$type))]
)
dev.off()

# betweenness centrality
betweenness=as.data.frame(betweenness(df4, directed = FALSE, normalized=T),stringsAsFactors=F)
colnames(betweenness)="betweenness"
dfout=cbind(dfout,betweenness)
V(df4)$betweenness <- betweenness(df4, directed = F) # assignment
pdf("ppi_betweenness1.pdf",width=10,height=10)
plot(df4, 
     vertex.label.cex = .6, 
     vertex.label.color = "black", 
     vertex.size = V(df4)$betweenness) # sized by betweenness
dev.off()
pdf("ppi_betweenness2.pdf",width=10,height=10)
plot(df4,
     vertex.label.cex = .6, 
     vertex.label.color = "black", 
     vertex.size = V(df4)$betweenness/max(V(df4)$betweenness) * 20) # Normalize sizes by dividing by the maximum and multiplying by some scalar when plotting.
dev.off()

# eigenvector
eigenvector=as.data.frame(evcent(df4)$vector,stringsAsFactors=F)
colnames(eigenvector)="eigenvector"
dfout=cbind(dfout,eigenvector)
V(df4)$eigenvector <- evcent(df4)$vector
plot(df4,
     vertex.label.cex = .6, 
     vertex.label.color = "black", 
     vertex.size = V(df4)$eigenvector/max(V(df4)$eigenvector) * 20)

# export
dfdegree=dfout[order(dfout$degree,decreasing=T),c("prot","degree")]
colnames(dfdegree)[1]="protein"
write.table(dfdegree,"ppi_degreeCentrality.tsv",row.names=F,quote=F,sep="\t")

dfeigenvector=dfout[order(dfout$eigenvector,decreasing=T),c("prot","eigenvector")]
colnames(dfeigenvector)[1]="protein"
write.table(dfeigenvector,"ppi_eigenvectorCentrality.tsv",row.names=F,quote=F,sep="\t")

dfbetweenness=dfout[order(dfout$betweenness,decreasing=T),c("prot","betweenness")]
colnames(dfbetweenness)[1]="protein"
write.table(dfbetweenness,"ppi_betweennessCentrality.tsv",row.names=F,quote=F,sep="\t")

dfout=dfout[,c("prot","degree","eigenvector","betweenness")]
colnames(dfout)[1]="protein"
write.table(dfout,"ppi_igraphRes.tsv",row.names=F,quote=F,sep="\t")
