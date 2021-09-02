# make working directory
mkdir -p PATH/TO/WORKIING/DIR
cd PATH/TO/WORKIING/DIR
# Under this directory, make a 'fastq' directory containing fastq files.

# trimgalore
cd fastq
for dir in `ls | grep ^Y`;do
  cd ${dir}
  R1=`ls *1.fq`
  R2=`ls *2.fq`
  trim_galore --paired ${R1} ${R2}
  gzip ${R1}
  gzip ${R2}
  cd ..
done

# make assembly
cd PATH/TO/WORKIING/DIR/fastq
list=`find . -name "*val_1.fq"`
left=`echo ${list}|sed 's/ /,/g'`
list=`find . -name "*val_2.fq"`
right=`echo ${list}|sed 's/ /,/g'`
Trinity --seqType fq --max_memory 60G --left ${left} --right ${right} --CPU 6

# quantification
mkdir -p PATH/TO/WORKIING/DIR/ref
cd PATH/TO/WORKIING/DIR/ref
ln -s PATH/TO/Trinity.fasta .
kallisto index -i kallisto.idx Trinity.fasta
cd ..
mkdir -p quant
cd quant
find ../fastq -name "*val_1.fq"|while read fn;do
  left=$fn
  right=`echo $fn|sed 's/L2_1_val_1/L2_2_val_2/g'`
  sample=`echo $fn|cut -f3 -d"/"` # Extract sample name (e.g., Y1, Y2, etc). Modify upon your directory structure.
  cmd="kallisto quant -i ../ref/kallisto.idx -o ${sample} --bootstrap-samples=100 --threads=4 --pseudobam ${left} ${right}"
  echo ${cmd}
done

# DEG using edgeR
mkdir -p PATH/TO/WORKIING/DIR/deg
cd PATH/TO/WORKIING/DIR/deg
## run deg.R

# extract deg trinity ids (FDR-p<=0.05)
awk 'BEGIN{OFS="\t"}{if($6 <= 0.05 && $2 >= 3){print $1}}' ../LRT.edgeR.res > transcriptID_fdr0.05edgeR.upregulated.list
awk 'BEGIN{OFS="\t"}{if($6 <= 0.05 && $2 <= -3){print $1}}' ../LRT.edgeR.res > transcriptID_fdr0.05edgeR.downregulated.list

# extract sequence of deg contigs
seqtk subseq PATH/TO/Trinity.fasta transcriptID_fdr0.05edgeR.upregulated.list > transcript_fdr0.05edgeR.upregulated.fasta.tmp
seqtk subseq PATH/TO/Trinity.fasta transcriptID_fdr0.05edgeR.downregulated.list > transcript_fdr0.05edgeR.downregulated.fasta.tmp
cut -f1 -d" " transcript_fdr0.05edgeR.upregulated.fasta.tmp > transcript_fdr0.05edgeR.upregulated.fasta
cut -f1 -d" " transcript_fdr0.05edgeR.downregulated.fasta.tmp > transcript_fdr0.05edgeR.downregulated.fasta
rm transcript_fdr0.05edgeR.*.fasta.tmp

# predict ORF of degs
TransDecoder.LongOrfs -t transcript_fdr0.05edgeR.upregulated.fasta
TransDecoder.LongOrfs -t transcript_fdr0.05edgeR.downregulated.fasta
# final coding region predictions
TransDecoder.Predict -t transcript_fdr0.05edgeR.upregulated.fasta --single_best_only --cpu 0
TransDecoder.Predict -t transcript_fdr0.05edgeR.downregulated.fasta --single_best_only --cpu 0

# blastp
# make geneID-geneName list
# download X.tropicalis cdna data
makeblastdb -in PATH/TO/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.pep.all.fa -out Xenopus_tropicalis.Xenopus_tropicalis_v9.1.pep.all -dbtype prot -parse_seqids
blastp -db ~/analysis/pubdata/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.pep.all -query transcript_fdr0.05edgeR.upregulated.fasta.transdecoder.pep -max_target_seqs 1 -outfmt 6 -num_threads 8 -out transcript_fdr0.05edgeR.upregulated.transdecoder.blastp.res
blastp -db ~/analysis/pubdata/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.pep.all -query transcript_fdr0.05edgeR.downregulated.fasta.transdecoder.pep -max_target_seqs 1 -outfmt 6 -num_threads 8 -out transcript_fdr0.05edgeR.downregulated.transdecoder.blastp.res
# reshape
for type in up down;do
 # blastx
 cut -f1 transcript_fdr0.05edgeR.${type}regulated.transdecoder.blastp.res | uniq | while read contig;do
  grep ${contig} transcript_fdr0.05edgeR.${type}regulated.transdecoder.blastp.res | awk 'BEGIN{OFS="\t"}{if($11 < 1e-5){print}}' | head -n1
 done > transcript_fdr0.05edgeR.${type}regulated.transdecoder.blastp.Eval1e-5.tophit.res
done

# prepare gene list
grep ">" PATH/TO/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.pep.all.fa | while read line;do
 res=`echo "${line}" | tr " " "\n" | grep "ENSXETP\|gene:\|transcript:\|gene_symbol:" | sed -e 's/\>//g' -e 's/gene://g' -e 's/gene_symbol://g' -e 's/transcript://g' | tr "\n" "\t"`
 echo "${res}"
done | awk 'BEGIN{OFS="\t"}{if(NF==4){print $1,$2,$3,$4}else{print}}' > proteinID_geneID_transcriptID_geneName.tsv

# PPI downstream analysis
## download PPI datasets (simple tabular text output and protein annotations) from STRING
# reshape protein annotation file
cat string_protein_annotations.tsv | cut -f1-2 | sed -e 's/8364\.//g' -e 's/\#//g'> string_protein_annotations2.tsv
# run PPI.R

