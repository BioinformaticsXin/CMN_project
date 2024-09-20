#Use CNVkit for WES data
trim_galore -j 4 --gzip -o /Mutation/2_trim_galore/D2_CA/ --paired '/WES/01.Raw_Data/D2_CA_HC3CMDSXY_1.fq.gz' '/WES/01.Raw_Data/D2_CA_HC3CMDSXY_2.fq.gz'
trim_galore -j 4 --gzip -o /Mutation/2_trim_galore/D2_CA/ --paired '/WES/01.Raw_Data/D2_CA_HCJWFDSXY_1.fq.gz' '/WES/01.Raw_Data/D2_CA_HCJWFDSXY_2.fq.gz'
trim_galore -j 4 --gzip -o /Mutation/2_trim_galore/D2_CJ/ --paired '/WES/01.Raw_Data/D2_CJ_HC3CMDSXY_1.fq.gz' '/WES/01.Raw_Data/D2_CJ_HC3CMDSXY_2.fq.gz'
trim_galore -j 4 --gzip -o /Mutation/2_trim_galore/D2_CJ/ --paired '/WES/01.Raw_Data/D2_CJ_HCJWFDSXY_1.fq.gz' '/WES/01.Raw_Data/D2_CJ_HCJWFDSXY_2.fq.gz'
trim_galore -j 4 --gzip -o /Mutation/2_trim_galore/D2_PB/ --paired '/WES/01.Raw_Data/D2_PB_HCJWFDSXY_1.fq.gz' '/WES/01.Raw_Data/D2_PB_HCJWFDSXY_2.fq.gz'
trim_galore -j 4 --gzip -o /Mutation/2_trim_galore/D3_CA/ --paired '/WES/01.Raw_Data/D3_CA_HC3CMDSXY_1.fq.gz' '/WES/01.Raw_Data/D3_CA_HC3CMDSXY_2.fq.gz'
trim_galore -j 4 --gzip -o /Mutation/2_trim_galore/D3_CA/ --paired '/WES/01.Raw_Data/D3_CA_HCJWFDSXY_1.fq.gz' '/WES/01.Raw_Data/D3_CA_HCJWFDSXY_2.fq.gz'
trim_galore -j 4 --gzip -o /Mutation/2_trim_galore/D3_CJ/ --paired '/WES/01.Raw_Data/D3_CJ_HCJWFDSXY_1.fq.gz' '/WES/01.Raw_Data/D3_CJ_HCJWFDSXY_2.fq.gz'

#BWA
bwa index -a bwtsw -p gatk_hg38 /data/Homo_sapiens_assembly38.fasta
for i in D2_CA D2_CJ D2_PB D3_CA D3_CJ;do
  ll /Mutation/2_trim_galore/${i}/ |grep -v "txt">>/clean/cong.txt
done

for i in D2_CA D2_CJ D2_PB D3_CA D3_CJ;do
  cp /Mutation/2_trim_galore/${i}/*gz /clean/
done

INDEX=/ref/GENOME/hg38.fa
cat /clean/cong.txt |while read id
do
  arr=($id)
  fq1=${arr[1]}
  fq2=${arr[2]}
  sample=${arr[0]}
  echo ${sample} ${fq1} ${fq2} 
  bwa mem -t 5 -R "@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:Illumina" ${INDEX} ${fq1} ${fq2}  | samtools sort -@ 5 -o /align/${sample}.bam -   
done 

gunzip -c gencode.v40.annotation.gtf.gz | grep 'transcript_type "protein_coding"' |  awk '($3=="exon") {printf("%s\t%s\t%s\n",$1,int($4)-1,$5);}' | sort -T . -t $'\t' -k1,1 -k2,2n |  bedtools merge > hg38_EXON.bed

cd /align
for sample in D2_PB D2_CJ
do 
cnvkit.py batch *D2_CA.bam    \
    --normal *${sample}.bam     \
    --output-reference /CNVkit/${sample}/my_reference.cnn \
    --output-dir /CNVkit/${sample} \
    --targets /ref/hg38_EXON.bed           \
    --fasta /ref/GENOME/hg38.fa           \
    -p 20               \
    --method amplicon             \
    --drop-low-coverage         \
    --scatter --diagram         \
    1>/CNVkit/${sample}/cnvkit.log 2>&1

cnvkit.py heatmap *D2_CA.cns >${sample}.pdf	
done

cnvkit.py batch *D3_CA.bam    \
    --normal *D3_CJ.bam     \
    --output-reference /CNVkit/D3_CJ/my_reference.cnn \
    --output-dir /CNVkit/D3_CJ \
    --targets /ref/hg38_EXON.bed           \
    --fasta /ref/GENOME/hg38.fa           \
    -p 20               \
    --method amplicon             \
    --drop-low-coverage         \
    --scatter --diagram         \
    1>/CNVkit/D3_CJ/cnvkit.log 2>&1
cnvkit.py heatmap *D3_CA.cns >D3_CA.pdf	


cnvkit.py batch *D2_CA.bam    \
    --normal *D2_PB.bam *D2_CJ.bam    \
    --output-reference /CNVkit/D2_PB_CJ/my_reference.cnn \
    --output-dir /CNVkit/D2_PB_CJ \
    --targets /ref/hg38_EXON.bed           \
    --fasta /ref/GENOME/hg38.fa           \
    -p 20               \
    --method amplicon             \
    --drop-low-coverage         \
    --scatter --diagram         \
    1>/CNVkit/D2_PB_CJ/cnvkit.log 2>&1
cnvkit.py heatmap *D2_CA.cns >D2_CA.pdf	
