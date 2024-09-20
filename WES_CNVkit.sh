#使用CNVkit对WES数据
#去接头
trim_galore -j 4 --gzip -o /Mutation/2_trim_galore/HRX_CA/ --paired '/WES/01.Raw_Data/HRX_CA_HC3CMDSXY_1.fq.gz' '/WES/01.Raw_Data/HRX_CA_HC3CMDSXY_2.fq.gz'
trim_galore -j 4 --gzip -o /Mutation/2_trim_galore/HRX_CA/ --paired '/WES/01.Raw_Data/HRX_CA_HCJWFDSXY_1.fq.gz' '/WES/01.Raw_Data/HRX_CA_HCJWFDSXY_2.fq.gz'
trim_galore -j 4 --gzip -o /Mutation/2_trim_galore/HRX_CJ/ --paired '/WES/01.Raw_Data/HRX_CJ_HC3CMDSXY_1.fq.gz' '/WES/01.Raw_Data/HRX_CJ_HC3CMDSXY_2.fq.gz'
trim_galore -j 4 --gzip -o /Mutation/2_trim_galore/HRX_CJ/ --paired '/WES/01.Raw_Data/HRX_CJ_HCJWFDSXY_1.fq.gz' '/WES/01.Raw_Data/HRX_CJ_HCJWFDSXY_2.fq.gz'
trim_galore -j 4 --gzip -o /Mutation/2_trim_galore/HRX_PB/ --paired '/WES/01.Raw_Data/HRX_PB_HCJWFDSXY_1.fq.gz' '/WES/01.Raw_Data/HRX_PB_HCJWFDSXY_2.fq.gz'
trim_galore -j 4 --gzip -o /Mutation/2_trim_galore/YMZ_CA/ --paired '/WES/01.Raw_Data/YMZ_CA_HC3CMDSXY_1.fq.gz' '/WES/01.Raw_Data/YMZ_CA_HC3CMDSXY_2.fq.gz'
trim_galore -j 4 --gzip -o /Mutation/2_trim_galore/YMZ_CA/ --paired '/WES/01.Raw_Data/YMZ_CA_HCJWFDSXY_1.fq.gz' '/WES/01.Raw_Data/YMZ_CA_HCJWFDSXY_2.fq.gz'
trim_galore -j 4 --gzip -o /Mutation/2_trim_galore/YMZ_CJ/ --paired '/WES/01.Raw_Data/YMZ_CJ_HCJWFDSXY_1.fq.gz' '/WES/01.Raw_Data/YMZ_CJ_HCJWFDSXY_2.fq.gz'

#BWA比对得到bam
bwa index -a bwtsw -p gatk_hg38 /data/Homo_sapiens_assembly38.fasta
for i in HRX_CA HRX_CJ HRX_PB YMZ_CA YMZ_CJ;do
  ll /Mutation/2_trim_galore/${i}/ |grep -v "txt">>/clean/cong.txt
done
#修改文件
for i in HRX_CA HRX_CJ HRX_PB YMZ_CA YMZ_CJ;do
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

#其中 --targets  bed指定bed文件，对于WES数据，即捕获区间。−−fasta{GENOME} 指定参考基因组 fasta 文件，-p nthread指定线程数。而−−annotate{refFlat} 参数是可选的，不加也不影响，其指定注释文件。
#如果是全基因组测序数据，用 batch --method wgs；如果是捕获基因组测序，包括全外显子，就用 batch --method amplicon；然后一定要提供捕获区域的bed文件，一般是外显子加上其侧翼上下游的50bp长度。
gunzip -c gencode.v40.annotation.gtf.gz | grep 'transcript_type "protein_coding"' |  awk '($3=="exon") {printf("%s\t%s\t%s\n",$1,int($4)-1,$5);}' | sort -T . -t $'\t' -k1,1 -k2,2n |  bedtools merge > hg38_EXON.bed

cd /align
for sample in HRX_PB HRX_CJ
do 
cnvkit.py batch *HRX_CA.bam    \
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

cnvkit.py heatmap *HRX_CA.cns >a.pdf	
done

cnvkit.py batch *YMZ_CA.bam    \
    --normal *YMZ_CJ.bam     \
    --output-reference /CNVkit/YMZ_CJ/my_reference.cnn \
    --output-dir /CNVkit/YMZ_CJ \
    --targets /ref/hg38_EXON.bed           \
    --fasta /ref/GENOME/hg38.fa           \
    -p 20               \
    --method amplicon             \
    --drop-low-coverage         \
    --scatter --diagram         \
    1>/CNVkit/YMZ_CJ/cnvkit.log 2>&1
cnvkit.py heatmap *YMZ_CA.cns >a.pdf	


cnvkit.py batch *HRX_CA.bam    \
    --normal *HRX_PB.bam *HRX_CJ.bam    \
    --output-reference /CNVkit/HRX_PB_CJ/my_reference.cnn \
    --output-dir /CNVkit/HRX_PB_CJ \
    --targets /ref/hg38_EXON.bed           \
    --fasta /ref/GENOME/hg38.fa           \
    -p 20               \
    --method amplicon             \
    --drop-low-coverage         \
    --scatter --diagram         \
    1>/CNVkit/HRX_PB_CJ/cnvkit.log 2>&1
cnvkit.py heatmap *HRX_CA.cns >a.pdf	


#CNVkit 结合 GISTIC2  
cnvkit.py export seg HRX_CA_HCJWFDSXY.cns HRX_CA_HC3CMDSXY.cns -o Samples.seg
sed -i 's/_YMZ_CA//' Samples.seg
cnvkit.py export seg *bqsr.cns -o gistic.segments

./install -mode silent -agreeToLicense yes -destinationFolder /tools/GISTIC2/MATLAB_Compiler_Runtime
cd /tools/GISTIC2

#!/bin/sh
## run example GISTIC analysis
## output directory
echo --- creating output directory ---
basedir=/HSCR-raw/sunjie_work/Melanoma/result/4_work/WES/test/result
#mkdir -p $basedir 

echo --- running GISTIC ---
## input file definitions
segfile=/HSCR-raw/sunjie_work/Melanoma/result/4_work/WES/test/result/Samples.seg
#markersfile=`pwd`/examplefiles/markersfile.txt
refgenefile=`pwd`/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat
#alf=`pwd`/examplefiles/arraylistfile.txt
#cnvfile=`pwd`/examplefiles/cnvfile.txt
## call script that sets MCR environment and calls GISTIC executable 
./gistic2 -b $basedir -seg $segfile -refgene $refgenefile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme

