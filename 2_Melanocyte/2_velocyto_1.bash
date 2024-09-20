# Velocyto data preparation using Linux commands
#####################################################################################
# 10X
for donor in D5_Nevus D5_Malignant_Nevus D5_Melanoma;do
  velocyto run10x -m /velocyto_data/repetitive_element/hs38_rmsk.gtf /2_Melanocyte/2_velocyto/${donor} /gene_annotation/Homo_sapiens.GRCh38.99.gtf
done

# BD
for donor in D1 D2 D3 D4;do
  cd /2_Melanocyte/2_velocyto/${donor}/
  samtools view -HS /${donor}/${donor}_R_1_final.BAM > head.txt
  samtools view /${donor}/${donor}_R_1_final.BAM | grep "MA:Z:*" | sed "s/MA:Z:/UB:Z:/" > temp.sam
  cat head.txt temp.sam | samtools view -Sb > changed.bam
  velocyto run -e ${donor} -o /2_Melanocyte/2_velocyto/${donor} -m /velocyto_data/repetitive_element/hs38_rmsk.gtf changed.bam /gene_annotation/Homo_sapiens.GRCh38.99.gtf
done



