# Get gene annotation from Homo_sapiens.GRCh38.94.gtf
1 Processed by linux  
```
awk -F "\t" '$3 == "gene" { print $9 }' Homo_sapiens.GRCh38.94.gtf |awk -F "; " '{print $3 "\t" $5}' > type.txt
```   
2 Processed by r  
```
gene2type = read.table( "type.txt" ,stringsAsFactors = F )  
gene2type = gene2type[,c(2,4)]  
lincRNA = gene2type[gene2type[,2]=="lincRNA",]
```
