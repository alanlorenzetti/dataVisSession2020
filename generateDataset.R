# alorenzetti 202008

# description ####
# this script will generate/parse
# our example dataset for this session

# loading libs ####
library(tidyverse)
library(rtracklayer)

# loading datasets ####
funcat = read_tsv("dataSource/funcat.tsv") %>% 
  select(-HL, -cai, -locus_tag, -arCOG_ID, -arCOGcode, - arCOGproduct) %>% 
  rename(pfeiLocusTag = "locus_tag")
annot = rtracklayer::import("dataSource/Hsalinarum-gene-annotation-pfeiffer2019.gff3")
exp = read_tsv("dataSource/exp.tsv")

# parsing
annotCDS = annot[annot$type == "CDS",]
annotCDS$length = width(annotCDS)

# getting length dataset
df = tibble(locus_tag = annotCDS$locus_tag,
            length = annotCDS$length)

# merging length dataset with funcat dataset
funcat = left_join(funcat, df, by = "locus_tag")

# merging exp dataset with funcat dataset
final = left_join(funcat, exp, by = "locus_tag") %>% 
  rename(locus_tag = "ID",
         pfeiProduct = "protein_product",
         arCOG = "biological_class",
         mean_abundance_rna_total_TP2 = "mRNA_expression",
         mean_abundance_protein_lysate_TP2 = "protein_expression") %>% 
  select(ID,
         protein_product,
         biological_class,
         mRNA_expression,
         protein_expression,
         length,
         GC) %>% 
  mutate(mRNA_expression = log10(mRNA_expression),
         protein_expression = log10(protein_expression))

# writing final dataset
openxlsx::write.xlsx(x = final, file = "data/haloExpression.xlsx", keepNA = T)
write_tsv(x = final, path = "data/haloExpression.tsv")
