# alorenzetti 202008

# description ####
# this script will generate/parse
# our example dataset for this session

# loading libs ####
library(tidyverse)
library(rtracklayer)

# loading and removing duplicate entries of ppi file
ppi = read_tsv("dataSource/ppi.sif", col_names = F) %>% 
  distinct()

write_tsv(ppi, path = "data/ppi.sif", col_names = F)

# loading and filtering cytoscape classification dataset ####
ppiFunCat = read_csv("dataSource/ppiFunCat.csv") %>% 
  select(ID = name,
         label,
         is_bait,
         product = productPfeiffer2019,
         class = KEGGpathway) %>% 
  mutate(label = case_when(is.na(label) ~ ID,
                           TRUE ~ as.character(label)),
         is_bait = case_when(is.na(is_bait) ~ FALSE,
                             TRUE ~ is_bait),
         product = case_when(is.na(product) ~ "Product Unknown",
                             TRUE ~ as.character(product)),
         class = case_when(is.na(class) ~ "Class Unknown",
                           class == "hal03022 Basal transcription factors" ~ "Transcription",
                           class == "tc" ~ "Transcription",
                           class == "hal03010 Ribosome" ~ "Translation",
                           class == "tl" ~ "Translation",
                             TRUE ~ as.character(class))) %>% 
  mutate(class = sub(pattern = "^hal..... ", replacement = "", x = class)) %>% 
  mutate(class = sub(pattern = "; .*$", replacement = "", x = class))



# saving cytoscape dataset
write_tsv(x = ppiFunCat, path = "data/ppiFunCat.tsv")

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
