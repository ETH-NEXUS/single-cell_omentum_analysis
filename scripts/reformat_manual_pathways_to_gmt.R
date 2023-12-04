'%&%' = function(a,b) paste(a,b,sep="")
library(reshape2)
library(xlsx)
path = "gfb_omentum_2020/data/selected_genes_pathways/2021-05_manual_pathways_to_gmt/"
#suffix = "data/general/Genesets_Manual_Literature/"
#xlfiles = list.files(path %&% suffix, pattern="xlsx$", recursive = F, full.names = F)
xlfiles = list.files(path, pattern="xlsx$", recursive = F, full.names = F)
s_names = gsub("(.*)\\.xlsx$", "\\1", xlfiles)
gs_manual_lit = vector("list", length(xlfiles))

# loop over all xlsx files
for(ii in seq_along(xlfiles)){
  # 1. open xlsx file
  # ii = 1
  dd = read.xlsx2(path %&% xlfiles[ii], sheetIndex = 1, header = F)
  # 2. get pathway name
  s_name = s_names[ii]
  # 3. get gene list
  s_genes = sort(dd[, 2])
  s_genes = s_genes[s_genes != ""]
  # write *.gmt
  out = c(s_name, s_name, s_genes)
  write.table(t(out), path %&% "gmt_files/" %&% s_name %&% ".gmt",
              row.names = F, col.names = F, sep="\t", quote = FALSE)
  gs_manual_lit[[ii]] = s_genes
  names(gs_manual_lit)[ii] = s_name
}


#saveRDS(gs_manual_lit, path %&% suffix %&% "gs_manual_lit.RDS")

out = matrix("", nrow=length(xlfiles), ncol=max(sapply(gs_manual_lit, length))+2)
dim(out)
out[, 1] = s_names
head(out)
out[, 2] = s_names
for(ii in seq_along(s_names)){
  genes = gs_manual_lit[[ii]]
  out[ii, 3:(length(genes)+2)] = genes
}  
write.table(out, path %&% "gs_manual_lit.gmt", 
            row.names = F, col.names = F, sep = "\t", quote = FALSE)
