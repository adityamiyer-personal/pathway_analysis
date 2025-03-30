#use this for any pathway to extract pathway-specific DEGs from the comparison where enriched [empty lists discarded]
mylist <- map(split(list_rbind(df33, names_to = "group"), f = list(list_rbind(df33, names_to = "group")$group, 
                                               list_rbind(df33, names_to = "group")$species)),function(x) {x %>%
    filter(ID == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    pull(geneID)})
mylist <- purrr::compact(mylist)
matrix = from_list(mylist)
matrix$gene_name = rownames(matrix)

#Upset plot with gene names- open with Zoom to look at the gene names
#look at SuperExactTest for upset plots
ComplexUpset::upset(
    matrix,
    intersect=colnames(matrix)[-length(colnames(matrix))],
    base_annotations=list(
        'Intersection size'=(
            intersection_size(
                bar_number_threshold=1,
                color='grey9',
                fill='grey80'
            )
            + ggfittext::geom_bar_text(
                mapping=aes(label=gene_name),
                min.size=0,
                position='stack',
                contrast=FALSE,
                vjust=1.1,
            )
        )
    ),
    width_ratio=0.15,
    height_ratio=1/4
)

cat(paste(matrix %>% 
  remove_rownames() %>%
  column_to_rownames("gene_name") %>%
  #filter_at(vars(-gene_name), all_vars(isTRUE(.))) %>%
  #filter(colnames(matrix)("-gene_name"))
  filter_all(all_vars(. == TRUE)) %>% 
  rownames()), sep = "\n")

pathway_genes = matrix %>% 
  remove_rownames() %>%
  column_to_rownames("gene_name") %>%
  #filter_at(vars(-gene_name), all_vars(isTRUE(.))) %>%
  #filter(colnames(matrix)("-gene_name"))
  filter_all(all_vars(. == TRUE)) %>% 
  rownames()
