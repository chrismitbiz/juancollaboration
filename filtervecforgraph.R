filtervecforgraph <- c("D_5__uncultured bacterium",                       "D_5__Candidatus Udaeobacter",                    
"D_5__Chthoniobacter"     ,                        "D_5__Candidatus Nitrocosmicus",                  
 "D_5__Candidatus Nitrososphaera"  ,                "D_5__Ramlibacter"            ,                   
"D_5__Burkholderia-Caballeronia-Paraburkholderia", "D_5__Sphingomonas"            ,                  
"D_5__Pseudolabrys"              ,                 "D_5__Steroidobacter"          ,                  
"D_5__Nitrosospira"               ,                "D_5__uncultured"              ,                  
"D_5__Bradyrhizobium"             ,                "D_5__Devosia"                 ,                  
"D_5__Microvirga"                  ,               "D_5__Skermanella"             ,                  
"D_5__Acidibacter"                  ,              "D_5__uncultured soil bacterium",                 
"D_5__uncultured planctomycete"      ,             "D_5__uncultured Planctomycetales bacterium" ,    
"D_5__Pirellula"                    ,              "D_5__Gemmatimonas"             ,                 
"D_5__Bacillus"                     ,              "D_5__Tumebacillus"            ,                  
"D_5__uncultured Chloroflexi bacterium",           "D_5__Flavisolibacter"         ,                  
"D_5__Streptomyces"                 ,              "D_5__Catenulispora"           ,                  
"D_5__Blastococcus"                 ,              "D_5__Geodermatophilus"        ,                  
"D_5__Pseudonocardia"               ,              "D_5__RB41")

write.csv(filtervecforgraph, "filtervecforgraph.csv")


#barplot
barplot.rel <- phyloseq::transform_sample_counts(physeqBSoil1, function(x) x / sum(x) )

figvec <- read.csv("~/Documents/Study/LaTrobe/Research/phD/Juan_project/juancollaboration/filtervecforgraph.csv", header = TRUE)
subs <- phyloseq::subset_taxa(barplot.rel, Genus %in% figvec$x)
subs <- psmelt(subs)
barplot.relsoil1 <- subs  %>% ggbarplot(data = ., x = "Genus", y = "Abundance", add = "mean_se", facet.by = c("Genus") , fill = "Treatments", width = 0.5, xlab = "", position = position_dodge(0.7), scales = "free", ylab = "Relative abundances") + theme(axis.title.x=element_blank(),             axis.text.x=element_blank(), axis.ticks.x=element_blank()) # theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
barplot.relsoil1 



## Vertisol - Sign. diff. on genus level (Wilcoxon test *p* < 0.001 holm corrected)
```{r vertsigdiff, echo=FALSE, fig.width=15, message=FALSE, warning=FALSE, cache=FALSE}


#Soil2 data
physeqBSoil2 <-  prune_samples(sample_data(physeqB.flt)$Soiltypes == "Vertisol", physeqB.flt)

#Soil1
#wilcoxon test
barplot.rel <- microbiome::transform(physeqBSoil2, "clr")
barplot.rel.melt <- psmelt(barplot.rel) %>%  arrange(desc(Phylum))
wilcoxotest <-  barplot.rel.melt %>% compare_means(Abundance ~ Treatments, data = ., group.by = "Genus", method = "wilcox.test", paired = FALSE, ref.group = "Control", symnum.args = list(), p.adjust.method = "holm") %>% as.data.frame() %>% filter(p.adj <= 0.001)
kable(wilcoxotest %>% filter(Genus != ""),booktabs=TRUE) %>% 
  kable_styling(latex_options="scale_down")


#barplot
vec <- unique(wilcoxotest %>% filter(Genus != "") %>% dplyr::select(Genus))  
barplot.rel <- phyloseq::transform_sample_counts(physeqBSoil2, function(x) x / sum(x) )
subs <- subset_taxa(barplot.rel, Genus %in% vec$Genus )
subs <- psmelt(subs)
barplot.relsoil2 <- subs  %>% ggbarplot(data = ., x = "Genus", y = "Abundance", add = "mean_se", facet.by = c("Genus") , fill = "Treatments", width = 0.5, xlab = "", position = position_dodge(0.7), scales = "free") # theme(axis.text.x = element_text(angle = 45, hjust = 1))

barplot.relsoil2 
```

