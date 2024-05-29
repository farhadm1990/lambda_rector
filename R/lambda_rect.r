#' Correcting for 16S rRNA copy-number based on Lambda phage standard.
#' @param ps a raw phyloseq object.
#' @param lamba_id the taxa name assinged to lambda phage in your tax_table of ps file
#' @param singletone_threshold number of samples in which one particular taxa should be presented to be considered non-spurious 
#' @param out_path path to output files
#' @param negative_cont If you have used negative controls, a vector of negative control sample names. 
#' @param negative_filt Wheather data should be filtered based on negative control reads
#' @param rare_depth Rarefaction depth
#' @param taxa_level Taxonomic level to which analysis should be done.
#' @param std_threshold The standard deviation threshold coefficient to find outlier samples withing each technical replicate group.
#' @param rel_abund_threshold The uper and lower end of the abundance for lambda internal standarda and Mock sequences as a two-value vector, c(lower, upper); default value is c(0.01,0.99). 
#' @export


lambda_rector <- function(ps, 
                          lambda_id = "Lambda",
                          singletone_threshold = 1, 
                          out_path = "./", 
                          negative_cont = NULL,
                          negative_filt = TRUE, 
                          rare_depth = 10000, 
                          taxa_level = "Kingdom", 
                          std_threshold = 1.48,
                          rel_abund_threshold = c(0.01,0.99)
                          ){
                          

###############################################################################################
# ps: a raw phyloseq object
# lambda_id: the taxa name you have given to lambda during creating the phyloseq object (required)
# singletone_threshold: number of samples in which one particular taxa should have presented to be considered non-spurious
# Out_path: path to the output files
# negative_cont: a vector of the barcode name of the negative control samples.
# rare_depth: rarefaction depth.
# taxa_level: taxonomic level on which your analysis will be based
# std_threshold: standard deviation threshold coefficient for detecting outlier wihitn technical replicate groups
###############################################################################################

if(
  !dir.exists(out_path)
){
  dir.create(out_path)
}

#loading packages
if(!require(pacman)){
  install.packages(pacman)
}

pacman::p_load(tidyverse, glue, tidyr, ggplot2, phyloseq, readr, ape, decontam, patchwork, reshape2, ggrepel, S4Vectors)

if (!interactive()) {
  print("This R session is not interactive.")
} else {
  print("This R session is interactive.")
}

# Checking if the count table is not in relative abundance 
if (sum(colSums(otu_table(ps)))/ncol(otu_table(ps)) == 100 ) {
  stop("You must provide me with raw count phyloseq object!")
}


################################# Internal functions #############################

out.ASV = function(phyloseq, threshold =1, binwidth = 0.01) {
  
#This function requires phyloseq, tidyverse and glue packages to be loaded. 
    if (sum(colSums(otu_table(phyloseq)))/ncol(otu_table(phyloseq)) == 100 ) {#making the relative abundance table
                    rel_abund = as(t(otu_table(phyloseq)), "matrix")
    } else if (sum(colSums(otu_table(phyloseq)))/ncol(otu_table(phyloseq)) == 1) {
                    rel_abund = as(t(otu_table(phyloseq)), "matrix")
                    } else {
                    rel_abund = as(t(apply(otu_table(phyloseq), 
                    ifelse(taxa_are_rows(phyloseq), 1,2), 
                    function(x) x/sum(x))), "matrix")  
                    } 
                      
                      
      names.single = apply(rel_abund, 1, function(x){ifelse(x == threshold, TRUE, 
                                                    ifelse(x == sum(x), TRUE, FALSE))}) %>% reshape2::melt() %>% 
                                                    filter(value == TRUE) %>% 
                                                    dplyr::select(2) %>% 
                                                    pull%>% as.vector()
                      
    if (length(names.single) == 0 ) {
          print(glue("WOW! {length(names.single)} singletones detected in this dataset"))
                 qplot.noSing = qplot(rel_abund, geom = "histogram", binwidth = binwidth, 
                  show.legend = F, main = "Frequency count of relative abundance, no singletones detected") +
         xlab ("Relative abundance in samples") + 
                   ylab("Frequency") + theme_bw()
      
    return(structure(list(qplot.noSing)))
                            
                        } else { 
                             
single.ASV = rel_abund[rownames(rel_abund) %in% names.single,]
single.ASV[single.ASV == 0] <- NA # A separate dataset for annotation of singletones on the barplot
                            
         qplot.withSing = qplot(rel_abund, geom = "histogram", binwidth = binwidth, 
        main = "Frequency count of relative abundance with singletones") +
        geom_bar(aes(single.ASV), fill = "red",  color = NA, width = binwidth)+
                       xlab ("Relative abundance in samples") + ylab("Frequency") + 
                       geom_label_repel(aes(x = 1, y =length(rel_abund)/5), 
                       label.padding =  unit(0.55, "lines"), 
                       label = glue("{length(names.single)}\n Singletones"), color = "black") + 
           theme_bw()
                            
                       qplot.rmSing = qplot(rel_abund[!rownames(rel_abund) %in% names.single, ], geom = "histogram",
                       binwidth = binwidth, main = "Frequency count of relative abundance without singletones") +
                       xlab ("Relative abundance in samples") + ylab("Frequency")+ theme_bw()
                            
                       print(glue('Oh no..! {length(names.single)} singletones detected in the dataset'))
                       return(structure(list(qplot.withSing, qplot.rmSing, unlist(names.single))) )
                    
                        }                        
    
                             
        }
# A color picker 
disc_color = function(n, method = "brewer", seed = 1990){
    if(method == "rainbow"){
        return(rainbow(n))
    } else if(method == "brewer"){
        
        
        if(n<=11){
            cols = RColorBrewer::brewer.pal(n = 11, name = "Set3")[1:n]
          

             
        } else if( n>11 && n <=74){
            sets_name  =  RColorBrewer::brewer.pal.info %>% data.frame() %>% rownames_to_column("names") %>% filter(category == "qual") %>% select(names) %>% pull
            sets_num = RColorBrewer::brewer.pal.info %>% data.frame() %>% rownames_to_column("names") %>% filter(category == "qual") %>% select(maxcolors) %>% pull 


            pooled = list()
            for(i in 1:length(sets_name)){
                pooled[[i]] <- RColorBrewer::brewer.pal(n = as.numeric(sets_num[i]), name = sets_name[i])
                
            }
            cols = do.call(c, pooled)[1:n]

        } else if( n > 74){
            sets_name  =  RColorBrewer::brewer.pal.info %>% data.frame() %>% rownames_to_column("names") %>% select(names) %>% pull
            sets_num = RColorBrewer::brewer.pal.info %>% data.frame() %>% rownames_to_column("names")  %>% select(maxcolors) %>% pull 
            pooled = list()
            for(i in 1:length(sets_name)){
                pooled[[i]] <- RColorBrewer::brewer.pal(n = as.numeric(sets_num[i]), name = sets_name[i])
                
            }
            cols = do.call(c, pooled)
            set.seed(seed)
            cols = cols[permute::shuffle(cols)][1:n]

        } else if(!is.numeric(n)){
            stop(
                "Please use a valid number!"
            )
        }
    } else if (method == "viridis"){
        cols = viridis::viridis(n = n, option = "D")
    }
return(cols)
}

# Gloomer
gloomer = function(ps = data, taxa_level = taxa_level, NArm = "TRUE"){
    rank.names = c('Kingdom','Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
    

#====================Sometimes in genus level, we might have multiple uncultured organisms, which if we want to make unique out of them for the species level it won't work====
    #since adding uncultured to uncultered is sill duplication. therefore if the taxa_level is set to species we first make a unique genus and then we go further to the speices===#

#Removing unculured Family
ps = subset_taxa(ps, !Family %in% c("uncultured", "NA", "uncategorized", "unassigend", "", " "))
    
if(taxa_level == "Species") {

    ps = subset_taxa(ps, !Genus %in% NA)#we remove genus tagged NA
tax_table(ps)[, taxa_level] <- ifelse(is.na(tax_table(ps)[, taxa_level]), paste0("unknown"), paste(tax_table(ps)[, taxa_level]))#convert NA in species into unknown
    
  physeq = tax_glom(physeq = ps, taxrank = taxa_level, NArm = NArm)
  taxdat = tax_table(physeq)[, seq_along(rank.names[1:which(rank.names == taxa_level)])]

   taxdat = taxdat[complete.cases(taxdat),] %>% as.data.frame
    otudat = otu_table(physeq)
   
#first take care of the uncultured genus
taxdat[,6] = ifelse(taxdat[,6] %in% c("uncategorized", NA, "uncultured", "unassigend", "", " "),
       paste0("[", taxdat[,length(rank.names[1:which(rank.names=="Genus")])-1], "]", "_", taxdat[,6]), taxdat[,6])
    
spec1 = taxdat[, taxa_level] %>% as.vector
spec2  = taxdat[, taxa_level] %>% as.vector

    uni  = matrix(NA, ncol = length(spec2), nrow = length(spec1))
    for(i in seq_along(spec1)){
        for(j in seq_along(spec2)){
    uni[i, j] = ifelse(spec1[i] == spec2[j] , "TRUE", "FALSE")
    }
        }

rownames(uni) <-spec1
colnames(uni) <- spec2   
uni[upper.tri(uni, diag = TRUE)] = 0 #get rid of diagonals and upper triangle

duplis = uni %>% reshape2::melt() %>% filter(value == "TRUE") 

if(dim(duplis)[[1]] > 0) {
duplis = uni %>% reshape2::melt() %>% filter(value == "TRUE") %>% dplyr::select(1) %>% unique() %>% unlist %>% as.vector
taxdat = taxdat %>% mutate( uni= ifelse(taxdat[, taxa_level] %in% duplis,
                    paste0("[", taxdat[,length(rank.names[1:which(rank.names==taxa_level)])-1], "]", "_", taxdat[,taxa_level]), taxdat[,taxa_level]))

#check if all the names are unique at species level, otherwise we will bring family instead of genus
   dupies <-  taxdat[duplicated(taxdat[,"uni"]), "uni"] 
    if(length(dupies)>0) {
        taxdat = taxdat %>% data.frame %>% mutate( uni2= ifelse(taxdat[, "uni"] %in% dupies,
                    paste0("[", taxdat[,length(rank.names[1:which(rank.names==taxa_level)])-2], "]", "_", taxdat[,"uni"]), taxdat[,"uni"]))
        
        taxdat[, taxa_level] = taxdat[, "uni2"]
        taxdat[, "uni"] <- NULL
        taxdat[, "uni2"] <- NULL
        taxdat <- as(taxdat, "matrix")   
        rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
        rownames(taxdat) <- taxdat[, taxa_level]
        taxdat <- tax_table(taxdat)
        taxa_names(physeq) <- taxa_names(taxdat)
        tax_table(physeq) <- taxdat
        otu_table(physeq) <- otudat
        
    }
    else 
    {
        
taxdat[, taxa_level] = taxdat[, "uni"]
taxdat[, "uni"] <- NULL
taxdat <- as(taxdat, "matrix")   
rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[, taxa_level]
taxdat <- tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat
           }
    
} else {
    
taxdat <- as.matrix(taxdat) 
taxdat <- tax_table(taxdat)
rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[, taxa_level]
taxdat <- tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat
    
}
       
    
#==========================================# 
} else if (taxa_level == "Genus") {
    
    physeq = tax_glom(physeq = ps, taxrank = taxa_level, NArm = NArm)
    taxdat = tax_table(physeq)[, seq_along(rank.names[1:which(rank.names == taxa_level)])]
    
   taxdat = taxdat[complete.cases(taxdat),] %>% as.data.frame
    otudat = otu_table(physeq)
    
# take care of the uncultured genus
taxdat[,6] = ifelse(taxdat[,6] %in% c("uncategorized", NA, "uncultured", "unassigend", "", " "),
       paste0("[", taxdat[,length(rank.names[1:which(rank.names==taxa_level)])-1], "]", "_", taxdat[,taxa_level]), taxdat[,taxa_level])
    
gen1 = taxdat[, taxa_level] %>% as.vector
gen2  = taxdat[, taxa_level] %>% as.vector

    uni  = matrix(NA, ncol = length(gen2), nrow = length(gen1))
    for(i in seq_along(gen1)){
        for(j in seq_along(gen2)){
    uni[i, j] = ifelse(gen1[i] == gen2[j] , "TRUE", "FALSE")
    }
        }

rownames(uni) <-gen1
colnames(uni) <- gen2   
uni[upper.tri(uni, diag = TRUE)] = 0 #get rid of diagonals and upper triangle

duplis = uni %>% reshape2::melt() %>% filter(value == "TRUE")

        if(dim(duplis)[[1]] > 0){#if there is not duplications, we can simply use the taxa names as the row name
    
        duplis = uni %>% reshape2::melt() %>% filter(value == "TRUE") %>% dplyr::select(1)%>% unique() %>% unlist %>% as.vector
        taxdat = taxdat %>% mutate( uni= ifelse(taxdat[, taxa_level] %in% duplis, 
                    paste0("[", taxdat[,length(rank.names[1:which(rank.names==taxa_level)])-1], "]", "_", taxdat[,taxa_level]), taxdat[,taxa_level]))
    
        taxdat[, taxa_level] = taxdat[, "uni"]
        taxdat[, "uni"] <- NULL

        taxdat <- as(taxdat, "matrix")
 
        rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
        rownames(taxdat) <- taxdat[taxdat[,taxa_level] %in% rownames(otudat), taxa_level]
        taxdat <- as.matrix(taxdat) 
        taxdat <- tax_table(taxdat)
        taxa_names(physeq) <- taxa_names(taxdat)
        tax_table(physeq) <- taxdat
        otu_table(physeq) <- otudat
 
        } else {

        taxdat <- as.matrix(taxdat) 
        taxdat <- tax_table(taxdat)
        rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
        rownames(taxdat) <- taxdat[, taxa_level]
        taxdat <- tax_table(taxdat)
        taxa_names(physeq) <- taxa_names(taxdat)
        tax_table(physeq) <- taxdat
        otu_table(physeq) <- otudat
       }   
    
} else {
    
    
physeq = tax_glom(physeq = ps, taxrank = taxa_level, NArm = TRUE)
    taxdat = tax_table(physeq)[, seq_along(rank.names[1:which(rank.names == taxa_level)])]
    
taxdat = taxdat[complete.cases(taxdat),] %>% as.data.frame
otudat = otu_table(physeq)
    
spec1 = taxdat[, taxa_level] %>% as.vector
spec2  = taxdat[, taxa_level] %>% as.vector

    uni  = matrix(NA, ncol = length(spec2), nrow = length(spec1))
    for(i in seq_along(spec1)){
        for(j in seq_along(spec2)){
    uni[i, j] = ifelse(spec1[i] == spec2[j] , "TRUE", "FALSE")
    }
        }

rownames(uni) <-spec1
colnames(uni) <- spec2   
uni[upper.tri(uni, diag = TRUE)] = 0 #get rid of diagonals and upper triangle

duplis = uni %>% reshape2::melt() %>% filter(value == "TRUE")

if(dim(duplis)[[1]] > 0){#if there is not duplications, we can simply use the taxa names as the row name
    
    duplis = uni %>% reshape2::melt() %>% filter(value == "TRUE") %>% dplyr::select(1)%>% unique() %>% unlist %>% as.vector
taxdat = taxdat %>% mutate( uni= ifelse(taxdat[, taxa_level] %in% duplis, 
                    paste(taxdat[,length(rank.names[1:which(rank.names==taxa_level)])-1], "_", taxdat[,taxa_level]), taxdat[,taxa_level]))
    
taxdat[, taxa_level] = taxdat[, "uni"]
taxdat[, "uni"] <- NULL
taxdat <- as.matrix(taxdat)   
rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[, taxa_level]
taxdat <- tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat
} else {

taxdat <- as.matrix(taxdat) 
taxdat <- tax_table(taxdat)
rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[, taxa_level]
taxdat <- tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat
}
#ps = phyloseq(otu_table(otudat, taxa_are_rows = T), tax_table(as.matrix(taxdat)), sample_data(physeq))
 

}
return(physeq) 
    }




# Filtering
## removing zero counts of lambda and mock

ps <- subset_taxa(ps, !Kingdom %in% c("NA", "unassigned", "Archaea", "Eukaryota") & !Phylum %in% c("NA", "unassigned", "uncultured") & !is.na(Phylum))

suspect <- character(0)
if(negative_filt == TRUE){
        # decontatamination by decontam
        sample_data(ps)$is_neg <- as.logical(ifelse(sample_names(ps) %in% negative_cont, TRUE, FALSE))

        contam <<- isContaminant(ps, method = "prevalence", neg = "is_neg", threshold = 0.1) %>% filter(contaminant == TRUE) %>% rownames()

        negative_otu <<-  rownames(tax_table(ps)[colnames(ps@otu_table) %in% negative_cont,])
        
        suspect <- c(contam, negative_otu) %>% as.character()

       suspect <<- data.frame(suspected_otu = suspect, Species = tax_table(ps)[taxa_names(ps) %in% suspect,7]) %>% filter(Species !="Lambda") %>% select(suspected_otu) %>% pull() %>% as.character()
}


message(glue("\n\nHere is the list for  {length(suspect)} detected suspected OTUs: \n"))
        
df = data.frame(suspected_otu = suspect, Species = tax_table(ps)[taxa_names(ps) %in% suspect,7])

print(df)
flush.console()


consent1 = readline(prompt = "Do you want me to remove them? (y/n) ")

while(!consent1 %in% c("y", "n")){
  consent1 <- readline(prompt = "Please enter 'y' for yes and 'n' for no: ")
  
}

if(consent1 == "y"){
  
  print(suspect)
  ps <- phyloseq::subset_taxa(ps, !taxa_names(ps) %in% suspect)
  print("Suspected taxa have been removed. Now I go to the next step! ") 
}

#removing negative control and empty samples
negative_cont <<- negative_cont
ps <- subset_samples(ps, !sample_names(ps) %in% negative_cont & sample_data(ps)$lambda_ng_ul > 0 & sample_data(ps)$mock_ng_ul > 0)

# Prevalance filtering
prev <- apply(ps@otu_table, 1, function(x){sum(x>0)})
df <- data.frame(asv_prev = prev, 
                tax_abund = taxa_sums(ps),
                tax_table(ps))

df = plyr::ddply(df, "Phylum", function(x){
              cbind(means = round(mean(x$asv_prev), 2),
                    sums = round(sum(x$asv_prev), 2))
}) %>% mutate(sanity = ifelse(means == sums, TRUE, FALSE))  

sus <<- df %>% filter(sanity == TRUE) %>% select(Phylum) %>% pull()

print(df)
flush.console()

      if(!is.null(sus)){

        consent2 <- readline(prompt = paste0(sus, " phyla are suspected, please double check and tell if I shoudl remove them? (y/n) "))
        
        while(!consent2 %in% c("y", "n")){
          consent2 <- readline(prompt = "Please enter 'y' for yes and 'n' for no: ")
        }

        if(consent2 == "y"){

        ps <- subset_taxa(ps, !Phylum %in% sus)


      }
      }

# removing negative column

sample_data(ps)$is_neg <- NULL
# Removing singletones                      
if(length((single.test = out.ASV(phyloseq = ps, threshold = singletone_threshold, binwidth = 0.001))) > 1){
 

sings <<- single.test[[3]]


  const3 <- readline(prompt = glue("I found {length(sings)} singletones in your dataset, do you want to see their disttribution plot and their names? (y/n) "))
  while(!const3 %in% c("y", "n")){
    const3 <- readline("You must enter y for Yes and n for No: ")
  }
  if(const3 == "y"){
    singletones %>% data.frame()
    flush.console()
    #to show the plot with singletones
    cat("Here is the plot with singletones! ")
    single.test[[1]]
    
    
    message("Here is the plot without singletones. ")
    single.test[[2]]

    ps <- subset_taxa(ps, !taxa_names(ps) %in% sings)
    
  }
}





# Rarefaction
print(glue("Rarefying data to {rare_depth} depth! \n"))

ps_rar <<- rarefy_even_depth(ps, sample.size = rare_depth, replace = FALSE)

rel <- transform_sample_counts(ps_rar, function(x){x/sum(x)})

#removing taxa with 0.1% and more thatn 98% abundance. 

filt <- prune_taxa( taxa_sums(rel) > rel_abund_threshold[1] &  taxa_sums(rel) < rel_abund_threshold[2], ps_rar)

print(glue("Agglomerating the taxa to {taxa_level} taxonomy level! \n"))
glomed_ps <- gloomer(filt, taxa_level = taxa_level, NArm = T)



print("Transforming sample counts to relative abundance!\n")




ps_rel <- transform_sample_counts(glomed_ps, function(x){x/sum(x)})


# Removing suspected samples within each lambda concentration group

std_dev_threshold <- std_threshold


sus_samp <<- psmelt(ps_rel) %>% select(Kingdom, Sample, lambda_ng_ul, mock_ng_ul, Abundance) %>% filter(Abundance <=rel_abund_threshold[2], Abundance > rel_abund_threshold[1], mock_ng_ul > 0.002, lambda_ng_ul > 0, Kingdom == lambda_id) %>% 
 group_by( lambda_ng_ul, mock_ng_ul, Kingdom) %>% 
 arrange(desc(lambda_ng_ul)) %>% 
 mutate( mean_abundance = mean(Abundance, na.rm = TRUE),
          sd_abundance = sd(Abundance, na.rm = TRUE), 
          is_outlier = abs(Abundance - mean_abundance) > (std_dev_threshold * sd_abundance)
                ) %>% 
        data.frame() %>% 
        filter(is_outlier) %>% 
        select(Sample)  %>% 
        pull() %>% as.character()

print(glue("{sus_samp} samples seem to be suspicous and I am going to remove them!\n But first take a look at them in {out_path}"))

p_with = psmelt(ps_rel) %>% select(Kingdom, Sample, lambda_ng_ul, mock_ng_ul, Abundance) %>% filter(Abundance <=rel_abund_threshold[2], Abundance > rel_abund_threshold[1], mock_ng_ul > 0.002, lambda_ng_ul > 0 ) %>% group_by(Sample, lambda_ng_ul, mock_ng_ul, Kingdom) %>% summarise(mean = sum(Abundance), .groups = "drop")%>% mutate(Kingdom = ifelse(Kingdom ==  lambda_id , "Lambda", "Sample")) %>% ggplot() + 
geom_col(aes(x = Sample, y = mean, fill = Kingdom), position = "stack", show.legend = T) + 
facet_grid(mock_ng_ul ~ lambda_ng_ul , scales = "free")  + theme_bw()+ 
theme(axis.text.x = element_text(angle = 45, hjust = 1 )) + 
ggtitle("With suspected samples")

ggsave(plot = p_with, paste0(out_path, "/plot_with_bad_samples.jpeg"), height = 12, width = 10, device = "jpeg", dpi = 300)

print("Here is how your data looks like wihtout those samples.")

p_without <- psmelt(ps_rel) %>% select(Kingdom, Sample, lambda_ng_ul, mock_ng_ul, Abundance) %>% filter(Abundance <=rel_abund_threshold[2], Abundance > rel_abund_threshold[1], mock_ng_ul > 0.002, lambda_ng_ul > 0, !Sample %in% sus_samp ) %>% group_by(Sample, lambda_ng_ul, mock_ng_ul, Kingdom) %>% summarise(mean = sum(Abundance), .groups = "drop")%>% mutate(Kingdom = ifelse(Kingdom == lambda_id, "Lambda", "Sample")) %>% ggplot() + 
geom_col(aes(x = Sample, y = mean, fill = Kingdom), position = "stack", show.legend = T) + 
facet_grid(mock_ng_ul ~ lambda_ng_ul , scales = "free")  + theme_bw()+ 
theme(axis.text.x = element_text(angle = 45, hjust = 1 )) + 
ggtitle("Without suspected samples")

ggsave(plot = p_without, paste0(out_path, "/plot_without_bad_samples.jpeg"), height = 12, width = 10, device = "jpeg", dpi = 300)

consent4 = readline("Should I remove those samples? y/n")
while(!consent4 %in% c("y", "n")){
  consent4 = readline("Please enter y for yes or n for no! ")
}

glomed_ps <<- glomed_ps
ps_rel <<- ps_rel

# Removing suspected samples upon user's consent
if(consent4 == "y"){

  glomed_ps <- subset_samples(glomed_ps, !sample_names(glomed_ps) %in% sus_samp)
  ps_rel  <- subset_samples(ps_rel, !sample_names(ps_rel) %in% sus_samp)
}



# creating copy number ps

lambda_matrix = as(
  otu_table(glomed_ps)[taxa_names(glomed_ps) == lambda_id,] %>% 
  t(), 
  "matrix"
  )

samp_matrix = as(
  otu_table(glomed_ps)[taxa_names(glomed_ps) != lambda_id,] %>% t(), 
  "matrix"
  )

copy_loaded = as(
  sample_data(glomed_ps) %>% 
  data.frame() %>% 
  select(loaded_copy_lambda) %>% 
  dplyr::rename( "Lambda" = loaded_copy_lambda), 
  "matrix")


temp <- apply(samp_matrix, 2, function(x){round(x * copy_loaded/ lambda_matrix)}) 

rownames(temp) <- rownames(samp_matrix)


temp <- as(cbind(temp, copy_loaded), "matrix")
temp <-  na.omit(temp)
temp <- as(temp, "matrix")

copy_corrected_ps <- glomed_ps




otu_table(copy_corrected_ps) <- otu_table(t(temp), taxa_are_row = TRUE)

# Visualizing the results relative abundance


melt_df <- psmelt(ps_rel) %>% 
filter(lambda_ng_ul > 0, mock_ng_ul > 0.002) %>% 
select(contains(taxa_level), Sample, Abundance, lambda_ng_ul, mock_ng_ul) %>% 
group_by(.data[[taxa_level]], Sample, lambda_ng_ul, mock_ng_ul) %>% 
summarise(abund = sum(Abundance), .groups = "drop") 



melt_df$mock_ng_ul <- factor(
  melt_df$mock_ng_ul, 
  levels = c(20,2,0.2,0.02), 
  labels = unique(melt_df$mock_ng_ul)
  )

melt_df$lambda_ng_ul <- factor(
  melt_df$lambda_ng_ul, 
  levels = c(1e-04, 1e-05, 1e-06), 
  labels = unique(melt_df$lambda_ng_ul)
  )


lambda_con_1 <- unique(melt_df$lambda_ng_ul)[1]
lambda_con_2 <- unique(melt_df$lambda_ng_ul)[2]
lambda_con_3 <- unique(melt_df$lambda_ng_ul)[3]

df1 = melt_df %>% filter(lambda_ng_ul == lambda_con_1)
df2 = melt_df %>% filter(lambda_ng_ul == lambda_con_2)
df3 = melt_df %>% filter(lambda_ng_ul == lambda_con_3)

n = length(unique(melt_df[[taxa_level]]))
set.seed(1990)

#cols = sample(colors(), size = n, replace = F)
cols = disc_color(n = n)

p1 = df1 %>% ggplot() +
    geom_col(
      aes(x = Sample, 
      y = abund, 
      fill = .data[[taxa_level]]), 
      show.legend = F
      ) + 
    # facet_grid(mock_ng_ul ~ lambda_ng_ul, scales = "free", shrink = TRUE) +
    facet_grid(
      glue("Lambda con: {lambda_ng_ul} ng") ~ glue("Mock con: {mock_ng_ul} ng"), 
      scales =  "free") + 
    theme_bw( ) +
    theme(axis.text.x = element_text(angle = 0), strip.text = element_text(color = "white", face = "bold"), strip.background = element_rect(fill = "#5e5e5e"))  +
    scale_fill_manual(values = cols) +
    labs(y = "", x = "") +  # Label for color legend
    guides(fill = guide_legend(ncol = 2, override.aes = list(size = 4, alpha = 0.7))) 

p2 = df2 %>% ggplot() +
    geom_col(aes(x = Sample, y = abund, fill = .data[[taxa_level]])) + 
    # facet_grid(mock_ng_ul ~ lambda_ng_ul, scales = "free", shrink = TRUE) +
    facet_grid(glue("Lambda con: {lambda_ng_ul} ng")~ glue("Mock con: {mock_ng_ul} ng"), scales =  "free") + 
    theme_bw( ) +
    theme(axis.text.x = element_text(angle = 0), strip.text = element_text(color = "white", face = "bold"), strip.background = element_rect(fill = "#5e5e5e"))  +
    scale_fill_manual(values = cols) +
    labs(fill = taxa_level, y = "Relative aboundance", x = "") +  # Label for color legend
    guides(fill = guide_legend(ncol = 2, override.aes = list(size = 4, alpha = 0.7))) 


p3 = df3 %>% ggplot() +
    geom_col(aes(x = Sample, y = abund, fill = .data[[taxa_level]]), show.legend = F) + 
    # facet_grid(mock_ng_ul ~ lambda_ng_ul, scales = "free", shrink = TRUE) +
    facet_grid(glue("Lambda con: {lambda_ng_ul} ng")~ glue("Mock con: {mock_ng_ul} ng"), scales =  "free") + 
    theme_bw( ) +
    theme(axis.text.x = element_text(angle = 0), strip.text = element_text(color = "white", face = "bold"), strip.background = element_rect(fill = "#5e5e5e"))  +
    scale_fill_manual(values =cols) +
    labs(y = "") +  # Label for color legend
      guides(fill = guide_legend(ncol = 2, override.aes = list(size = 4, alpha = 0.7))) 

pl_rel <- p1 + p2 + p3 + plot_layout(nrow = 3, widths = 10) + plot_annotation("Proportion of reads for loaded sample and Lambda in different mock and lambda concentrations.")

ggsave(plot =pl_rel, paste0(out_path,"/", taxa_level, "_relative.jpeg"), width = 15, height = 10, dpi =300)


# Copy_corrected with  Lambda

melt_df <- psmelt(copy_corrected_ps) %>% 
filter(lambda_ng_ul > 0, mock_ng_ul > 0.002) %>% 
select(contains(taxa_level), Sample, Abundance, lambda_ng_ul, mock_ng_ul) %>% 
group_by(.data[[taxa_level]], Sample, lambda_ng_ul, mock_ng_ul)  %>% 
summarise(abund = mean(Abundance), .groups = "drop")


melt_df$mock_ng_ul <- factor(
  melt_df$mock_ng_ul, 
  levels = c(20,2,0.2,0.02), 
  labels = unique(melt_df$mock_ng_ul)
  )

melt_df$lambda_ng_ul <- factor(
  melt_df$lambda_ng_ul, 
  levels = c(1e-04, 1e-05, 1e-06), 
  labels = unique(melt_df$lambda_ng_ul)
  )


lambda_con_1 <- unique(melt_df$lambda_ng_ul)[1]
lambda_con_2 <- unique(melt_df$lambda_ng_ul)[2]
lambda_con_3 <- unique(melt_df$lambda_ng_ul)[3]

df1 = melt_df %>% filter(lambda_ng_ul == lambda_con_1)
df2 = melt_df %>% filter(lambda_ng_ul == lambda_con_2)
df3 = melt_df %>% filter(lambda_ng_ul == lambda_con_3)

n = length(unique(melt_df[[taxa_level]]))
set.seed(1990)

# cols = sample(colors(), size = n, replace = F)
cols = disc_color(n = n)

p1 = df1 %>% ggplot() +
    geom_col(aes(x = Sample, y = abund, fill = .data[[taxa_level]]), show.legend = F) + 
    # facet_grid(mock_ng_ul ~ lambda_ng_ul, scales = "free", shrink = TRUE) +
    facet_grid(glue("Lambda con: {lambda_ng_ul} ng")~ glue("Mock con: {mock_ng_ul} ng"), scales =  "free") + 
    theme_bw( ) +
    theme(axis.text.x = element_text(angle = 0), strip.text = element_text(color = "white", face = "bold"), strip.background = element_rect(fill = "#5e5e5e"))  +
    scale_fill_manual(values = cols) +
    labs(y = "", x = "") +  # Label for color legend
    guides(fill = guide_legend(ncol = 2, override.aes = list(size = 4, alpha = 0.7))) 

p2 = df2 %>% ggplot() +
    geom_col(aes(x = Sample, y = abund, fill = .data[[taxa_level]])) + 
    # facet_grid(mock_ng_ul ~ lambda_ng_ul, scales = "free", shrink = TRUE) +
    facet_grid(glue("Lambda con: {lambda_ng_ul} ng")~ glue("Mock con: {mock_ng_ul} ng"), scales =  "free") + 
    theme_bw( ) +
    theme(axis.text.x = element_text(angle = 0), strip.text = element_text(color = "white", face = "bold"), strip.background = element_rect(fill = "#5e5e5e"))  +
    scale_fill_manual(values = cols) +
    labs(fill = taxa_level, y = "16S rRNA Copy-number of loaded Lambda and Sample", x = "") +  # Label for color legend
    guides(fill = guide_legend(ncol = 2, override.aes = list(size = 4, alpha = 0.7))) 


p3 = df3 %>% ggplot() +
    geom_col(aes(x = Sample, y = abund, fill = .data[[taxa_level]]), show.legend = F) + 
    # facet_grid(mock_ng_ul ~ lambda_ng_ul, scales = "free", shrink = TRUE) +
    facet_grid(glue("Lambda con: {lambda_ng_ul} ng")~ glue("Mock con: {mock_ng_ul} ng"), scales =  "free") + 
    theme_bw( ) +
    theme(axis.text.x = element_text(angle = 0), strip.text = element_text(color = "white", face = "bold"), strip.background = element_rect(fill = "#5e5e5e"))  +
    scale_fill_manual(values = cols) +
    labs( y = "") +  # Label for color legend
      guides(fill = guide_legend(ncol = 2, override.aes = list(size = 4, alpha = 0.7))) 

pl_cp = p1 + p2 + p3 + plot_layout(nrow = 3, widths = 10) + plot_annotation("16S rRNA Copy-number for loaded sample and Lambda in different mock and lambda concentrations.")

ggsave(plot =pl_cp, paste0(out_path,"/", taxa_level, "_copy_number.jpeg"), width = 15, height = 10, dpi =300)

# Copy_corrected without  Lambda

melt_df <- psmelt(copy_corrected_ps) %>% filter(lambda_ng_ul > 0, mock_ng_ul > 0.002, OTU !="Lambda") %>% select(contains(taxa_level), Sample, Abundance, lambda_ng_ul, mock_ng_ul) %>% 
group_by(.data[[taxa_level]], Sample, lambda_ng_ul, mock_ng_ul)  %>% summarise(abund = mean(Abundance), .groups = "drop")

melt_df$mock_ng_ul <- factor(
  melt_df$mock_ng_ul, 
  levels = c(20,2,0.2,0.02), 
  labels = unique(melt_df$mock_ng_ul)
  )

melt_df$lambda_ng_ul <- factor(
  melt_df$lambda_ng_ul, 
  levels = c(1e-04, 1e-05, 1e-06), 
  labels = unique(melt_df$lambda_ng_ul)
  )



lambda_con_1 <- unique(melt_df$lambda_ng_ul)[1]
lambda_con_2 <- unique(melt_df$lambda_ng_ul)[2]
lambda_con_3 <- unique(melt_df$lambda_ng_ul)[3]

df1 = melt_df %>% filter(lambda_ng_ul == lambda_con_1)
df2 = melt_df %>% filter(lambda_ng_ul == lambda_con_2)
df3 = melt_df %>% filter(lambda_ng_ul == lambda_con_3)

n = length(unique(melt_df[[taxa_level]]))
set.seed(1990)
# cols = sample(colors(), size = n, replace = F)
cols = disc_color(n = n)

p1 = df1 %>% ggplot() +
    geom_col(aes(x = Sample, y = abund, fill = .data[[taxa_level]]), show.legend = F) + 
    # facet_grid(mock_ng_ul ~ lambda_ng_ul, scales = "free", shrink = TRUE) +
    facet_grid(glue("Lambda con: {lambda_ng_ul} ng")~ glue("Mock con: {mock_ng_ul} ng"), scales =  "free") + 
    theme_bw( ) +
    theme(axis.text.x = element_text(angle = 0), strip.text = element_text(color = "white", face = "bold"), strip.background = element_rect(fill = "#5e5e5e"))  +
    scale_fill_manual(values = cols) +
    labs(y = "", x = "") +  # Label for color legend
   guides(fill = guide_legend(ncol = 2, override.aes = list(size = 4, alpha = 0.7))) 

p2 = df2 %>% ggplot() +
    geom_col(aes(x = Sample, y = abund, fill = .data[[taxa_level]])) + 
    # facet_grid(mock_ng_ul ~ lambda_ng_ul, scales = "free", shrink = TRUE) +
    facet_grid(glue("Lambda con: {lambda_ng_ul} ng")~ glue("Mock con: {mock_ng_ul} ng"), scales =  "free") + 
    theme_bw( ) +
    theme(axis.text.x = element_text(angle = 0), strip.text = element_text(color = "white", face = "bold"), strip.background = element_rect(fill = "#5e5e5e"))  +
    scale_fill_manual(values = cols) +
    labs(fill = taxa_level, y = "16S rRNA Copy-number of loaded Sample", x = "") +  # Label for color legend
    guides(fill = guide_legend(ncol = 2, override.aes = list(size = 4, alpha = 0.7))) 

p3 = df3 %>% ggplot() +
    geom_col(aes(x = Sample, y = abund, fill = .data[[taxa_level]]), show.legend = F) + 
    # facet_grid(mock_ng_ul ~ lambda_ng_ul, scales = "free", shrink = TRUE) +
    facet_grid(glue("Lambda con: {lambda_ng_ul} ng")~ glue("Mock con: {mock_ng_ul} ng"), scales =  "free") + 
    theme_bw( ) +
    theme(axis.text.x = element_text(angle = 0), strip.text = element_text(color = "white", face = "bold"), strip.background = element_rect(fill = "#5e5e5e"))  +
    scale_fill_manual(values = cols) +
    labs( y = "") +  # Label for color legend
      guides(fill = guide_legend(ncol = 2, override.aes = list(size = 4, alpha = 0.7))) 

pl_cp = p1 + p2 + p3 + plot_layout(nrow = 3, widths = 10) + plot_annotation("16S rRNA Copy-number for loaded sample in different mock and lambda concentrations.")

ggsave(plot =pl_cp, paste0(out_path,"/", taxa_level, "_copy_number_noLambda.jpeg"), width = 15, height = 10, dpi =300)



  structure(list(copy_corrected_ps = copy_corrected_ps, ps_gloomed = glomed_ps, ps_rel = ps_rel, raw_ps = ps, copy_corrected_matrix = temp))

}
