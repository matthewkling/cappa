
library(tidyverse)
library(ape)
library(raster)
library(phytools)
library(ggtree)

select <- dplyr::select

indir <- "e:/phycon/002_sxt/pred_s25km"
outdir <- "e:/phycon/results/pred_s25km"

setwd(outdir)

r <- readRDS("optimization_results.rds")

d <- r[sapply(r, class)=="data.frame"] %>%
      do.call("rbind", .) %>%
      group_by(taxonomy, lambda, algorithm, currency, era, slate) %>%
      mutate(rank = ifelse(algorithm=="backward", max(rank, na.rm=T)-rank+1, rank),
             priority = rank/max(rank)) %>%
      arrange(desc(rank)) %>%
      mutate(slope = lag(prop) - prop,
             slope = ifelse(is.na(slope), 0, slope),
             slope = scales::rescale(slope),
             marginal = -marginal)

# add data on land area
template <- stack("E:/phycon/data/cpad_cced/cpad_cced_raster_15km.tif")
land <- template[[2]] %>% rasterToPoints() %>% as.data.frame()
names(land)[3] <- "land_area"
d <- left_join(d, land)

land <- left_join(land, d %>% ungroup() %>% dplyr::select(x, y) %>% distinct() %>% mutate(data=T)) %>%
      mutate(data=ifelse(is.na(data), F, T),
             land_area=ifelse(data, land_area, NA))
template <- template[[2]]
template[!is.na(template)] <- land$land_area
writeRaster(template, "e:/phycon/shiny/cappa/raster_template.tif", overwrite=T)

d <- d %>%
      group_by(taxonomy, lambda, algorithm, currency, era, slate) %>%
      arrange(rank) %>%
      mutate(cum_area = cumsum(land_area))


lambda <- 1



for(taxonomy in c("chrono", "phylo", "clade", "species")){
      
      tax <- taxonomy
      
      pd <- d %>% 
            ungroup() %>%
            filter(algorithm=="forward", lambda==1, currency=="endemism", 
                   taxonomy==tax, slate=="existing") %>%
            arrange(id)
      cell <- pd$id[pd$rank==1]
      
      intactness <- raster("E:/phycon/data/intactness/intactness_15km.tif")
      
      con <- stack("E:/phycon/data/cpad_cced/cpad_cced_raster_15km.tif")
      names(con) <- c("con", "land")
      land <- con$land
      con <- con$con[]
      xy <- as.data.frame(coordinates(land))
      
      security <- function(x, lambda=1){
            lambda <- 2^lambda
            (1-(1-x)^lambda)^(1/lambda)
      }
      
      
      
      setwd(indir)
      if(taxonomy=="endemic species"){
            diversity_mat <- readRDS("site_by_species.rds")
            diversity_mat <- diversity_mat[,readRDS("site_by_species_endemic.rds")]
            branch_lengths <- rep(1, ncol(diversity_mat))
            tree <- starTree(colnames(diversity_mat))
      }
      if(taxonomy=="species"){
            diversity_mat <- readRDS("site_by_species.rds")
            branch_lengths <- rep(1, ncol(diversity_mat))
            tree <- starTree(colnames(diversity_mat))
      }
      if(taxonomy=="otu"){
            diversity_mat <- readRDS("site_by_otu.rds")
            branch_lengths <- rep(1, ncol(diversity_mat))
            tree <- starTree(colnames(diversity_mat))
      }
      if(taxonomy=="chrono"){
            diversity_mat <- readRDS("site_by_branch.rds")
            tree <- readRDS("site_by_chrono_phylogeny.rds")
            branch_lengths <- tree$edge.length / sum(tree$edge.length)
      }
      if(taxonomy=="phylo"){
            diversity_mat <- readRDS("site_by_branch.rds")
            tree <- readRDS("site_by_phylo_phylogeny.rds")
            branch_lengths <- tree$edge.length / sum(tree$edge.length)
      }
      if(taxonomy=="clade"){
            diversity_mat <- readRDS("site_by_branch.rds")
            tree <- readRDS("site_by_phylo_phylogeny.rds")
            tree$edge.length <- rep(1, length(tree$edge.length))
            branch_lengths <- rep(1, length(tree$edge.length))
      }
      
      #diversity_mat <- readRDS(paste0(indir, "/site_by_branch.rds"))
      #tree <- readRDS(paste0(indir, "/site_by_chrono_phylogeny.rds"))
      #branch_lengths <- tree$edge.length / sum(tree$edge.length)
      
      
      
      if(taxonomy %in% c("clade", "chrono", "phylo")){
            
            # label internal edges by distinctive OTU pair
            label_nodes <- function(tree){
                  
                  for(i in 1:tree$Nnode){
                        clade <- extract.clade(tree, length(tree$tip.label)+i)
                        root_node <- setdiff(1:nrow(clade$edge), clade$edge[,2])
                        if(length(clade$tip.label)==2) root_node <- 3
                        children <- clade$edge[,2][clade$edge[,1]==root_node] %>%
                              sapply(function(x){
                                    if(length(x)==0){
                                          clade$tip.label
                                    }else if(x<=length(clade$tip.label)){
                                          clade$tip.label[x]
                                    }else{
                                          sort(extract.clade(clade, x)$tip.label)[1]
                                    }
                              }) %>%
                              sort()
                        
                        label <- paste(paste(children, collapse=" & ")) %>%
                              paste0(" (", length(clade$tip.label), "-OTU clade)")
                        tree$node.label[i] <- label
                  }
                  
                  return(tree)
            }
            
            
            tree <- label_nodes(tree)
            tree$edge.label <- c(tree$tip.label, tree$node.label)[tree$edge[,2]]
            
            colnames(diversity_mat) <- tree$edge.label
      }
      
      
      ############## tree data ############
      
      td <- fortify(tree)
      
      # parents
      for(i in 1:nrow(td)){
            td$x0[i] <- td$x[td$node==td$parent[i]][1]
            td$y0[i] <- td$y[td$node==td$parent[i]][1]
      }
      
      # elbows
      td$x1 <- td$x0
      td$y1 <- td$y
      
      # radius & angle
      td$r <- td$x
      td$r0 <- td$x0
      td$r1 <- td$x1
      td$a <- td$y / max(td$y) * 2 * pi
      td$a0 <- td$y0 / max(td$y) * 2 * pi
      td$a1 <- td$y1 / max(td$y) * 2 * pi
      
      # polar coordinates
      td$xp <- td$r * cos(td$a)
      td$xp0 <- td$r0 * cos(td$a0)
      td$xp1 <- td$r1 * cos(td$a1)
      td$yp <- td$r * sin(td$a)
      td$yp0 <- td$r0 * sin(td$a0)
      td$yp1 <- td$r1 * sin(td$a1)
      
      td <- select(td, label, xp:yp1)
      
      td1 <- select(td, label, xp, yp, xp1, yp1)
      td2 <- select(td, label, xp0, yp0, xp1, yp1)
      names(td2) <- names(td1)
      td <- rbind(td1, td2)
      
      
      
      
      ##########  asdf  ###########
      
      inbounds <- apply(diversity_mat, 1, function(x) !all(is.na(x)))
      land <- land[inbounds]
      con <- con[inbounds]
      xy <- xy[inbounds,]
      intactness <- intactness[inbounds]
      
      # normalize occurrence probabilities by land area to correct edge effects
      diversity_mat <- apply(diversity_mat[inbounds,], 2, function(x) x*land)
      
      # normalize by range size to generate endemism
      endemism_mat <- apply(diversity_mat, 2, function(x) x/sum(x, na.rm=T))
      
      #ncol <- 100
      #ramp <- colorRampPalette(c("black", "blue", "red", "yellow"))(ncol)
      #plot(tree, edge.color=ramp[cut(cd$value, ncol)], show.tip.label=F, type="fan")
      #plot(rnorm(100), rnorm(100), col=ramp[cut(rnorm(100), ncol)])
      
      
      # Marginal conservation benefit
      R <- apply(diversity_mat, 2, sum, na.rm=T) # range sizes
      V <- branch_lengths
      C <- apply(diversity_mat, 2, function(p){
            x <- cbind(con, p) %>% na.omit()
            weighted.mean(x[,1], x[,2])})
      B <- security(C, lambda=lambda)
      MB <- apply(diversity_mat, 2, function(p) p*(1-con)) %>% # p gain
            apply(1, function(x) x/R) %>% t() %>% # gain as percentage of range size
            apply(1, function(x) C + x) %>% t() %>% # resulting C
            apply(c(1,2), security, lambda=lambda) %>% # resulting B
            apply(1, function(x) x-B) %>% t() %>%
            apply(1, function(x) x*V) %>% t()
      
      marginal <- apply(MB, 1, sum)
      
      
      ##### DIVERSITY STATISTICS #####
      
      # D: Diversity
      D <- apply(diversity_mat, 1, sum, na.rm=T)
      
      # PD: Phylogenetic diversity
      #V <- branch_lengths
      V[V==Inf] <- max(V[V!=Inf])
      PD <- apply(diversity_mat, 1, function(p) sum(p * V, na.rm=T))
      
      # E: Endemism, i.e. total endemic diversity, aka WE
      #R <- apply(diversity_mat, 2, sum, na.rm=T) # range sizes
      E <- apply(diversity_mat, 1, function(p) sum(p / R, na.rm=T))
      
      # PE: Phylogenetic endemism
      PE <- apply(diversity_mat, 1, function(p) sum(p * V / R, na.rm=T))
      
      # Em: Mean endemism
      # The following derivation is equivalent to ED / D
      Em <- apply(diversity_mat, 1, function(p) weighted.mean(1/R, w=p, na.rm=T))
      
      # PDm: Mean phylogenetic diversity, i.e. branch length of mean resident
      # The following derivation is equivalent to PD / D
      PDm <- apply(diversity_mat, 1, function(p) weighted.mean(V, w=p, na.rm=T))
      
      # PEm: Mean phylogenetic endemism, i.e. branch length / range size of mean resident
      # The following derivation is equivalent to PE / D
      PEm <- apply(diversity_mat, 1, function(p) weighted.mean(V / R, w=p, na.rm=T))
      
      # Mean branch length of the endemics
      BEm <- apply(diversity_mat, 1, function(p) weighted.mean(V, w=p/R, na.rm=T))
      
      pd <- pd %>%
            ungroup() %>%
            select(id, x, y, land_area, con, rank) %>%
            mutate(int=intactness, marginal=marginal,
                   D=D, PD=PD, E=E, PE=PE, Em=Em, PDm=PDm, PEm=PEm, BEm=BEm)
            
      data <- list(tree=tree, 
                   tree_data=td,
                   diversity_mat=diversity_mat, 
                   pd=pd, 
                   R=R, V=V, C=C, B=B, MB=MB)
      
      saveRDS(data, paste0("e:/phycon/shiny/cappa/data/facets/data_", taxonomy, ".rds"))
      #saveRDS(data, paste0("e:/phycon/shiny/sandbox/data_", taxonomy, ".rds"))
      
      # split this one to get under github size limit
      if(taxonomy == "species"){
            saveRDS(data[1:8], paste0("e:/phycon/shiny/cappa/data/facets/data_", taxonomy, ".rds"))
            saveRDS(data[9], paste0("e:/phycon/shiny/cappa/data/facets/data_", taxonomy, "2.rds"))
            
      }
}



