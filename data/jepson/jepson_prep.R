
library(tidyverse)

source("e:/chilefornia/chilefornia/jepson_parse.r")

taxa <- readRDS("e:/phycon/shiny/cappa/facets/data_species.rds")$diversity_mat %>%
      colnames() %>%
      data.frame(gs=.) %>%
      separate(gs, c("genus", "species"), sep="_", remove=F)

genera <- unique(taxa$genus)
urls <- sapply(genera, genus_urls)
gurls <- sapply(urls, function(x) ifelse(length(x)>0, x, NA))
g <- data.frame(genus=genera, genus_url=gurls, stringsAsFactors=F)
 
surls <- lapply(g$genus_url, function(x) if(is.na(x)) return(NA) else(return(try(species_urls(x)))))

s <- data.frame()
for(i in 1:length(surls)){
      surl <- surls[[i]]
      if(class(surl) != "character") surl <- NA
      
      sd <- data.frame(genus=g$genus[i],
                       genus_url=g$genus_url[i],
                       species_url=surl,
                       stringsAsFactors=F)
      s <- rbind(s, sd)
}


species_info <- function(species_url){
      #species_url <- "http://ucjeps.berkeley.edu/eflora/eflora_display.php?tid=49279"
      #species_url <- "http://ucjeps.berkeley.edu/eflora/eflora_display.php?tid=71178"
      if(is.na(species_url)) return(data.frame())
      
      html <- read_html(species_url)
      
      vars <- html %>%
            html_nodes("table+ table div") %>%
            html_children()
      
      taxon <- vars[1] %>%
            html_children() %>%
            html_text() %>%
            paste(collapse=" var. ")
      
      
      pics <- read_html(species_url) %>%
            html_nodes("img") %>%
            as.character() %>%
            str_split(" ") %>%
            sapply(function(x) x[grepl("src=", x)]) %>%
            gsub("\"|src=", "", .)
      pics <- pics[grepl("calphotos|illustrations", pics)]
      
      if(length(pics)==0) pics <- NA
      
      data.frame(species=taxon, species_url=species_url, pic_urls=pics)
}


d <- lapply(s$species_url, function(x) try(species_info(x)))


ds <- d[sapply(d, class)=="data.frame"] %>%
      Reduce("rbind", .) %>%
      full_join(s)

ds <- ds %>%
      mutate(species=as.character(species),
             species_url=as.character(species_url),
             pic_urls=as.character(pic_urls)) %>%
      separate(species, c("species", "subspecies"), sep=" var. ")



# add data on clades
spcl <- data.table::fread("E:/phycon/data/occurrences/California_Species_clean_All_epsg_3310.csv") %>%
      dplyr::select(clade, current_name_binomial) %>%
      mutate(clade = ifelse(clade=="Solanum", "Solanum_Californian", clade)) %>% # corrects an error in the underlying data
      distinct() %>% 
      mutate(clade=gsub(" ", "_", clade)) %>%
      rename(species=current_name_binomial)

p <- full_join(spcl, ds)

#p <- select(p, species, species_url, pic_urls) %>% distinct()

saveRDS(p, "e:/phycon/shiny/cappa/jepson/jepson_dir.rds")
