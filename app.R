

library(shiny)
library(shinydashboard)
library(rintrojs)
library(tidyverse)
library(raster)
library(rgdal)
library(DT)
library(ape)
library(phangorn)
library(viridis)
library(leaflet)
library(plotly)
library(scales)

select <- dplyr::select


metrics <- expand.grid(phylo=c(T, F),
                       diversity=c("community total", "community mean"),
                       endemism=c("ignore range size", "upweight small-range taxa"),
                       stringsAsFactors=F) %>%
      mutate(stat=c("PD", "D", "PDm", "Dm", "PE", "E", "PEm", "Em"),
             description=c("phylogenetic diversity (PD, the sum of the presence probabilities of the phylogenetic branch segments at a site, weighted by their branch lengths).", 
                           "diversity (D, the sum of the presence probabilities of the species at a site).", 
                           "mean phylogeneic diversity (PDm, the average of the branch lengths of the phylogenetic branch segments at a site, weighted by their presence probabilities).", 
                           "mean branch length (a meaningless uniform map, since species richness ignores branch length).",
                           "phylogenetic endemism (PE, the sum of the presence probabilities of the phylogenetic branch segments at a site, weighted by their branch lengths and by the inverse of their range sizes).", 
                           "endemism (E, the sum of the presence probabilities of the species at a site, weighted by the inverse of their range sizes).", 
                           "mean phylogeneic endemism (PEm, the average of the branch lengths divided by range sizes of the phylogenetic branch segments at a site, weighted by their presence probabilities).", 
                           "mean endemism (Em, the average of the inverse range sizes of the species at a site, weighted by their presence probabilities).")
      )

stats <- data.frame(var=c("biodiversity", "rank", "marginal", 
                          "con", "int",
                          "range"),
                    label=c("biodiversity", "conservation priority", "marginal value", 
                            "conservation status", "landscape intactness",
                            "taxon range"),
                    description=c("The spatial phylodiversity metric is defined by the two selections above and the branch length metric at left. Currently displaying",
                                  "Overall conservation priority for this diversity facet. 
                                  High priority sites contain concentrations of long-branch taxa with small ranges that are poorly protected by current preserves or other top-priority sites.",
                                  "Marginal value quantifies the amount of benefit that would derive from adding a site to the CURRENT set of protected lands, ignoring complementarity with other high-priority sites.", 
                                  "Current site conservation status, based on public land ownership and management, and private conservation easements. 
                                  1 = fully protected (e.g. national park), 0 = unprotected private land.",
                                  "Landscape intactness description", 
                                  "Distribution of the selected taxon within California, modeled based on herbarium specimen locations, climate, and landscape intactness."),
                    stringsAsFactors=F)



facets <- data.frame(var=c("time", "anagenesis", "cladogenesis", "species"),
                     label=c("survival time", "divergence", "diversification", "species richness"),
                     description=c("Measuring branch lengths in years highlights the total amount of time that independent lineages survived to give rise to a set of extant taxa.",
                                   "Measuring branch lengths in accumulated mutations highlights genetic divergence among taxa.", 
                                   "Measuring all brach segments as equal length emphasizes assemblages representing many speciation events.", 
                                   "Defining biodiversity as species richness ignores evolutionary relationships, treating each species as equally distinctive as represented by a 'star tree' with no internal branches and identical terminal branch lengths."),
                     file=c("data/facets/data_chrono.rds", "data/facets/data_phylo.rds",
                            "data/facets/data_clade.rds", "data/facets/data_species.rds"),
                     stringsAsFactors=F)

if(!file.exists("data/facets/data_species.rds")){
      c("data/facets/data_species1.rds", "data/facets/data_species2.rds") %>%
            lapply(readRDS) %>%
            do.call("c", .) %>%
            saveRDS("data/facets/data_species.rds")
}

lineages <- data.frame(label=c("branch length", "endemism", "range unprotection", 
                               "presence at selected site", "benefit at selected site", "highlight selected taxon"),
                       var=c("brl", "end", "opp", 
                             "div", "value", "selected"),
                       description=c("The length of a branch segment represents how much unique evolutionary heritage is shared by its descendant lineages.",
                                     "Endemism is the inverse of the combined California range size of a clade, a measure of extinction vulnerability.",
                                     "Range unprotection reflects the proportion of a taxon's California range that is not conserved",
                                     "The modeled probability that each taxon is found in the grid cell selected on the map.",
                                     "The relative contribution of each clade to the conservation benefit that would come from fully protecting the selected grid cell.",
                                     "The currently selected clade, which can be changed by clicking the tree or the table."),
                       stringsAsFactors=F)

palettes <- list(heat=c("white", "yellow", "#dd4b39"),
                 verdant=c("white", "darkolivegreen1", "darkgreen"),
                 nebula=c("black", "darkblue", "red", "yellow"),
                 grayscale=c("white", "black"))

# spatial data for grid cells 
r <- raster("data/raster_template.tif")
g <- rasterToPolygons(r) %>%
      spTransform(crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
g$id <- paste0("cell", 1:nrow(g))


#blurb <- readLines("blurb.txt")
jepson <- readRDS("data/jepson/jepson_dir.rds")

#button_style <- "white-space:normal; text-align:center; background-color:#fff; color:#000000; border-color:#fff; width:100%; height:150px;"
button_style <- "white-space:normal; text-align:center; background-color:#fff; color:#000000; border-color:#fff; width:100%;"


###### USER INTERFACE ######

ui <- dashboardPage(skin="blue",
                    
                    
                    dashboardHeader(disable=T),
                    dashboardSidebar(disable=T),
                    
                    dashboardBody(
                          
                          introjsUI(),
                          
                          tags$head(
                                # import css styles
                                tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
                                
                                # create an input$dimension variable representing window size
                                tags$script('
                                var dimension = [0, 0];
                                $(document).on("shiny:connected", function(e) {
                                    dimension[0] = window.innerWidth;
                                    dimension[1] = window.innerHeight;
                                    Shiny.onInputChange("dimension", dimension);
                                });
                                $(window).resize(function(e) {
                                    dimension[0] = window.innerWidth;
                                    dimension[1] = window.innerHeight;
                                    Shiny.onInputChange("dimension", dimension);
                                });
                            ')
                          ),
                          
                          fluidRow(
                                column(width=12,
                                       introBox(
                                             box(width=NULL, title=span("California Plant Phylodiversity Atlas", span("[beta]", style="color:#dd4b39; font-size:15px")),  
                                                 solidHeader=TRUE, status="danger", collapsible=T, collapsed=F,
                                                 fluidRow(
                                                       column(2,
                                                              actionButton("teaser", 
                                                                           span(h4(strong("Explore the evolutionary heritage and conservation landscape of a world floristic biodiversity hotspot."))), 
                                                                           style=button_style)),
                                                       column(2,
                                                              actionButton("tour", 
                                                                           span(h5(strong("New here?")), h5(strong("Take a tour."))), 
                                                                           icon=icon("globe", "fa-3x", lib="glyphicon"),
                                                                           style=button_style)),
                                                       column(2,
                                                              actionButton("citation", 
                                                                           span(h5("This site accompanies the paper",
                                                                                   "\"Kling M, Mishler B, Thornhill A, Baldwin B, Ackerly D. 2018 Facets of phylodiversity: evolutionary diversification, divergence and survival as conservation targets.",
                                                                                   strong(em("Phil. Trans. R. Soc. B\"")))), 
                                                                           onclick="window.open('http://dx.doi.org/10.1098/rstb.2017.0397', '_blank')",
                                                                           style=button_style)),
                                                       column(2,
                                                              actionButton("paper", 
                                                                           span(h5(strong("Want details?")), h5(strong("Read the paper."))), 
                                                                           icon=icon("file", "fa-3x", lib="glyphicon"),
                                                                           onclick="window.open('http://dx.doi.org/10.1098/rstb.2017.0397', '_blank')",
                                                                           style=button_style)),
                                                       column(2,
                                                              actionButton("credits", 
                                                                           span(h5("Part of the California Plant Phylodiversity Project, funded by the National Science Foundation.",
                                                                                   "App created by Matthew Kling. Site still under development so please pardon any buggy behavior.")), 
                                                                           style=button_style)),
                                                       column(2,
                                                              tags$button(
                                                                    id = "jepson",
                                                                    class = "btn action-button", 
                                                                    style=button_style,
                                                                    onclick="window.open('http://ucjeps.berkeley.edu/', '_blank')",
                                                                    tags$img(src='jepson.jpeg', width="75%", align="center", href="http://ucjeps.berkeley.edu/")
                                                              ))
                                                 )
                                             ),
                                             data.step=1,
                                             data.intro=span(h4(strong("Biodiversity under threat"), style="color:#dd4b39;"),
                                                             h5("Two criteria make California a designated global biodiversity hotspot: an exceptional number of native species, and widespread human degradation of natural ecosystems.",
                                                                "This interactive atlas allows you to explore patterns in the state's 5000+ native plant species from multiple angles, including their evolutionary relationships, geographic distributions, and community composition, as well as how protected each is across its range.",
                                                                "The tool accompanies a journal article focused on identifying future land conservation priorities---poorly-protected locations that are home to concentrations of distinctive and vulnerable native plants."))
                                       )
                                )
                          ),
                          
                          fluidRow(
                                column(width=12,
                                       introBox(
                                             box(width=NULL, 
                                                 title = "Taxon", #textOutput("selected_taxon"),
                                                 solidHeader=FALSE, status="danger", collapsible=F, collapsed=F,
                                                 #fluidRow(#column(1, actionButton("random_taxon", "", icon=icon("refresh", lib="glyphicon"),
                                                 #         #                       style="background-color:#fff; color:#000000; border-color:#fff; height:0px; align:center;")),
                                                 #         column(12, ))),
                                                 fluidRow(column(12, (uiOutput("taxon_description")))),
                                                 fluidRow(column(12, htmlOutput("photos")))   
                                             ),
                                             data.step=5,
                                             data.intro=span(h4(strong("Who are these plants anyway?"), style="color:#dd4b39;"),
                                                             h5("This panel shows images of the the taxon you've selected (or a random one to start).",
                                                                "Photos of species in the taxon are drawn randomly from the Jepson eFlora;",
                                                                "click images to view botanical and ecological information on the eFlora website."),
                                                             h5("Click the phylogeny or the community table below to select a specific taxon.", #, or click the refresh button to select a random taxon.",
                                                                "To map the range of the selected clade, change map variable to 'taxon range';",
                                                                "to view it in the evolutionary tree, change the branch color variable to 'highlight selected taxon'."))
                                       )
                                )
                          ),
                          
                          fluidRow(
                                
                                
                                column(width=4,
                                       introBox(
                                             box(width=NULL, title="Phylogeny", solidHeader=FALSE, status="danger",
                                                 
                                                 fluidRow(column(6, selectInput("facet", "Biodiversity facet", facets$label, selected="survival time")),
                                                          column(6, selectInput("phylocolor", "Branch color", lineages$label, selected="presence at selected site"))),
                                                 fluidRow(column(6, h5(textOutput("facet_caption"))),
                                                          column(6, h5(textOutput("phylocolor_caption")))),
                                                 fluidRow(column(12, plotlyOutput("tree", height=500))),
                                                 fluidRow(column(6, selectInput("palette_phylo", "Color palette", names(palettes), selected="verdant")),
                                                          column(6, selectInput("phytrans", "Color scale", c("linear", "log", "rank"))))
                                             ),
                                             data.step = 2,
                                             data.position="right",
                                             data.intro = span(h4(strong("What is biodiversity?"), style="color:#dd4b39;"),
                                                               h5("This tool compares four different facets of biodiversity, each based on a different concept of how related species are to one another.",
                                                                  "The traditional definition of simply counting species does not consider their relationships at all,",
                                                                  "while three alternative phylogenetic definitions each quantify diversity as the length of branches connecting species on the evolutionary tree.",
                                                                  "Using evolutionary distances recognizes that closely related species tend to share traits and histories,",
                                                                  "making species with few close relatives more distinctive and irreplaceable."),
                                                               h5("Pan, zoom, hover, and click to explore the tree and select taxa.",
                                                                  "Changes you make here will update the other figures as well."))
                                       )
                                       
                                ),
                                
                                
                                column(width=4, 
                                       introBox(
                                             box(width=NULL, title="Geography", solidHeader=FALSE, status="danger",
                                                 
                                                 fluidRow(column(12,
                                                                 selectInput("variable", "Map variable", 
                                                                             stats$label, selected="taxon range"))),
                                                 conditionalPanel(condition="input.variable == 'biodiversity'",
                                                                  fluidRow(
                                                                        column(6, selectInput("endemism", NULL, c("ignore range size", "upweight small-range taxa"), 
                                                                                              selected="upweight small-range taxa")),
                                                                        column(6, selectInput("diversity", NULL, c("community total", "community mean"), 
                                                                                              selected="community total")))),
                                                 fluidRow(column(12,
                                                                 h5(textOutput("stat_caption")),
                                                                 leafletOutput("mymap", height=600))),
                                                 fluidRow(column(4, selectInput("trans", "Color scale", c("linear", "log", "rank"))),
                                                          column(4, selectInput("palette", "Color palette", names(palettes), selected="verdant")),
                                                          column(4, sliderInput("opacity", "Color opacity", 0, 1, .5, step=.05))),
                                                 fluidRow(column(11, downloadButton("download_map", "download mapped data as CSV")))
                                             ),
                                             data.step = 3,
                                             data.position="right",
                                             data.intro = span(h4(strong("Where is biodiversity and conservation opportunity?"), style="color:#dd4b39;"),
                                                               h5("This map shows how biodiversity and conservation patterns vary across California,",
                                                                  "from individual species ranges to flora-wide biodiversity patterns to future conservation priorities."),
                                                               h5("Change the drop-down menu to select a variable to display.",
                                                                  "Click the map to change the selected site."))
                                       )
                                ),
                                
                                
                                column(width=4, 
                                       introBox(
                                             box(width=NULL, title="Community", solidHeader=FALSE, status="danger",
                                                 DT::dataTableOutput("table"),
                                                 h5("NOTE: Every taxon is listed for every cell, but only those with high 'presence' values influence the results.",
                                                    "All variables are rescaled as a fraction of the maximum value for any taxon."),
                                                 fluidRow(column(11, downloadButton("download_community", "download as CSV")))
                                             ),
                                             data.step = 4,
                                             data.position="left",
                                             data.intro = span(h4(strong("What is the local flora?"), style="color:#dd4b39;"), 
                                                               h5("Here you'll see a list of plants found in the cell selected on the map, along with details about how much total conservation 'benefit' each one contributes to the urgency of conserving the cell.",
                                                                  "This benefit comes from a combination of the likelihood the taxon is found in the cell (presence), how much unique evolutionary history it represents under the selected biodiversity facet (branch length)",
                                                                  "how small its range is in California (endemism), and how poorly protected it is across that range (range unprotection)."),
                                                               h5("Sort columns or search for taxa. Click on a taxon to select it.",
                                                                  "(If you want to search for a species, make sure 'species' is selected as the biodiversity facet at far left;",
                                                                  "otherwise the list of taxa represents every clade on the evolutionary tree, which captures relationships but is not fully resolved to the species level.)"))
                                       )
                                       
                                )
                                
                          )
                    )
)




##### SERVER-SIDE DATA PROCESSING #####

server <- function(input, output, session) {
      
      data <- reactive({ readRDS(facets$file[facets$label==input$facet]) })
      
      variable <- reactive({ metrics$stat[metrics$diversity==input$diversity &
                                                metrics$endemism==input$endemism &
                                                metrics$phylo==(input$facet != "species richness")] })
      stat_caption <- reactive({
            caption <- stats$description[stats$label==input$variable]
            if(input$variable=="biodiversity") caption <- paste(caption, paste0(input$facet, "-based"), metrics$description[metrics$stat==variable()])
            if(input$variable=="taxon range" & !is.null(selected_taxon())) caption <- sub("the selected taxon", selected_taxon(), caption)
            if(input$variable=="taxon range" & is.null(selected_taxon())) caption <- "Click a taxon in the phylogeny or the table to view its geographic range within California."
            caption
      })
      output$stat_caption <- renderText(stat_caption())
      
      
      
      #### define selected taxon based on inputs from alternative sources
      selected_taxon <- reactiveVal()
      
      observeEvent(is.null(selected_taxon()), {if(is.null(selected_taxon())) selected_taxon(sample(colnames(data()$diversity_mat), 1))}, priority=105)
      observeEvent(input$random_taxon, {selected_taxon(sample(colnames(data()$diversity_mat), 1))}, priority=105)
      observeEvent(input$table_rows_selected, {selected_taxon(as.character(comp()$name[input$table_rows_selected]))}, priority=60)
      observeEvent(event_data("plotly_click"), {selected_taxon(data()$tree_data$label[as.integer(event_data("plotly_click")$key)])}, priority=50)
      
      #observeEvent(input$facet, {data()$pd[1,]}, priority=100)
      observeEvent(input$facet, {if(! selected_taxon() %in% colnames(data()$diversity_mat)) selected_taxon(sample(colnames(data()$diversity_mat), 1)) }, priority=95)
      #observeEvent(input$facet, {selected_taxon(sample(colnames(data()$diversity_mat), 1))})
      
      
      output$selected_taxon <- renderText({ selected_taxon() })
      
      
      ########## MAP PANEL COMPONENTS #######
      
      # values and color ramp for map grid cells
      map_data <- reactive({
            var <- stats$var[stats$label==input$variable]
            if(var=="biodiversity") var <- variable()
            
            if (#!is.null(selected_taxon()) & 
                  var=="range"){
                  map_values <- data()$diversity_mat[,selected_taxon()] / data()$pd$land_area
            } else if (is.null(selected_taxon()) & var=="range"){
                  map_values <- rep(0, nrow(data()$diversity_mat))
            } else if (var=="Dm"){
                  map_values <- rep(1, nrow(data()$diversity_mat))
            } else {
                  map_values <- unlist(data()$pd[,var])
            }
            
            
            if(input$trans=="log") map_values <- log10(map_values)
            if(input$trans=="rank") map_values <- rank(map_values)
            
            if(input$variable=="conservation priority") map_values <- max(map_values) - map_values
            
            grid <- g
            grid$values <- map_values
            pal <- colorNumeric(palettes[[input$palette]], domain=map_values)
            #pal <- colorNumeric(palettes[["heat"]], domain=map_values)
            
            list(grid=grid, pal=pal)
      })
      
      
      # basemap
      output$mymap <- renderLeaflet({
            leaflet(g) %>%
                  setView(-120, 37, zoom = 6) %>%
                  #fitBounds(-124.5, 32, -114, 42) %>%
                  addMapPane("background_map", zIndex = 1) %>% 
                  addMapPane("polygons", zIndex = 2) %>% 
                  addMapPane("labels", zIndex = 3) %>%  
                  addMapPane("top", zIndex = 4) %>% 
                  addProviderTiles(providers$Esri.WorldTerrain,
                                   options = pathOptions(pane = "background_map")) %>%
                  addProviderTiles(providers$Stamen.TonerLabels,
                                   options = pathOptions(pane = "labels", opacity=.65)) %>%
                  addProviderTiles(providers$Stamen.TonerLines,
                                   options = pathOptions(pane = "labels", opacity=.5))
      })
      
      # identify user-clicked grid cell
      selected_cell <- reactive({
            clicked <- input$mymap_shape_click
            if(is.null(clicked)){
                  point <- 119 # highest-priority cell for chronogram diversity
            }else{
                  coords <- as.data.frame(cbind(clicked$lng, clicked$lat))
                  point <- SpatialPoints(coords)
                  proj4string(point) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
            }
            cell <- g[point,]
            return(list(shape=cell, id=cell$id))
      })
      
      # grid cells layer -- updates with changes to map data
      observe({
            pal <- map_data()$pal
            leafletProxy("mymap", data=map_data()$grid) %>%
                  addPolygons(layerId=map_data()$grid$id,
                              color="gray90", weight=.1,
                              fillColor= ~pal(values), fillOpacity=input$opacity,
                              highlightOptions = highlightOptions(color = "blue", weight = 6, bringToFront = TRUE),
                              options=pathOptions(pane="polygons"))
      })
      
      # selected cell -- updates with map click
      observe({
            leafletProxy("mymap") %>%
                  addPolygons(data=selected_cell()$shape,
                              layerId="selected",
                              color="cyan", weight=4, opacity=1, fillColor="transparent",
                              options=pathOptions(pane="top"))
      })
      
      
      output$download_map <- downloadHandler(
            filename = function() {
                  "map_data.csv"
            },
            content = function(con) {
                  map_data()$grid %>%
                        as.data.frame() %>%
                        cbind(data()$pd[,c("x", "y")]) %>%
                        dplyr::select(x, y, values) %>%
                        rename(cell_value=values) %>%
                        write.csv(con, row.names=F)
            }
      )
      
      
      
      
      
      ############ COMMUNITY PANEL COMPONENTS ##############
      
      cda <- reactive({
            scl <- function(x) x / sum(x)
            i <- selected_cell()$id %>% substr(5, 15) %>% as.numeric()
            data.frame(
                  name = colnames(data()$diversity_mat),
                  div = data()$diversity_mat[i,],
                  end = 1/data()$R,
                  brl = data()$V,
                  opp = 1-data()$B,
                  value = data()$MB[i,]) %>%
                  mutate_at(vars(div:value), funs(scl))
      })
      
      comp <- reactive({
            scl <- function(x) x / max(x, na.rm=T)
            cda() %>%
                  arrange(desc(value)) %>%
                  mutate_at(vars(div:value), funs(scl)) %>%
                  mutate(cum_value = cumsum(value/sum(value)),
                         rank = 1:nrow(.)) %>%
                  gather(var, prop, div:cum_value) %>%
                  mutate(var = factor(var, levels=c("div", "end", "brl", "opp", "value", "cum_value"),
                                      labels=c("presence", "endemism", "branch length", "range unprotection",
                                               "combined benefit", "cumulative benefit"))) %>%
                  filter(var %in% c("combined benefit", "presence", "endemism", "branch length", "range unprotection")) %>%
                  mutate(prop = round(prop, 4)) %>%
                  spread(var, prop) %>%
                  arrange(rank) %>%
                  dplyr::select("name", "combined benefit", "presence", "endemism", "branch length", "range unprotection")
      })
      
      output$table <- renderDataTable(comp() %>% DT::datatable(selection=list(mode="single",
                                                                              #selected=which(comp()$name==selected_taxon()), 
                                                                              target="row"),
                                                               options=list(scrollX=T, scrollY=T)))
      
      output$download_community <- downloadHandler(
            filename = function() {
                  "community_data.csv"
            },
            content = function(con) {
                  comp() %>%
                        write.csv(con, row.names=F)
            }
      )
      
      
      
      
      ########## PHYLOGENY PANEL COMPONENTS ########
      
      phylo_data <- reactive({
            
            pd <- data()$tree_data
            
            if(input$phylocolor=="highlight selected taxon"){
                  v <- dplyr::select(cda(), name)
                  
                  if(grepl("OTU", selected_taxon())){ # selected taxon is a node
                        kids <- extract.clade(data()$tree, selected_taxon())
                        v$value <- as.integer(v$name %in% c(kids$node.label, kids$tip.label))
                  }else{ # selected taxon is a tip
                        v$value <- 0
                        v$value[v$name==selected_taxon()] <- 1
                  }
                  names(v) <- c("label", "value")
            }else{
                  v <- cda()[,c("name", lineages$var[lineages$label==input$phylocolor])]
                  names(v) <- c("label", "value")
                  if(input$phytrans=="log") v$value <- log10(v$value)
                  if(input$phytrans=="rank") v$value <- rank(v$value)
                  v$value[!is.finite(v$value)] <- min(v$value[is.finite(v$value)])
            }
            
            pd <- left_join(pd, v)
            
            clrs <- sub("white", "gray90", palettes[[input$palette_phylo]])
            pd$color <- colorRampPalette(clrs)(10)[cut(pd$value, 10)]
            
            
            pd
      })
      
      output$tree <- renderPlotly({
            pd <- phylo_data()
            key <- row.names(pd)
            
            pd <- pd %>% 
                  #arrange(value) %>% # sorting to plot high values on top, but breaks taxon selection
                  select(yp, yp1, xp, xp1, color, label)
            
            p <- ggplot(pd, aes(x=yp1, y=xp1, xend=yp, yend=xp, key=key, text=label)) + 
                  geom_segment(color=pd$color) +
                  coord_fixed() +
                  theme_void() +
                  theme(panel.grid.major=element_blank(),
                        legend.position="none")
            ggplotly(p, tooltip=c("text"))
      })
      
      output$tree_height <- renderText({ round(input$dimension[1]/3) })
      
      output$facet_caption <- renderText(facets$description[facets$label==input$facet])
      
      output$phylocolor_caption <- renderText(lineages$description[lineages$label==input$phylocolor])
      
      
      
      
      ####### TAXON PANEL #######
      
      jep <- reactive({
            
            width <- input$dimension[1]
            height <- round(width / 12)
            
            max_pics <- 15
            
            sel_tax <- selected_taxon()
            #if(is.null(selected_taxon())) sel_tax <- sample(colnames(data()$diversity_mat), 1)
            
            imgs <- function(urls){
                  url <- c()
                  for(u in 1:nrow(urls)) url <- c(url, '<a href="', urls$species_url[u], '" target="_blank"><img src="', 
                                                  urls$pic_urls[u], '" title="', urls$species[u],
                                                  '" height="', height, 'px" hspace="5px">')
                  
                  url <- c('<div class="photos">', url, '</div>')
                  return(url)
            }
            
            if(input$facet=="species richness"){ # species
                  st <- sub("_", " ", sel_tax)
                  j <- filter(jepson, species==st)
                  
                  url <- select(j, species, species_url, pic_urls) %>%
                        unique() %>%
                        sample_n(min(nrow(.), max_pics))
                  
                  return(list(data=j,
                              level="species",
                              text=paste(st, "is a native California plant species."),
                              text2=NULL,
                              img=imgs(url)))
            }
            
            st <- str_trim(substr(sel_tax, 1, regexpr(" \\(", sel_tax)))
            if(nchar(st)==0) st <- sel_tax
            
            if(!grepl("OTU", sel_tax)){ # terminals
                  j <- filter(jepson, clade==st)
                  
                  url <- select(j, species, species_url, pic_urls) %>%
                        unique() %>%
                        sample_n(min(nrow(.), max_pics))
                  
                  txt <- paste(st, "is a terminal clade (OTU) containing", 
                               length(unique(j$species)), "California species.")
                  if(length(unique(j$species))==1) txt <- sub("species",
                                                              paste0("species (", unique(j$species), ")"), txt)
                  
                  list(data=j,
                       level="terminal",
                       text=txt,
                       text2=NULL,
                       img=imgs(url))
                  
            }else{ # internal branches
                  tips <- extract.clade(data()$tree, sel_tax)$tip.label
                  j <- filter(jepson, clade %in% tips)
                  
                  url <- select(j, species, species_url, pic_urls) %>%
                        unique() %>%
                        sample_n(min(nrow(.), max_pics))
                  
                  list(data=j,
                       level="internal clade",
                       text=paste(st, "is a clade containing", 
                                  length(unique(j$clade)), "terminal OTUs and",
                                  length(unique(j$species)), "California species."),
                       text2="Clades are named such that they represent the most recent common ancestor of the two named tips.",
                       img=imgs(url))
            }
      })
      
      output$taxon_description <- renderUI({
            tagList(strong(jep()$text),
                    jep()$text2,
                    "Photos are from", 
                    a("CalPhotos", href="https://calphotos.berkeley.edu/"), 
                    "and have clickable links to the corresponding description in the", 
                    a("Jepson eFlora.", href="http://ucjeps.berkeley.edu/eflora/"))
      })
      
      
      output$photos <- renderText({ jep()$img })
      
      
      
      
      ####### TOUR #######
      
      observeEvent(input$tour,
                   introjs(session, options = list("nextLabel"="next",
                                                   "prevLabel"="back",
                                                   "skipLabel"="done"))
      )
      
}

# Run the application 
shinyApp(ui = ui, server = server)

