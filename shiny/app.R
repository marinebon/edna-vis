# setup ----

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(stringr)
  library(shiny)
  library(shinydashboard)
  library(shinythemes)
  library(leaflet)
  library(scales)
  library(plotly)
  library(RColorBrewer)
  library(DT)
  select = dplyr::select
})


# load data
# setwd('shiny') # for debug, when not running app

if (!file.exists('data')){
  file.copy('../data', '.', recursive=T)
}

otu_csv   = 'data/otu.csv'
otl_csv   = 'data/otl.csv'   # otl_csv   = 'data/otl.csv'
sites_csv = 'data/sites.csv' # sites_csv = 'data/sites.csv'

# read in csv's
otu   = read_csv(otu_csv)
otl   = read_csv(otl_csv)
sites = read_csv(sites_csv)

# merge data
otu = otu %>%
  left_join(
    sites, 
    by=c('site'='site_code')) %>%
  left_join(
    otl, 
    by='DUP_ID')

n_otu_max = otu  %>%
  group_by(site) %>%
  summarize(
    n_otu = sum(count)) %>%
  .$n_otu %>% max()

# ui: user interface ----
ui <- dashboardPage(
  skin = 'blue', # theme = shinytheme("slate"), # themeSelector(),
                      
  dashboardHeader(
    titleWidth=250,
    title=span(tagList(icon('tint'), 'eDNA Explorer'))),
    # dropdownMenu(
    #   type = "notifications",
    #   messageItem(
    #     from = "Help",
    #     message = "Understand the user interface.",
    #     icon = icon("question")))),
  
  dashboardSidebar(
    width=250,
    
    # TODO: add link to Help with icon('question')
    
    sidebarMenu(
      menuItem(
        "Charts", tabName = "charts" , icon = icon("line-chart")),
      menuItem(
        "Tree"  , tabName = "tree"   , icon = icon("sitemap"))),
    
    hr(),
    strong('BIOM File Input:'),
    selectInput(
      'sel_file', label = 'Select existing:', width='100%',
      'COI_0316_json_obs_md.biom', multiple=T),
    
    fileInput('up_file', 'Upload your own:',
              accept=c('text/csv', 
                       'text/comma-separated-values,text/plain', 
                       '.csv')),
    hr(),
    strong('Filters:'),
    selectInput(
      'rank', label = 'Taxa, Rank:', width='100%',
      # paste(sprintf("'%s'='%s'", stringr::str_to_title(names(otu)), names(otu)), collapse=',')
      c('Phylum'='phylum','Class'='class','Order'='order','Family'='family','Genus'='genus','Species'='species'), 
      multiple=F), 
    selectInput(
      'taxa', label = 'Taxa, Values:', width='100%',
      unique(otu[['phylum']]), multiple=T),
    selectInput(
      'sites', label = 'Sites:', width='100%',
      with(sites, set_names(site_code, sprintf('%s: %s', site_code, site_name))),
      multiple=T),
    sliderInput(
      'date_range', label = 'Date:', width='100%',
      min = min(otu$date), max = max(otu$date), 
      value = c(min(otu$date), max(otu$date)),
      timeFormat='%Y-%m', animate=T)),
  
  dashboardBody(
    tags$head(tags$link(rel='stylesheet', type ='text/css', href='styles.css')),
    tabItems(
      tabItem(
        tabName = "charts",
        box(
          leafletOutput('map')),
        box(
          plotlyOutput('plot')),
        box(
          width=12, 
          DT::dataTableOutput('table'))),
      tabItem(
        tabName = "tree",
        box(
          width=12, height=665,
          plotOutput('tree'))))))

# server: backend functions ----
server <- function(input, output, session) {
  
#   showModal(modalDialog(
#     title = "Overview",
#     HTML(markdownToHTML(
#       text=
#         "
# - **Filter** on the left 
# - **Map** summarizes across time
# - **Plot** summarizes across space when more than one taxa chosen")))) #, easyClose=T, size='l'

  # taxa, update with rank ----
  observe({

    # Can use character(0) to remove all choices
    if (is.null(input$rank)){
      x <- character(0)
    } else {
      x <- unique(otu[[input$rank]])
    }
    
    # Can also set the label and select items
    updateSelectInput(
      session, "taxa",
      choices = x)
  })
  
  # data, filtered ----
  d_f = reactive({
    
    # filter by taxa
    if (!is.null(input$taxa)){
      idx = otu[[input$rank]] %in% input$taxa
      cat(file=stderr(), str(sum(idx)))
      otu_f      = otu[idx,]
      otu_f$rank = input$rank
      otu_f$taxa = otu_f[[input$rank]]
    } else{
      otu_f = otu %>%
        mutate(
          rank = 'ANY',
          taxa = 'ALL')
    }
    
    # filter by sites
    if (!is.null(input$sites)){
      cat(file=stderr(), str(input$sites))
      otu_f = filter(otu_f, site %in% input$sites)
    }
    
    otu_f  %>%
      # filter by date
      filter(
        date >= input$date_range[1],
        date <= input$date_range[2]) %>%
      # summarize
      group_by(rank, taxa, site, site_name, lon, lat, date) %>%
      summarize(
        n_otu = sum(count)) %>%
      ungroup() 
    })
  
  # map ----
  output$map <- renderLeaflet({
    
    d_m = d_f()  %>%
      group_by(site, site_name, lon, lat) %>%
      summarize(
        n_otu = sum(n_otu)) %>%
      ungroup()
    
    # color palette
    if (length(unique(d_m$site)) <= 1){
      col_vals = c(0, n_otu_max)
    } else {
      col_vals = d_m$n_otu
    }
    pal <- colorNumeric('YlOrRd', col_vals)
    
    # map of sites colored by OTUs
    leaflet(
      d_m, 
      options = leafletOptions(minZoom = 4, maxZoom = 10)) %>% 
      addProviderTiles('Esri.OceanBasemap') %>%
      addCircleMarkers(
        ~lon, ~lat,
        radius = ~rescale(n_otu, to=c(10, 50)),
        color = ~pal(n_otu),
        stroke = FALSE, fillOpacity = 0.5,
        popup = ~sprintf("%s (%s): %s", site_name, site, format(n_otu, big.mark=",", scientific=F)), 
        label = ~site_name) %>%
      addLegend(
        "topleft", pal = pal, values = col_vals,
        title = "OTUs", opacity = 1)
  })
  
  # plot ----
  output$plot <- renderPlotly({

    if (is.null(input$taxa) | length((input$taxa)) == 1){
    
      #saveRDS(d_f(), 'test_d_f4plot.rds')
      #setwd('shiny')
      #X = readRDS('test_d_f4plot.rds')
      
      # color by site, plot of n_otu over time
      plot_ly(d_f(), x = ~date, y = ~n_otu, color = ~site, type='scatter', mode='lines+markers') %>%
      #plot_ly(X, x = ~date, y = ~n_otu, color = ~site, type='scatter', mode='lines+markers') %>%
        layout(
          legend = list(x = 0.01, y = 0.99),
          xaxis = list(title='', tickformat='%Y-%m'), 
          yaxis = list(title='OTUs'))
      
    } else {
      # color by taxa, summarize across sites
      plot_ly(
        d_f()  %>%
          group_by(taxa, date) %>%
          summarize(
            n_otu = sum(n_otu)) %>%
          ungroup(), 
        x = ~date, y = ~n_otu, color = ~taxa, type='scatter', mode='lines+markers') %>%
        layout(
          legend = list(x = 0.1, y = 0.9),
          xaxis = list(title='', tickformat='%Y-%m'), 
          yaxis = list(title='OTUs'))
      
    }
  })
   
  # table ----
  output$table <- DT::renderDataTable({
    datatable(
      d_f() %>%
        mutate(
          Site = sprintf('%s: %s', site, site_name)) %>% 
        select(-site, -site_name, -lon, -lat) %>% 
        select(Rank = rank, Taxa = taxa, Site, Date=date, OTUs=n_otu),
      rownames=F, options = list(pageLength = 5), #, dom = 'tip', deferRender=TRUE, scrollY=300, scroller=TRUE),
      extensions="Scroller", style="bootstrap", class="compact", width="100%")
  })
  
  
  # tree ----
  
  output$tree <- renderPlot({
    
    library(ape)
    library(rotl)
    library(ggtree) # source("https://bioconductor.org/biocLite.R"); biocLite("BiocUpgrade"); biocLite("ggtree", type = "source")
    library(phylobase)
    
    # summarize data for taxonomic tree
    d_t = otu %>%
      mutate(
        ott_name = str_replace(unique_name, ' ', '_')) %>%
      group_by(ott_id, ott_name) %>%
      summarize(
        n_otu = sum(count))
    
    # color
    pal = col_numeric('Spectral', d_t$n_otu)
    d_t = d_t %>%
      # TODO: left_join(d_f()), but need unique_name in otu -> d_f() to get filtered n_otu
      mutate(
        col_otu = pal(n_otu))
    
    ott_notfound = c(5264367, 632176, 621380, 67823, 955367, 588763, 566119, 3634672, 1083518, 2841628)
    tree <- tol_induced_subtree(
      ott_ids = setdiff(d_t$ott_id, ott_notfound), label_format = 'name')
    #Error: HTTP failure: 400 The following OTT ids were not found: 
    #  [5264367, 632176, 621380, 67823, 955367, 588763, 566119, 3634672, 1083518, 2841628]
    #Warning: In collapse_singles(tr) :
    #  Dropping singleton nodes with labels: Polysiphonia, Amphibalanus, Creseis, Abylopsis
    
    # ggtree ----
    
    # join data using phylo4
    d_t4 = d_t %>% 
      bind_rows(data_frame(
        ott_name = tree$tip.label[!tree$tip.label %in% d_t$ott_name],
        col_otu  = 'gray')) %>%
      as.data.frame()
    rownames(d_t4) = d_t4$ott_name
    
    tree4 = phylo4d(as(tree, 'phylo4'), d_t4[tree$tip.label,])
    
    #ggtree(tree4) %<+% d_t4 + aes(color=I(col_otu))
    ggtree(tree4, layout='radial') + # 'rectangular', 'slanted', 'fan', 'circular', 'radial' or 'unrooted'
      #scale_color_continuous(low='darkgreen', high='red') +
      #theme(legend.position="right") +
      geom_tippoint(aes(color=I(col_otu))) +
      geom_tiplab(size=3, aes(angle=angle))
    
    
  }, width=650, height=650)
}

# Run the application 
shinyApp(ui = ui, server = server)

