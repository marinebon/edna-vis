# setup ----

# load libraries
suppressPackageStartupMessages({
  library(shiny)
  library(tidyverse)
  library(lubridate)
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

# read in operational taxonomic units (otus)
otu = read_csv(otu_csv) %>%
  left_join(
    read_csv(sites_csv), 
    by=c('site'='site_code')) %>%
  mutate(
    date  = as_date(date))

n_otu_max = otu  %>%
  group_by(site) %>%
  summarize(
    n_otu = sum(count)) %>%
  .$n_otu %>% max()

# read in open tree of life
otl = read_csv(otl_csv)

# ui: user interface ----
ui <- fluidPage(
  navbarPage(
    title='eDNA explorer',
    # tabsetPanel(
    #   tabPanel(
    #     'Space & Time', icon = icon('map'),
    fluidRow(
      column(
        6,
        selectInput(
          'phylum', label = 'Filter by Phylum:', width='100%',
          unique(otu$phylum), multiple=T)),
      column(
        6,
        sliderInput(
          'date_range', label = 'Filter by Date:', width='100%',
          min = min(otu$date), max = max(otu$date), 
          value = c(min(otu$date), max(otu$date)),
          timeFormat='%Y-%m', animate=T))),
    fluidRow(
      column(
        6,
        leafletOutput('map')),
      column(
        6,
        plotlyOutput('plot'))),
    fluidRow(
      column(
        12,
        br(),
        DT::dataTableOutput('table')))))
# tabPanel(
#   'Taxa', icon = icon('sitemap')))))

# server: backend functions ----
server <- function(input, output) {
  
  d_f = reactive({
    
    if (!is.null(input$phylum)){
      otu_phylum = otu %>%
        filter(phylum %in% input$phylum)
    } else{
      otu_phylum = otu %>%
        mutate(
          phylum = 'ALL')
    }
    
    otu_phylum  %>%
      filter(
        date >= input$date_range[1],
        date <= input$date_range[2]) %>%
      group_by(phylum, site, site_name, lon, lat, date) %>%
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

    if (is.null(input$phylum) | length((input$phylum)) == 1){
      # color by site, plot of n_otu over time
      plot_ly(d_f(), x = ~date, y = ~n_otu, color = ~site, type='scatter', mode='lines+markers') %>%
        layout(
          legend = list(x = 0.1, y = 0.9),
          xaxis = list(title='', tickformat='%Y-%m'), 
          yaxis = list(title='OTUs'))
      
    } else {
      # color by phylum, summarize across sites
      plot_ly(
        d_f()  %>%
          group_by(phylum, date) %>%
          summarize(
            n_otu = sum(n_otu)) %>%
          ungroup(), 
        x = ~date, y = ~n_otu, color = ~phylum, type='scatter', mode='lines+markers') %>%
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
        select(Phylum = phylum, Site, Date=date, OTUs=n_otu),
      rownames=F, options = list(pageLength = 10), #, dom = 'tip', deferRender=TRUE, scrollY=300, scroller=TRUE),
      extensions="Scroller", style="bootstrap", class="compact", width="100%")
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

