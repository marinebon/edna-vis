# edna-vis
visualization app for environmental DNA (eDNA), using R shiny

# future: blogdown site

[rstudio/blogdown: Create Blogs and Websites with R Markdown](https://github.com/rstudio/blogdown)

## Initialize

```r
# install
devtools::install_github('rstudio/blogdown')
blogdown::install_hugo()

# create a new site
blogdown::new_site()

# install theme
blogdown::install_theme('jrutheiser/hugo-lithium-theme')
blogdown::install_theme('christianmendoza/hugo-smpl-theme')
blogdown::install_theme('kakawait/hugo-tranquilpeak-theme')
blogdown::install_theme('jpescador/hugo-future-imperfect')
```

## Edit

```
# run server in background
options(servr.daemon = TRUE)
blogdown::serve_site()

#To stop the server, restart your R session or
servr::daemon_stop("4333969472")
```