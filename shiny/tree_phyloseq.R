
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
  
  otl_notfound = c(5264367, 632176, 621380, 67823, 955367, 588763, 566119, 3634672, 1083518, 2841628)
  tree <- tol_induced_subtree(
    ott_ids = setdiff(d_t$ott_id, otl_notfound), label_format = 'name')
  #Error: HTTP failure: 400 The following OTT ids were not found: 
  #  [5264367, 632176, 621380, 67823, 955367, 588763, 566119, 3634672, 1083518, 2841628]
  #Warning: In collapse_singles(tr) :
  #  Dropping singleton nodes with labels: Polysiphonia, Amphibalanus, Creseis, Abylopsis
  
# phyloseq -----

# tree_nwk  = '../data/tree.tre'
# tree <- tol_induced_subtree(setdiff(unique(d_t$ott_id), otl_notfound), label_format = 'name', file = tree_nwk)
# tree = ape::read.tree(tree_nwk)

tree$tip.label

otl = otl %>%
  mutate(
    ott_tip_label = str_replace(search_string, ' ', '_'))

otl %>% filter(DUP_ID == 'OTU_3626')


match(tree$tip.label, otl$ott_tip_label)

# ott_tip_label = str_replace(search_string, ' ', '_')) # "Jania"        "Cephalothrix" "Aurelia"
tree$tip.label[!tree$tip.label %in% otl$ott_tip_label] 
# ott_tip_label = str_replace(search_string, ' ', '_')) # "Lophopanopeus_lobipes" "Penaeus_aztecus"

View(tree$edge)
tree$tip.label %>% head
tree$node.label %>% head

# By convention, the tips of the tree are numbered 1 through n for n tips; and the nodes are numbered n + 1 through n + m for m nodes.

plot(tree, type = "cladogram") # edge.width = 2, label.offset = 0.1, 
nodelabels()
#tiplabels()

tr238 <- extract.clade(tree, 238)
plot(tr238)
nodelabels()

# phyloseq ----
library(phyloseq) # source('http://bioconductor.org/biocLite.R'); biocLite('phyloseq')

otu_txt   = '../data/OTU_table_taxa_all.txt'
otu_phyloseq_txt   = '../data/OTU_table_taxa_all_phyloseq.txt'

# setup columns for otu
otu_1 = read_tsv(otu_txt, n_max=1)
otu_cols = c('#OTU ID', names(otu_1)[-1], 'Consensus Lineage')

# read otu, in wide format with many 0s, and extra taxa vector column named
otu_w = read_tsv(otu_txt, col_names=otu_cols, skip=1) %>%
  mutate(
    `Consensus Lineage` = str_replace_all(`Consensus Lineage`, "([\\[\\]'])", "") %>% str_replace_all(",", ";"))
#s = "['k__Eukaryota','p__Arthropoda','c__Maxillopoda','o__Calanoida','f__Clausocalanidae','g__Clausocalanus','s__Clausocalanus furcatus']"
write_lines('# QIIME v1.2.1-dev OTU table', otu_phyloseq_txt)

# write otu file in format readable by phyloseq::import_qiime()
write_tsv(otu_w, otu_phyloseq_txt, append=T, col_names=T)

o = import_qiime(
  otufilename  = otu_phyloseq_txt) #, parseFunction = parse_taxonomy_greengenes) # parse_taxonomy_qiime) # parse_taxonomy_default)

# metaa
ntaxa(o)
nsamples(o)
sample_names(o)[1:5]
rank_names(o)
#sample_variables(o) # sam_data slot is empty.

otu_table(o)
tax_table(o) #tax_table(GlobalPatterns)

o_p = phyloseq(otu_table(o), tax_table(o))

o_s = subset_taxa(o_p, Family == "Mithracidae" & !is.na(Family))
plot_bar(o_s, fill = "Genus")

phy_tree(o_p)

taxa_sums(o_p)

otu_taxa = otu %>%
  group_by(kingdom, phylum, class, order, family, genus, species) %>%
  summarize(
    n_otu = n()) %>%
  arrange(kingdom, phylum, class, order, family, genus, species) %>%
  ungroup()

otu_taxa = otu_taxa %>%
  left_join(
    otu_taxa %>%
      group_by(kingdom, phylum, class, order, family, genus) %>%
      summarize(
        n_species_in_genus = n_distinct(species))) %>%
  left_join(
    otu_taxa %>%
      group_by(kingdom, phylum, class, order, family) %>%
      summarize(
        n_genus_in_family = n_distinct(genus))) %>%
  left_join(
    otu_taxa %>%
      group_by(kingdom, phylum, class, order) %>%
      summarize(
        n_family_in_order = n_distinct(family))) %>%
  left_join(
    otu_taxa %>%
      group_by(kingdom, phylum, class) %>%
      summarize(
        n_order_in_class = n_distinct(order))) %>%
  left_join(
    otu_taxa %>%
      group_by(kingdom, phylum) %>%
      summarize(
        n_class_in_phylum = n_distinct(class))) %>%
  left_join(
    otu_taxa %>%
      group_by(kingdom) %>%
      summarize(
        n_phylum_in_kingdom = n_distinct(phylum))) %>%
  select(
    kingdom, n_phylum_in_kingdom, phylum, n_class_in_phylum, class, n_order_in_class, order, n_family_in_order, family, n_genus_in_family, genus, n_species_in_genus, species, n_otu)

View(otu_taxa)  

data(GlobalPatterns)
GlobalPatterns
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 19216 taxa and 26 samples ]
## sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
## tax_table()   Taxonomy Table:    [ 19216 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 19216 tips and 19215 internal nodes ]
GP.chl = subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
#Map the sample environment ("SampleType") to point color, and taxonomic family to point shape. Additionally, label the tips with the Genus name and scale the point size by abundance
plot_tree(GP.chl, color = "SampleType", shape = "Family", label.tips = "Genus", 
          size = "abundance", plot.margin = 0.5, ladderize = TRUE)


# color.terminal.branches.R ----

source("color.terminal.branches.R")
color.terminal.branches(phy, states, breaks=4, cols=c("black","red"), edge.width=2, show.tip.label=TRUE)
color.terminal.branches(
  tree, set_names(d_t$n_otu, d_t$ott_name), type='cladogram', # "phylogram" (the default), "cladogram", "fan", "unrooted", "radial"
  breaks=10, cols=rev(brewer.pal(11,'Spectral')), edge.width=2, show.tip.label=TRUE)
