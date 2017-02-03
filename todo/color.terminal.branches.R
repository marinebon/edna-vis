

# phy: object of class phylo
# data: named numeric vector
# breaks: number of color breaks
# cols: character vector of length 2, start col, end col
# non.terminal.col: vector character length 1, color non terminal branches

color.terminal.branches <- function(phy, data, breaks=6, cols=c("black", "red"), non.terminal.col= "black", edge.width=1, col.bias=1, legend.title = "", show.tip.label=FALSE, alt.col.data=NULL, ...){

sapply(phy$tip.label, check.spp.search.duplicates, phy$tip.label)

# Number of tips in the phylogeny
n.tips <- length(phy$tip.label)

# This will be the main vector which colors will be assigned to
col.vector <- vector(mode="character",length=nrow(phy$edge))

# Non-terminal branches will be black
col.vector[phy$edge[,2]>n.tips] <- non.terminal.col


# Create a vector of colors 
colors <-  colorRampPalette(cols, bias=col.bias)(breaks)

# range of trait values present
if(is.null(alt.col.data)==TRUE){
range.vals <- range(data)
} else {
range.vals <- range(c(alt.col.data, data))
}
diff <- range.vals[2] - range.vals[1]

bin.size <- diff / breaks

bin.ref <- as.data.frame(matrix(NA, ncol=3, nrow=breaks))
colnames(bin.ref) <- c("Bin","Start","End")
bin.ref$Bin <- seq(1, breaks, 1)
bin.ref[1,2] <- range.vals[1]
bin.ref[1,3] <- range.vals[1] + bin.size

for(i in 2:nrow(bin.ref)){

bin.ref[i,2] <- bin.ref[i-1,3]+0.0001
bin.ref[i,3] <- bin.ref[i,2] + bin.size

}


## for legend
legend.text <- vector()

for(i in 1:nrow(bin.ref)){

t <- paste(as.character(bin.ref[i,2]), " - ",as.character(bin.ref[i,3]))

legend.text[i] <- t

}


# create a duplicate of phy$edge which can be manipulated as required
edge.data <- as.data.frame(phy$edge)


# Yes, its a loop, but its easier
for(i in 1:length(phy$tip.label)){

spp <- phy$tip.label[i]

spp.trait.val <- as.numeric(data[grep(spp,names(data))])

bin.num <- bin.ref[spp.trait.val >= bin.ref$Start & spp.trait.val <= bin.ref$End,1]

spp.col <- colors[bin.num]

edge.row <- as.numeric(rownames(edge.data[edge.data$V2==i,]))

col.vector[edge.row] <- spp.col

}

max.val <- range.vals[2]

#layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
l <-  matrix(c(1,1,2),ncol=3,nrow=1)
layout(l)
#par(mfrow=c(1,2))
par(mai=rep(0.1,4))
plot(phy, show.tip.label=show.tip.label, edge.col=col.vector, edge.width=edge.width, ...)

#image <- matrix(colors, ncol=1)
#image <- as.matrix(rev(legend_image))
#legend_image <- as.raster(image)
legend_image <- as.raster(matrix(rev(colors), ncol=1))

plot(c(0,.9),c(0,max.val),type = 'n', axes = F,xlab = '', ylab = '', main = legend.title, xlim=c(0,.2))

rasterImage(legend_image, 0, 0, .08,max.val)
text(x=0.11, y = seq(0,max.val,l=5), labels = seq(0,max.val,l=5))
}

###

check.spp.search.duplicates <- function(search.spp, all.spp){

returns <- grep(search.spp, all.spp)

if(length(returns)>1){

stop("The species  ",search.spp,"  will cause problems - searching for this species in all the species present will return multiple species, consider renaming!")

}

}


