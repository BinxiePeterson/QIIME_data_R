######################################################
# Loading tables in R
######################################################  

### OTU tables  

# Note: This the the rarefied one!
otu <- read.table("otutable_r10k.txt", 
                  comment="", 
                  header=TRUE, 
                  sep="\t",
                  skip=1, 
                  as.is=TRUE, 
                  check.names=F,
                  row=1)

# View first 2 lines using head()
head(otu, n=2)

# View dimensions
dim(otu)

# Print row names (which are OTU IDs)
row.names(otu)

#Print column names (which are samples IDs the taxonomy header)
colnames(otu)



### Alpha diversity tables

# Read in the alpha diversity table
alpha <- read.table("alpha_diversity.txt", 
                    sep='\t', 
                    header=TRUE, 
                    as.is=TRUE, 
                    check.names=FALSE,
                    row=1)


  

### Beta diversity tables  

# Load the beta diversity matrix, notice that we use read.table(),
# but then change from a dataframe to a matrix with as.matrix()
beta <- as.matrix(read.table("weighted_unifrac_dm.txt", 
                             sep = "\t", 
                             header=T, 
                             row = 1, 
                             as.is = T, 
                             check.names = F))




### Metadata files   

metadata <- read.table('mouse_mapfile.txt', 
                       header=T, 
                       sep='\t', 
                       check.names=F, 
                       comment='',
                       row=1)

# View the 'Genotype' column of the mapping file dataframe
metadata[,'Genotype']

metadata$Genotype

# The class function will also tell you whether your variable is 
# a factor, numeric, character, etc.
class(metadata[,'Genotype'])
```


######################################################
# Formatting the data
###################################################### 

# First, define all the samples in the OTU table.
# Remember, when we load in the OTU table, samples are columns
# Now let's see what the intersect is with the  metadata row names
IDs_Keep <- intersect(colnames(otu), rownames(metadata))

# Now let's filter the metadata to keep only those samples
# We tell R to make a new data frame that only has the rows we want
metadata <- metadata[IDs_Keep,]

# Now let's filter the OTU table to keep just the intersecting samples

# First, lets rotate the table so the samples are the rows (like the map and alpha)
# We will store it as a new otu table (incase we need the old one)
otu2 <- t(otu)

# This will also remove the taxonomy, because it's not a sample ID we want
otu2 <- otu2[IDs_Keep,]

# Now let's filer the alpha diversity table to keep only those samples
# Alpha diversity has the samples as row names
alpha <- alpha[IDs_Keep, ]

# Now let's filter the beta diversity table to keep only those samples
# Beta diversity has the samples as row names AND column names
# We must filter both the rows and columns
beta <- beta[IDs_Keep,IDs_Keep]

#Let's check to make sure the sample orders match  
rownames(metadata) == rownames(otu2)

#Let's see how many samples are in the otu table and mapping
nrow(otu2)
nrow(metadata)



######################################################
# Alpha diversity differences in R
###################################################### 

### Testing for differences

# Run a test using the shannon index in the alpha diversity table
wilcox.test(alpha$shannon ~ metadata$BMI_Type)

# Note that you would run a t-test in the same way, only using 
# the t.test() function

# Now we run a test using the shannon index in the alpha diversity table
kruskal.test(alpha$shannon ~ metadata$Genotype)

# Note that you would run an ANOVA in the same way, only using 
# the aov() function




### Plotting the data

#We can also see what box plots of the alpha diversity looks like.  

alpha2 <- alpha
alpha2$BMI_Type <- metadata$BMI_Type

library(ggplot2)

ggplot(data=alpha2, aes(x=BMI_Type, y= shannon)) + 
  geom_boxplot(outlier.color = NA) + # removes outlier points becuase we add in the jitter anyways
  geom_jitter(width= 0.1, aes(color=BMI_Type)) +
  theme_bw() +
  guides(color=F) #because the x-axis is already labeled


# We can also print this plot to a pdf with the `pdf()` function followed by `dev.off()` to close the pdf.

plot_output <- ggplot(data=alpha2, aes(x=BMI_Type, y= shannon)) + 
  geom_boxplot(outlier.color = NA) + 
  geom_jitter(width= 0.1, aes(color=BMI_Type)) +
  theme_bw() +
  guides(color=F) 

pdf("Alpha_Diversity.pdf", height=4, width=4)
plot(plot_output)
dev.off()



######################################################
# PCoA in R 
###################################################### 
 
library(ape)
library(vegan)
library(ggplot2)




### Principal coordinates analysis  (PCoA)

# Run the pcoa() function on the beta diversity table,
# and store the vectors generated as a dataframe 
PCOA <- data.frame(pcoa(beta)$vectors)

# If you look at the PCOA table, you'll see the column names 
# are the 'axes' and the row names are sample IDs. We want them to 
# be labeled "PC" instead of "axis"

# We will make a vector with place holders
new_names <- rep("", ncol(PCOA))

# Fill in first with PC followed by the number (e.g. PC1, PC2, PC3...)
for(i in 1:ncol(PCOA)){
  new_names[i] <- paste("PC",i, sep="")
}

# Replace the column names of PCOA
names(PCOA) <- new_names

# Create a column that is SampleIDS for PCOA
PCOA$SampleID <- rownames(PCOA)

#Create a column that is SampleIDs for the metadata
metadata$SampleID <- rownames(metadata)

# Merge the metadata and beta diversity
PCOA <- merge(PCOA, metadata, by = "SampleID")




### Plotting the PCoA  

# Note that geom_point() makes it a scatter plot where the points 
# are colored according to BodySite
ggplot(PCOA) + 
  geom_point(aes(x = PC1, y = PC2, color = Source)) + 
  labs(title="PCoA Plot")

# Now let's add some clusters.  This makes it look great, but can 
# also be misleading and make us think there are groups when there 
# aren't. Note that we are using BodySite to color the points and body 
# AREA to fill the clusters

ggplot(PCOA) + 
  geom_point(aes(x = PC1, y = PC2, color = Source)) + 
  labs(title="PCoA and Clusters") + 
  stat_ellipse(alpha = 0.3, geom="polygon", linetype="blank", aes(x = PC1, y = PC2, fill = Source))



### Changing plotting parameters

colors <- c("darkorchid", "deepskyblue", "darkred", "darkolivegreen")
colors2 <- c("lightblue", "azure3")
ggplot(PCOA, aes(x = PC1, y = PC2)) + 
    stat_ellipse(alpha = 0.5, geom="polygon", aes(fill=Source)) + 
    geom_point(alpha=0.65, size = 3, aes(color = Genotype)) + 
    labs(title="Mouse Microbiome Beta Diversity") + 
    scale_color_manual(values= colors) + 
    scale_fill_manual(values=colors2) + 
    theme(plot.title = element_text(size = 16), 
          axis.title = element_text(size = 14), 
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)) + 
    theme_bw()




### Testing for signifcant differences  

# Turn the beta table into resemblance matrix using as.dist() 
beta_dist = as.dist(beta)

# Test for a significant difference across all groups.  
# This will run an ADONIS test.
ad = adonis(beta_dist ~ metadata[,"Genotype"], data=metadata, permutations=999)
ad

# Now let's write our output to a file.   

# This takes just the analysis of variance table (aoc.tab) 
# from the output
a.table <- ad$aov.tab

# This writes it to a file
write.table(a.table, file="beta_analysis.txt", quote=FALSE, sep="\t", col.names = NA)


