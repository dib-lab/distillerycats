library(readr)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(ggpubr)

# read in comp and format  ----------------------------------------------
comp <- read_csv(snakemake@input[["comp"]])         # read in mat
colnames(comp) <- gsub("_filt", "", colnames(comp)) # remove _filt from name
rownames(comp) <-colnames(comp)                     # set rownames as samples


# read in and formate metadata -------------------------------------------

info <- read_tsv(snakemake@input[["info"]]) %>%
  filter(library_name %in% colnames(comp)) %>%
  group_by(library_name, study_accession, diagnosis, subject) %>%
  summarise(read_count = sum(read_count)) %>%
  as.data.frame()

# cacl vars for plotting ---------------------------------------------------

dist <- dist(comp)                                              # calc dist on comp mat
fit_all <- cmdscale(dist, eig = T)                              # calculate MDS
fit <- as.data.frame(fit_all$points)                            # abstract data
fit$sample <- rownames(fit)                                     # set rownames
colnames(fit) <- c("dim1", "dim2", "sample")                    # set column names
fit <- left_join(fit, info, by = c("sample" = "library_name"))  # join with metadata
var <- round(fit_all$eig*100/sum(fit_all$eig), 1)               # calc percent var


# plot by study  ------------------------------------------------------------

# make base plot
study_plt <- ggplot(fit, aes(x = dim1, y = dim2, color = study_accession, shape = diagnosis)) +
  geom_point() +
  theme_minimal() +
  ggtitle("IBD studies") +
  labs(x = paste0("PCo 1 (", var[1], "%)"),
       y = paste0("PCo 2 (", var[2], "%)")) +
  scale_shape(solid = F) +
  scale_color_brewer(palette = "Paired") +
  theme(plot.title = element_text(hjust = 0.5))

# remove the legend, as ggmarginal plots outside of legend
study_plt_final <- study_plt + theme(legend.position = "none")
# add histograms to edges of plt
study_plt_final <- ggExtra::ggMarginal(study_plt_final, type = "density",
                                       groupColour = T, groupFill = T)
# generate the legend from the original plot as an object
legend <- study_plt + 
  theme(legend.text=element_text(size=12),
        legend.title=element_text(size=13))
legend <- as_ggplot(get_legend(legend))

# join legend with plot
study_plt_final <- ggarrange(as_ggplot(study_plt_final), legend,
                             ncol = 2, heights = 3, widths = c(2, 1))

pdf(snakemake@output[['study']], height = 6, width = 9)
study_plt_final
dev.off()

# plot by diagnosis ------------------------------------------------------------

# make base plot
diagnosis_plt <- ggplot(fit, aes(x = dim1, y = dim2, shape = study_accession, color =diagnosis)) +
  geom_point() +
  theme_minimal() +
  ggtitle("IBD studies") +
  labs(x = paste0("PCo 1 (", var[1], "%)"),
       y = paste0("PCo 2 (", var[2], "%)")) +
  scale_shape(solid = F, name = "Study") +
  scale_color_viridis_d(name = "Diagnosis") + 
  theme(plot.title = element_text(hjust = 0.5))

# remove the legend, as ggmarginal plots outside of legend
diagnosis_plt_final <- diagnosis_plt + theme(legend.position = "none")
# add histograms to edges of plt
diagnosis_plt_final <- ggExtra::ggMarginal(diagnosis_plt_final, type = "density",
                                           groupColour = T, groupFill = T)
# generate the legend from the original plot as an object
legend <- diagnosis_plt + 
  theme(legend.text=element_text(size=12),
        legend.title=element_text(size=13))
legend <- as_ggplot(get_legend(legend))

# join legend with plot
diagnosis_plt_final <- ggarrange(as_ggplot(diagnosis_plt_final), legend,
                                 ncol = 2, heights = 3, widths = c(2, 1))

pdf(snakemake@output[['diagnosis']], height = 6, width = 9)
diagnosis_plt_final
dev.off()
