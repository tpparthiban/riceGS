

################################################################################
#############                 GS analysis pipeline                  ############
#############         IRRI Irrigated breeding Program               ############
#############  Parthiban Thathapalli Prakash, Jerome Bartholome     ############
#############                      July 2021                        ############
################################################################################


# The present script is divided into four main sections:
# - Description of the data
# - Training Set Selection for GS
# - Single location analysis
# - Genomic prediction analysis with multi-environment data

# load needed packages ---------------------------------------------------------

# data management and visualization packages
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Genomic prediction / mixed model packages

# We want to draw the attention of the user that the analysis can be done with
# asreml (license required) or with sommer.
# The sommer package is used to compute the genomic relationship matrix

# Choose the package for statistical analysis
# comment or un-comment the following lines to choose the package
# pkg <- "sommer"
pkg <- "asreml"

if (pkg == "asreml") {
  library(asreml)
} else if (pkg == "sommer") {
  library(sommer)
}
library(STPGA)


# Perform training set optimization ? -------------------------------------
# comment or un-comment the following lines to choose to do training set
# optimization via STPGA
# TSO <- FALSE
TSO <- TRUE

# Load data ---------------------------------------------------------------
load("IRRI_GS_data.RData")

# load functions ----------------------------------------------------------
source("IRRI_GS_functions_1.R")

# Description of the data -------------------------------------------------

# Phenotypic data

# The object pheno_data (data frame) contains phenotypic records for the training set (OYT class)
# The training set was evaluated in five location in Bangladesh using p-rep design.
# Three main traits are available : grain yield (t/ha) - ULD_TON
#                                   plant heihgt (cm) - HT_AVG
#                                   days to 50% flowering (days) - FLW50
# 366 lines were evaluated including X lines from the estimation set.
# The others lines are AYT lines and checks.

pheno_data %>% dplyr::select("STUDY",
                             "YLD_TON",
                             "HT_AVG",
                             "FLW50",
                             "MISS HILL") %>%
  tidyr::gather(key = "trait", value = "value", -STUDY) %>%
  group_by(STUDY, trait) %>%
  summarise(nb_obs = n(),
            mean = mean(value, na.rm = T))


# Plotting the data for visual inspection
# a preliminary cleanup step was performed
# to remove values related to input errors

pheno_data %>% dplyr::select(c("STUDY",
                               "YLD_TON",
                               "HT_AVG",
                               "FLW50",
                               "MISS HILL")) %>%
  gather(key = "trait", value = "value", -STUDY) %>%
  ggplot(aes(x = STUDY,
             y = value, )) +
  theme_light(base_size = 11) +
  geom_boxplot() +
  scale_x_discrete("Location") +
  scale_y_continuous('Trait value') +
  theme(
    plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "lines"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  facet_grid(trait ~ ., scales = "free")


# Genotypic data

# The object geno_data (matrix) contains the SNP data of the entire cohort.
# 1079 SNPs are available on 1722 genotypes.
# SNp markers were developed by Arbeleaz et al. 2019.


dim(geno_data)


# Training Set Selection for GS  ------------------------------------------

# For this part we use the package STPGA: Selection of Training Populations by Genetic Algorithm
# The main function in our case is GenAlgForSubsetSelection which is embedded in
# repeatgenalg function which allows for multiple starting points for genetic algorithm

# Important note: The training set that is selected here may be different from the one
# that was actually used in the breeding program. The script is given as an
# example of the strategy that was used.

# This part may take some time to run

if (TSO == T) {
  # Construct the relationship matrix
  geno <- t(geno_data)
  A <- sommer::A.mat(geno)

  res_TS <- OptiTS(A,
                   sTS = 300,
                   rep = 1)

  # Format the data for the graph
  TS <- as.data.frame(res_TS[[1]])
  names(TS)[1] <- paste("Gid")
  TS$Set <- "Training set"

  pc <- prcomp(x = A, scale. = F)
  pcs <- as.data.frame(pc$x)
  variance <-  summary(pc)$importance[2, ]
  cumvar <- cumsum(variance)
  pcs <- pcs[, 1:5]
  pcs <- data.frame("Gid" = as.character(row.names(pcs)), pcs)

  pcs_all <- merge(pcs, TS, by = "Gid", all.x = TRUE)
  pcs_all$Set[is.na(pcs_all$Set)] <- "Predicted set"

  # Plot the PCs to see the distribution of the training set

  ggplot(pcs_all, aes(x = PC1, y = PC2, group = Set)) +
    geom_point(aes(shape = Set,
                   color = Set),
               size = 2,
               alpha = 0.7) +
    theme_light(base_size = 12) +
    scale_color_manual(values = c("Training set" = "red", "Predicted set" = "grey50")) +
    xlab("PC1 (25.9%)") + ylab("PC2 (24.7%)") +
    theme(legend.position = c(0.85, 0.85),
          legend.title = element_blank())
}


# Single location analysis ------------------------------------------------

# Apply single trial function - for multiple trials and multiple traits

traits <- c("YLD_TON", "FLW50", "HT_AVG")

if (pkg == "asreml") {
  asreml.options(trace = F)
}

print(pkg)

res_list <- list()
model_list <- list()
a <- 0

# Loop on the locations
for (loc in unique(pheno_data$STUDY_ID)) {
  # Loop on the traits
  for (trait in traits) {
    a <- a + 1
    #subset the phenotypic data
    dataset <- pheno_data %>% dplyr::filter(STUDY_ID == loc) %>%
      dplyr::select(c(
        "LOCATION",
        "GID",
        "DESIGN",
        "REP",
        "DESIGN_X",
        "DESIGN_Y",
        trait
      )) %>%
      droplevels()

    colnames(dataset)[colnames(dataset) == trait] <-
      deparse(substitute(trait))

    cat("\n")
    print(paste("Location:", unique(dataset$LOCATION)))
    print(paste("Trait:", trait))
    # Check if too many data are missing
    if (sum(is.na(dataset)) > 100)
    {
      model <- NA
    }
    else
    {
      if (pkg == "asreml") {
        res <- single_trial_asreml(dataset)
      } else if (pkg == "sommer") {
        res <- single_trial_sommer(dataset)
        res[[2]]$gid <- gsub("GID", "", res[[2]]$gid)
      }

      res_list[[a]] <- res[[2]]
      model_list[[a]] <- res[[1]]
    }
  }
}

# Combine the output
res <- dplyr::bind_rows(res_list)


# Plot the broad sense heritability for all trials and traits

res %>% group_by(location, trait) %>%
  summarise(h2 = unique(h2)) %>%
  ggplot() +
  theme_light(base_size = 12) +
  geom_bar(aes(x = location,
               y = h2),
           stat = "identity") +
  scale_x_discrete("Location") +
  scale_y_continuous('Heritability') +
  theme(
    plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "lines"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  facet_grid(trait ~ .,)

# Boxplot of the BLUPs for each location

ggplot(res, aes(x = location,
                y = blups)) +
  theme_light(base_size = 12) +
  geom_boxplot() +
  scale_x_discrete("Location") +
  scale_y_continuous('BLUP value') +
  theme(
    plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "lines"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  ) +
  facet_grid(trait ~ ., scales = "free")

# Boxplot of the de-regressed BLUPs for each location

ggplot(res, aes(x = location,
                y = dblup)) +
  theme_light(base_size = 12) +
  geom_boxplot() +
  scale_x_discrete("Location") +
  scale_y_continuous('De-regressed BLUP value') +
  theme(
    plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "lines"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  ) +
  facet_grid(trait ~ ., scales = "free")



#  Multi location analysis  -----------------------------------------------

# Estimate G matrix from markers

G <- sommer::A.mat(t(geno_data))
diag(G) <- diag(G) + 0.001

# asreml package requires to compute the inverse of the G matrix
if (pkg == "asreml") {
  invG <- solve(G)
  #check the attribute rowNames
  attr(G, "rowNames") <- row.names(G)
  attr(invG, "rowNames") <- row.names(G)
  attr(invG, "INVERSE") <- TRUE
}

# Format the results from the single trial analysis
res_s <- res %>% dplyr::select(c("location",
                                 "gid",
                                 "trait",
                                 "blups",
                                 "dblup")) %>%
  pivot_wider(names_from = trait, values_from = c(blups, dblup)) %>%
  filter(gid %in% rownames(G)) %>%
  mutate(breedingzone = "Bangladesh") %>%
  mutate(location = factor(location)) %>%
  mutate(breedingzone = factor(breedingzone)) %>%
  mutate(gid = factor(gid)) %>%
  droplevels()


# Start the GS analysis by selecting the traits
traits <- c("dblup_YLD_TON", "dblup_FLW50", "dblup_HT_AVG")

if (pkg == "asreml") {
  asreml.options(trace = F)
}

print(pkg)

res_list <- list()
model_list <- list()
a <- 0

# Loop on the breeding zones
for (loc in unique(res_s$breedingzone)) {
  # Loop on the traits
  for (trait in traits) {
    a <- a + 1
    dataset <- res_s %>% filter(breedingzone == loc) %>%
      dplyr::select(c("location",
                      "breedingzone",
                      "gid",
                      trait)) %>%
      droplevels()

    colnames(dataset)[colnames(dataset) == trait] <-
      deparse(substitute(trait))

    cat("\n")
    print(paste("Breedingzone:", unique(dataset$breedingzone)))
    print(paste("Trait:", trait))

    if (sum(is.na(dataset)) > 1500)
    {
      model <- NA
    }
    else
    {
      if (pkg == "asreml") {
        res_bv <- gblup_asreml(dataset, invG)
      } else if (pkg == "sommer") {
        # add the other part of the cohort to the training set
        dataset <-
          dataset %>%  full_join(tibble(
            gid = rep(as.factor(rownames(G)), 5),
            location = as.factor(c(
              rep("Loc1", nrow(G)),
              rep("Loc2", nrow(G)),
              rep("Loc3", nrow(G)),
              rep("Loc4", nrow(G)),
              rep("Loc5", nrow(G))
            ))
          ),
          by = c("location", "gid")) %>%
          droplevels()
        res_bv <- gblup_sommer(dataset, G)
      }
      res_list[[a]] <- res_bv[[2]]
      model_list[[a]] <- res_bv[[1]]
    }
  }
}

res_gebv <- dplyr::bind_rows(res_list)

# Plot the narrow sense heritability for all traits in a breeding zone

res_gebv %>% group_by(trait) %>%
  summarise(h2 = unique(h2)) %>%
  ggplot() +
  theme_light(base_size = 11) +
  geom_bar(aes(x = trait,
               y = h2),
           stat = "identity") +
  scale_x_discrete("Trait") +
  scale_y_continuous('Heritability') +
  theme(
    plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "lines"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



# Multi location summary file
agr_traits <- res_gebv %>% dplyr::select(c("breedingzone",
                                           "gid",
                                           "trait",
                                           "BV",
                                           "rel")) %>%
  pivot_wider(names_from = trait, values_from = c(BV, rel))



# End --------------------------------------------------------------------------
