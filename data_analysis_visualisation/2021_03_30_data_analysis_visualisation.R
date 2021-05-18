# ================ #
# ---- METADATA ----
# ================ #

# Schenkel MA, Beukeboom LW & Pen I (2020). Epistatic interactions between sex chromosomes and 
# autosomes can affect the stability of sex determination systems. [JOURNAL, VOLUME, PAGES, DOI]

# Contact:
# Martijn Schenkel (m.a.schenkel@rug.nl; maschenkel@gmail.com)
# Ido Pen (i.r.pen@rug.nl)

# Application: 
# * Collect generated data from simulations into single large data set
# * Process data
# * Visualise data

# =========================== #
# ----Administrative setup ----
# =========================== #

# Clear workspace and load required packages
rm(list = ls())
library(tibble)
library(tidyverse)
library(viridis)
library(cowplot)
library(mgcv)

# set work directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Variables related to simulation run
Date <- "2020_08_03" # Date of simulation
Types <- c("COA", "DXD", "DOM") # Abbreviated epistasis types
Transition <- c("YM_AM", "YM_FD") # SD transitions (Y-A, Y-W)

# Range of simulations (N = 1000)
start <- 0
end   <- 1999

# Tracker for failed simulations
failed <- NA

# ======================= #
# ---- Data collection ----
# ======================= #

# Loop over epistasis types, SD transition types
for(ty in 1:length(Types)) for(tr in 1:length(Transition)) 
{
  # Define base tibble.
  d <- tibble(Seed = numeric(),
              Generation = numeric(),
              SexRatio = numeric(),
              Frequency = numeric(),
              Sex = character(),
              LinkageGroup = character(),
              Copy = character(),
              Allele = character(),
              Type = character(),
              SAY = numeric(),
              SAI = numeric(),
              SAII = numeric(),
              Epsilon = numeric(),
              Scenario = character(),
              Epistasis = character())
  
  # Determine which SA loci can interwact with EPI, SAY or SAI (SAA) for Y-A and SAY or SAII (SAW) for Y-W
  if(tr == 1){ Interact <- c("SAY", "SAI") }
  if(tr == 2){ Interact <- c("SAY", "SAII") }
  
  # Now loop over interacting loci
  for(i in 1:length(Interact))
  { 
    # Print current epistasis type + SD type + interacting locus
    print(paste(Types[ty], "_", Transition[tr], "_", Interact[i], sep = ""))
    
    # Loop independent simulations (n = 0 to n = 999)
    for(n in start:end)
    {
      # Define current folder holding data file (= simulation ID)
      dir <- paste0("home/p275703/2019_04_SD_epistasis/", Date, "/", Types[ty], "_", Transition[tr], "_", Interact[i], "_", n, "/")

      # If data file exists, load and add to "overall" data tibble d
      # If not, add simulation ID to list of failed simulations
      if(file.exists(paste0(dir, n, ".txt")))
      {
      dt <- read.table(file = paste(dir, n, ".txt", sep = ""),
                      header = T)

      d <- dt %>% full_join(d,  by = c("Seed", "Generation", "SexRatio", 
                                       "Frequency", "Sex", "LinkageGroup", 
                                       "Copy", "Allele", "Type", 
                                       "SAY", "SAI", "SAII", 
                                       "Epsilon", "Scenario", "Epistasis"))
      }
      else {
        failed <- c(failed, dir) # if failed, add to list of failed simulations
      }
    }
  }
  # Write data for epistasis type + SD transition to output. 
  write.table(x = d, file = paste0(Date, "_", Types[ty], "_", Transition[tr], ".txt"), row.names = F, col.names = T)
}
# warnings can be ignored, factor -> character conversion

# Define an overall tibble for ALL simulations
d2 <- tibble(Seed = numeric(),
            Generation = numeric(),
            SexRatio = numeric(),
            Frequency = numeric(),
            Sex = character(),
            LinkageGroup = character(),
            Copy = character(),
            Allele = character(),
            Type = character(),
            SAY = numeric(),
            SAI = numeric(),
            SAII = numeric(),
            Epsilon = numeric(),
            Scenario = character(),
            Epistasis = character())

# Put together data from all epistasis + SD transition considered.
for(ty in 1:length(Types)) for(tr in 1:length(Transition)) 
{
  dt <- read.table(paste0(Date, "_", Types[ty], "_", Transition[tr], ".txt"), T)
  d2 <- dt %>% full_join(d2)
}

write.table(x = d2, file = paste0(Date, "_data_collected_full.txt"), row.names = F, col.names = T)

# ======================= #
# ---- Data processing ----
# ======================= #

# Redefine tibble d as the one containing ALL data. 
d <- read.table(paste0(Date, "_data_collected_full.txt"), T)

# Correct reported frequencies by dividing by frequency of sex in which they are reported
d <- d %>% mutate(FreqCor = if(Sex == "M"){ Frequency / SexRatio } else { Frequency / (1 - SexRatio) },
                  Epistasis = factor(Epistasis),
                  Scenario = factor(Scenario),
                  ScenEpi = factor(paste(Scenario, Epistasis)))
with(d, table(ScenEpi))

# Split into dataset with epsilon-values above zero (randomly sampled) and those that are zero (independently set).
# We fit separate GAMs on these dataset a bit further down, because the original analysis showed that the fit for the
# dataset with non-zero values tends to be unable to deal with the qualitative difference between epsilon = 0 and 
# epsilon != 0.

dNonZero <- d %>% filter(Epsilon != 0)
dZero    <- d %>% filter(Epsilon == 0)

with(dZero, table(ScenEpi))
# For Y-A transitions:

# Filter out A (YM) frequencies on paternal copy in males and round to 0 or 1. Tiny deviations from these values 
# can exist where the frequency of an SD gene is very close to 0 or 1 (i.e. delta = 10^-16), which messes up 
# binomial analysis later on in the script.
dA <- dNonZero %>% filter(Sex == "M", Allele == "AM", Copy == "Paternal",
                   Scenario %in% c("YM_AM_SAY", "YM_AM_SAI")) %>% mutate(FC = round(FreqCor))

# For Y-W transitions: 

# Filter out W (FD) frequencies on maternal copies in females and round to 0 or 1 again.
dF <- dNonZero %>% filter(Sex == "F", Allele == "FD", Copy == "Maternal",
                   Epistasis != 0,
                   Scenario %in% c("YM_FD_SAY", "YM_FD_SAII")) %>% mutate(FC = round(FreqCor))


# Repeated for datasets with epsilon = 0
dAZ <- dZero %>% filter(Sex == "M", Allele == "AM", Copy == "Paternal",
                          Scenario %in% c("YM_AM_SAY", "YM_AM_SAI")) %>% mutate(FC = round(FreqCor))

# For Y-W transitions: 

# Filter out W (FD) frequencies on maternal copies in females and round to 0 or 1 again.
dFZ <- dZero %>% filter(Sex == "F", Allele == "FD", Copy == "Maternal",
                          Scenario %in% c("YM_FD_SAY", "YM_FD_SAII")) %>% mutate(FC = round(FreqCor))

# ===================== #
# ---- Data analysis ----
# ===================== #

# Raw data plots.
dA %>% ggplot(aes(SAY, SAI, color = FC)) + geom_point() + facet_grid(cut(Epsilon, breaks = 5)~Epistasis~Scenario)
dF %>% ggplot(aes(SAY, SAII, color = FC)) + geom_point() + facet_grid(cut(Epsilon, breaks = 5)~Epistasis~Scenario)
dAZ %>% ggplot(aes(SAY, SAI, color = FC)) + geom_point() + facet_grid(Epistasis~Scenario)
dFZ %>% ggplot(aes(SAY, SAII, color = FC)) + geom_point() + facet_grid(Epistasis~Scenario)

# Define GAMs. Used here is a tensor smooth for SAY, SAI or SAII, and Epsilon, estimated per group (SA locus involved * epistasis type) 
# with a model I setup cf. Pedersen et al. (no common trend, group-specific wiggliness)

# Pedersen EJ, Miller DL, Simpson GL, Ross N. 2019. Hierarchical generalized additive models in ecology: an introduction with mgcv. 
# PeerJ 7:e6876 https://doi.org/10.7717/peerj.6876

# GAM for Y-A transitions. We analyze the presence/absence (1/0) of A on the paternal allele in males.
# if present, A is the SD gene; if absent, Y is the SD gene.
GA <- gam(FC ~ ScenEpi + t2(SAY, SAI, Epsilon, k = 4, m = 2, bs = "ts", by = ScenEpi),
          data = dA,  full = T, family = "binomial", file = paste0(Date,"_GAM_YA.rds"))

# GAM for Y-W transitions. We analyze the presence/absence (1/0) of W on the maternal allele in females.
# if present, W is the SD gene; if absent, Y is the SD gene.
GF <- gam(FC ~ ScenEpi + t2(SAY, SAII, Epsilon, k = 4, m = 2, bs = "ts", by = ScenEpi),
          data = dF,  full = T, family = "binomial")

# REPEAT FOR EPSILON = 0 DATASET #

# GAM for Y-A transitions. We analyze the presence/absence (1/0) of A on the paternal allele in males.
# if present, A is the SD gene; if absent, Y is the SD gene.
GAZ <- gam(FC ~ ScenEpi + t2(SAY, SAI, k = 4, m = 2, bs = "ts", by = ScenEpi),
          data = dAZ,  full = T, family = "binomial")

# GAM for Y-W transitions. We analyze the presence/absence (1/0) of W on the maternal allele in females.
# if present, W is the SD gene; if absent, Y is the SD gene.
GFZ <- gam(FC ~ ScenEpi + t2(SAY, SAII, k = 4, m = 2, bs = "ts", by = ScenEpi),
          data = dFZ,  full = T, family = "binomial")

# Save GAM fits
saveRDS(object = GA, file = paste0(Date,"_GAM_YA.rds"))
saveRDS(object = GF, file = paste0(Date,"_GAM_YW.rds"))
saveRDS(object = GAZ, file = paste0(Date,"_GAM_YAZ.rds"))
saveRDS(object = GFZ, file = paste0(Date,"_GAM_YWZ.rds"))

# Load GAM fits
GA <- readRDS(paste0(Date,"_GAM_YA.rds"))
GF <- readRDS(paste0(Date,"_GAM_YW.rds"))
GAZ <- readRDS(paste0(Date,"_GAM_YAZ.rds"))
GFZ <- readRDS(paste0(Date,"_GAM_YWZ.rds"))

# k-value validation
gam.check(GA) # Good, only p<0.05 for SAI (SA^A) + coadaptation. Still, edf << k, so should not be
              # problematic.Fitted values in graphs later on also seem to fit data quite well.
gam.check(GAZ) # Great
gam.check(GF) # Great
gam.check(GFZ) # Great
# Generate dataframes with predictions
dpa <- expand.grid(Scenario = levels(as.factor(dA$Scenario))[1:2],
                   Epistasis = levels(as.factor(dA$Epistasis)),
                   SAY = seq(0, 0.05, length.out = 200),
                   SAI = seq(0, 0.05, length.out = 200),
                   Epsilon = seq(0.01, 0.05, length.out = 5))
dpa <- dpa %>% mutate(ScenEpi = paste(Scenario, Epistasis))
dpa <- dpa %>% transform(Frequency = predict(GA, newdata = dpa, type = "response"))

# epsilon = 0 predictions
dpaZ <- expand.grid(Scenario = levels(as.factor(dAZ$Scenario))[1:2],
                    Epistasis = levels(as.factor(dAZ$Epistasis)),
                    SAY = seq(0, 0.05, length.out = 200),
                    SAI = seq(0, 0.05, length.out = 200))
dpaZ <- dpaZ %>% mutate(ScenEpi = paste(Scenario, Epistasis))
dpaZ <- dpaZ %>% transform(Frequency = predict(GAZ, newdata = dpaZ, type = "response"))

# Repeated for Y-W transitions
dpf <- expand.grid(Scenario = levels(as.factor(dF$Scenario))[3:4],
                   Epistasis = levels(as.factor(dF$Epistasis)),
                   SAY = seq(0, 0.05, length.out = 200),
                   SAII = seq(0, 0.05, length.out = 200),
                   Epsilon = seq(0.01, 0.05, length.out = 5))
dpf <- dpf %>% mutate(ScenEpi = paste(Scenario, Epistasis))
dpf <- dpf %>% transform(Frequency = predict(GF, newdata = dpf, type = "response"))

# epsilon = 0 predictions
dpfZ <- expand.grid(Scenario = levels(as.factor(dFZ$Scenario))[3:4],
                   Epistasis = levels(as.factor(dFZ$Epistasis)),
                   SAY = seq(0, 0.05, length.out = 200),
                   SAII = seq(0, 0.05, length.out = 200))
dpfZ <- dpfZ %>% mutate(ScenEpi = paste(Scenario, Epistasis))
dpfZ <- dpfZ %>% transform(Frequency = predict(GFZ, newdata = dpfZ, type = "response"))


dpa <- dpaZ %>% mutate(Epsilon = 0) %>% full_join(dpa) %>% as_tibble()
dpf <- dpfZ %>% mutate(Epsilon = 0) %>% full_join(dpf) %>% as_tibble()
# Save generated predictions to file, it takes quite a while to run this stuff.
write.table(x = dpa, file = paste0(Date, "_AM_predicted_freq_YM_AM.txt"), row.names = F, col.names = T)
write.table(x = dpf, file = paste0(Date, "_FD_predicted_freq_YM_FD.txt"), row.names = F, col.names = T)


# ================ #
# ---- Graphics ----
# ================ #

# # Define strip labels (ugly but it works).

dpa$fS[dpa$Scenario=="YM_AM_SAY"] <- "SA^Y"
dpa$fS[dpa$Scenario=="YM_AM_SAI"] <- "SA^A"

dpa$fE[dpa$Epistasis=="COA"] <- "Coadaptation"
dpa$fE[dpa$Epistasis=="DOM"] <- "Dominance"
dpa$fE[dpa$Epistasis=="DXD"] <- "Overdominance"


dpf$fS[dpf$Scenario=="YM_FD_SAY"] <- "SA^Y"
dpf$fS[dpf$Scenario=="YM_FD_SAII"] <- "SA^W"

dpf$fE[dpf$Epistasis=="COA"] <- "Coadaptation"
dpf$fE[dpf$Epistasis=="DOM"] <- "Dominance"
dpf$fE[dpf$Epistasis=="DXD"] <- "Overdominance"

# Reorder loci labels so that SA^Y is first, this puts it on the top row in the plots.
# dpy <- dpy %>% mutate(fS = factor(fS, levels = c("SA^Y", "SA^A")))
dpa <- dpa %>% mutate(fS = factor(fS, levels = c("SA^Y", "SA^A")))
dpf <- dpf %>% mutate(fS = factor(fS, levels = c("SA^Y", "SA^W")))

# Define plot function. Note that we use geom_contour to draw a contour line for those parameter 
# values where the GAM predicts the frequency will be 1/2. This is in essence the transition 
# boundary between Y as the SD gene and A as the SD gene. We know from the earlier raw data plots 
# that the frequency of A below this line is <1/2 and above the line >1/2, so that we can annotate 
# it with some labels to illustrate which areas in parameter space yield a transition and which don't.

baseplot <- function(p){
  
  p +  
  
  # Draw contour line at frequency = 1/2.
  geom_contour(aes(z = Frequency, color = factor(Epsilon)), binwidth = 1/2, size = 1) + 
    
    # Generate separate facet for each combination of interacting SA locus + epistasis type
    facet_grid(fS~fE, labeller = label_parsed) + 
    
    # X-Y axes. Note that we label with 1-5 instead of 0.01-0.05, which makes it easier to read.
    scale_x_continuous(expand = c(0,0), limits = c(0,0.05), breaks = seq(0, 0.05, 0.01), labels = 0:5) + 
    scale_y_continuous(expand = c(0,0), limits = c(0,0.05), breaks = seq(0, 0.05, 0.01), labels = 0:5) +
    
    # Overall plot makeup
    theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2), color = rgb(0.2,0.2,0.2)),
          strip.text = element_text(color = "white"),
          plot.margin = margin(0,0,0,0, "lines"),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0.75, "lines"),
          text = element_text(size = 10),
          legend.position = "right",
          legend.box.margin = margin(0,0,0,0),
          legend.margin = margin(1,1,1,1)) +
    scale_color_viridis(option = "magma", begin = 0.1, end = 0.9, discrete = T) +
    
    # Define legend.
    guides(color = guide_legend(title.position = "top",
                                title.hjust = 0.5,
                                ncol = 1)) +
    labs(x = bquote("Effect of SA"^"Y"~"("*italic(S)[XY]%*%10^2*")"),
         # y = bquote("Effect of SA"^"A"~"("*italic(S)[A]%*%10^2*")"),
         color = bquote("Strength of\nepistasis ("*epsilon*")"))
  
}

# Generate plots.

pya <- dpa %>%
  ggplot(aes(SAY, SAI, fill = Epsilon)) %>% 
  baseplot() + 
  
  # Axis/legend titles
  labs(y = bquote("Effect of SA"^"A"~"("*italic(S)[A]%*%10^2*")")) +
  
  # Annotate which areas have Y vs. which have A as the SD gene
  annotate("text", label = "A", parse = T, x = 0.01, y = 0.04, color = rgb(0.2,0.2,0.2)) + 
  annotate("text", label = "Y", parse = T, x = 0.04, y = 0.01, color = rgb(0.2,0.2,0.2))

# Y->W transition
pyf <- dpf %>% 
  ggplot(aes(SAY, SAII, fill = Epsilon)) %>% 
  baseplot() + 
  
  # Axis/legend titles
  labs(y = bquote("Effect of SA"^"W"~"("*italic(S)[W]%*%10^2*")")) + 
  
  # Annotate which areas have Y vs. which have A as the SD gene
  annotate("text", label = "W", parse = T, x = 0.01, y = 0.04, color = rgb(0.2,0.2,0.2)) + 
  annotate("text", label = "Y", parse = T, x = 0.04, y = 0.01, color = rgb(0.2,0.2,0.2))

# Save to PDF
pdf(file = paste0(Date, "_YM_AM_transition_(Figure_2).pdf"), height = 4, width = 6)
pya
dev.off()

pdf(file = paste0(Date, "_YM_FD_transition_(Figure_3).pdf"), height = 4, width = 6)
pyf
dev.off()


# Save to PNG
png(file = paste0(Date, "_YM_AM_transition_(Figure_2).png"), height = 4, width = 6, res = 1000, units = "in")
pya
dev.off()

png(file = paste0(Date, "_YM_FD_transition_(Figure_3).png"), height = 4, width = 6, res = 1000, units = "in")
pyf
dev.off()

dpa <- dpa %>% mutate(fEps = paste0("epsilon == ", Epsilon))
pya_coa <- dpa %>% filter(Epistasis == "COA") %>% 
  ggplot(aes(SAY, SAI, fill = Frequency > 1/2)) + 
  geom_tile() +
  scale_fill_viridis(option = "magma", begin = 0.05, end = 0.95, discrete = T, labels = c(0, 1)) +
  facet_grid(fS ~ fEps, labeller = label_parsed) + 
  labs(x = bquote("Effect of SA"^"Y"~"("*italic(S)[XY]%*%10^2*")"),
       y = bquote("Effect of SA"^"A"~"("*italic(S)[A]%*%10^2*")"),
       fill = "Frequency of A") + 
  
  # X-Y axes. Note that we label with 1-5 instead of 0.01-0.05, which makes it easier to read.
  scale_x_continuous(expand = c(0,0), limits = c(0,0.05), breaks = seq(0, 0.05, 0.01), labels = 0:5) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,0.05), breaks = seq(0, 0.05, 0.01), labels = 0:5) + 
  
  # Overall plot makeup
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2), color = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        plot.margin = margin(0,0,0,0, "lines"),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.75, "lines"),
        text = element_text(size = 10),
        
        legend.position = "bottom",
        legend.box.margin = margin(0,0,0,0),
        legend.margin = margin(1,1,1,1)) + 
  guides(fill = guide_legend(title.position = "top",
                               title.hjust = 0.5))


dpf <- dpf %>% mutate(fEps = paste0("epsilon == ", Epsilon))
pyf_coa <- dpf %>% filter(Epistasis == "COA") %>% 
  ggplot(aes(SAY, SAII, fill = Frequency > 1/2)) + 
  geom_tile() +
  scale_fill_viridis(option = "magma", begin = 0.05, end = 0.95, discrete = T, labels = c(0, 1)) +
  facet_grid(fS ~ fEps, labeller = label_parsed) + 
  labs(x = bquote("Effect of SA"^"Y"~"("*italic(S)[XY]%*%10^2*")"),
       y = bquote("Effect of SA"^"W"~"("*italic(S)[W]%*%10^2*")"),
       fill = "Frequency of W") + 
  
  # X-Y axes. Note that we label with 1-5 instead of 0.01-0.05, which makes it easier to read.
  scale_x_continuous(expand = c(0,0), limits = c(0,0.05), breaks = seq(0, 0.05, 0.01), labels = 0:5) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,0.05), breaks = seq(0, 0.05, 0.01), labels = 0:5) + 
  
  # Overall plot makeup
  theme(strip.background = element_rect(fill = rgb(0.2,0.2,0.2), color = rgb(0.2,0.2,0.2)),
        strip.text = element_text(color = "white"),
        plot.margin = margin(0,0,0,0, "lines"),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.75, "lines"),
        text = element_text(size = 10),
        legend.position = "bottom",
        legend.box.margin = margin(0,0,0,0),
        legend.margin = margin(1,1,1,1)) + 
  guides(fill = guide_legend(title.position = "top",
                             title.hjust = 0.5))

# Save to PDF
pdf(file = paste0(Date, "_YM_AM_transition_COA_(SuppFig1).pdf"), height = 4, width = 6)
pya_coa
dev.off()

pdf(file = paste0(Date, "_YM_FD_transition_COA_(SuppFig2).pdf"), height = 4, width = 6)
pyf_coa
dev.off()

# Save to PNG
png(file = paste0(Date, "_YM_AM_transition_COA_(SuppFig1).png"), height = 4, width = 6, res = 1000, units = "in")
pya_coa
dev.off()

png(file = paste0(Date, "_YM_FD_transition_COA_(SuppFig2).png"), height = 4, width = 6, res = 1000, units = "in")
pyf_coa
dev.off()
