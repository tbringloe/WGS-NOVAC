# Info --------------------------------------------------------------------
#
# Overview: Analysis of ddRAD and lcWGS datasets in Delphinapterus leucas (Beluga)
# 
# Author: Audrey Bourret, modifications by Trevor Bringloe and Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory
# Location: Maurice Lamontagne Institute
# Date: 2023-05-25
#

# Library -----------------------------------------------------------------
library(here)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(viridis)

# PCA plots --------------------------------------------------------------

#load and plot eigenvec values
eigenval <- scan(file.path(here::here(), "02_Results/01b_ddRAD_Bringloe/01_PCA/", "MOBELS_S_20_00703_Beluga_SNPf_final_27iv23.eigenval"))

#get percentages of variation explained by vectors
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

#plot to visualise variance explained by vector
gg.ein <- ggplot(pve[tail(order(pve$pve), 20), ], aes(PC, pve, fill=PC)) + geom_bar(stat="identity") + geom_line(linetype="dashed") + scale_fill_viridis() + 
  ylab("Percentage variance explained") + xlab("PC axis") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position = "none", axis.text.x=element_blank(), axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18)) + 
  theme(axis.text = element_text(size = 18, colour = "black"), axis.title.y.right = element_text(angle = 90), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1))

ggsave(filename = file.path(here::here(), "02_Results/01b_ddRAD_Bringloe/01_PCA/", "lcWGS_PCA_eigenvalues.png"), 
       plot = gg.ein,
       height = 3.5, width = 4, units = "in")   

#plot pca with % values on axes
pca <- read_csv(file.path(here::here(), "02_Results/01b_ddRAD_Bringloe/01_PCA/", "MOBELS-S_20_00703-Beluga-lcWGS-PCA-9v23.csv"), col_names = TRUE)
pca$Region <- factor(pca$Region, levels = c("Saint Lawrence Estuary","Cumberland Sound","Frobisher Bay","Ungava Bay","South Hudson Strait","North Hudson Strait","North-East Hudson Bay","Resolute Passage","South-East Hudson Bay","Belcher Islands","James Bay","Western Hudson Bay"))
pca$Season <- factor(pca$Season, levels = c("Spring","Summer","Fall","Winter","Unknown"))

gg.pca <- ggplot(pca, aes(PC1, PC2, col = Region, shape = Season, size = Season)) + geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(stroke = 1, alpha=0.65) + # added stroke for presentation
  scale_color_manual(name = "Harvest region", values = c("deepskyblue","springgreen4","plum2","purple3","orchid3","palevioletred2","salmon","lightblue","red3","chocolate3","orange","royalblue1","grey"),
                     labels = c("Saint Lawrence estuary", "Cumberland Sound", "Frobisher Bay", "Ungava Bay", "South Hudson Strait", "North Hudson Strait", "North-East Hudson Bay", "Resolute Passage", "South-East Hudson Bay", "Belcher Islands", "James Bay", "Western Hudson Bay")) +
  scale_size_manual(values = c(3,7.5,3,3,3)) + # 3 for spring-fall and NA, 9 for summer
  scale_shape_manual(values = c(6,19,2,0,5)) + # reversed triangle Spring, dots Summer, tringle Fall, squadre Winter, losange Unknown
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = guide_legend(title = "Harvest region",
                               label.position = "right",
                               title.position = "top", title.hjust = 0,
                               override.aes = list(size = 5),
                               order = 1),
         size = guide_legend(title = "Harvest season"),
         shape = guide_legend(title = "Harvest season")) +
  theme(legend.position = c(0.85,0.7),
        legend.spacing.y = unit(0.25, "cm"),
        legend.title = element_text(size = 17, face = "bold"),
        legend.text = element_text(size = 16),
        legend.background = element_blank()) +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 1), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 2), "%)")) +
  labs(colour = "Harvest region") + guides(shape = FALSE) + guides(size=FALSE)
#if plotting regions on seperate plots
facet(gg.pca, facet.by = "Region") + guides(colour = FALSE) + theme(legend.position = c(0.32,0.1))

ggsave(filename = file.path(here::here(), "02_Results/01b_ddRAD_Bringloe/01_PCA/", "lcWGS_PCA_1+2.png"), 
       plot = gg.pca,
       height = 3.5, width = 4, units = "in")   

# Dataset -----------------------------------------------------------------

# population data
pop.data <- read_csv(file.path(here::here(), "02_Results/01b_ddRAD_Bringloe/02_Admixture", "beluga_lcWGS_2v23.csv"), col_names = TRUE)

# Admixture ---------------------------------------------------------------

# Original SNP
bed.file <- file.path(here::here(), "02_Results/01b_ddRAD_Bringloe/02_Admixture", "MOBELS_S_20_00703_Beluga_SNPf_final_27iv23.bed")
file.exists(bed.file)
fam.file <- bed.file %>% str_replace(".bed", ".fam")
fam <- read.table(fam.file)

#Code to run ADMIXTURE
for(k in 1:10){
  
  print(k)  
  
  setwd(file.path(here::here(), "02_Results/01b_ddRAD_Bringloe/02_Admixture/") ) 
  
  cmd <- paste("--cv", # to perform cross-validation in the log file 
               bed.file,
               k, # the number of K
               #"-B999",
               "-j8"#
  )
  
  A <- system2("admixture", cmd, stdout = T, stderr = T) 
  
  cat(file = paste0("Beluga_lcWGS.k",k, ".log"),
      "\n", cmd, "\n",
      A, # what to put in my file
      append= F, sep = "\n")
  
  setwd(here::here())
  
}

# Cross-validation results:

CV.res <- data.frame(k = 1:10,
                     CV = NA,
                     stringsAsFactors = F)


for(i in 1:nrow(CV.res)){
  # Which k
  k <- CV.res[i, "k"]
  
  # Extract from the log file
  temp <- readLines(file.path(here::here(),"02_Results/01b_ddRAD_Bringloe/02_Admixture/", paste0("Beluga_lcWGS",k, ".log")))
  CV.temp <- temp %>% str_subset("CV error")
  CV <- sapply(str_split(CV.temp, ":"), `[`, 2) %>% str_remove_all(" ")
  
  # Add to my data.frame
  CV.res[i, "CV"] <- CV
  
}

CV.res$CV <- as.numeric(as.character(CV.res$CV))

CV.res %>% arrange(CV)

plot(CV.res$CV)

gg.CV <- CV.res %>% mutate(color = ifelse(k == 1, "red", "black")) %>% 
  ggplot(aes(x = factor(k), y = CV)) + 
  geom_point(size = 2, aes(col = color)) +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "K", y = "Cross-validation error") +
  theme_bw() +
  theme(legend.position = "none")
gg.CV

ggsave(filename = file.path(here::here(), "02_Results/01b_ddRAD_Bringloe/02_Admixture/", "Admixture.test.CV.png"), 
       plot = gg.CV,
       height = 3.5, width = 4, units = "in")   


k <- 10
Q.k2.res <-  read.table(file.path(here::here(), "02_Results/01b_ddRAD_Bringloe/02_Admixture/", paste0("Beluga_lcWGS.",2,".Q")))
Q.k3.res <-  read.table(file.path(here::here(), "02_Results/01b_ddRAD_Bringloe/02_Admixture/", paste0("Beluga_lcWGS.",3,".Q")))
Q.k4.res <-  read.table(file.path(here::here(), "02_Results/01b_ddRAD_Bringloe/02_Admixture/", paste0("Beluga_lcWGS.",4,".Q")))
Q.k5.res <-  read.table(file.path(here::here(), "02_Results/01b_ddRAD_Bringloe/02_Admixture/", paste0("Beluga_lcWGS.",5,".Q")))
Q.k6.res <-  read.table(file.path(here::here(), "02_Results/01b_ddRAD_Bringloe/02_Admixture/", paste0("Beluga_lcWGS.",6,".Q")))

Q.res <- bind_rows(cbind(fam$V1, Q.k6.res, K = 6),
                   cbind(fam$V1, Q.k5.res, K = 5),
                   cbind(fam$V1, Q.k4.res, K = 4),
                   cbind(fam$V1, Q.k3.res, K = 3),
                   cbind(fam$V1, Q.k2.res, K = 2))


head(Q.res)

#Q.fas.res <- cbind(fam.fas$V1, Q.fas.res)

names(Q.res) <- c("ID_GQ", paste0("Q", 1:6), "K")

#reorder(ID, Qvalue, FUN = function(x) min(x))
pop.data$Membership <- factor(pop.data$Membership, levels = c("WHB","PAN","RES","LGR","CBS","JAM-BEL","SLE"))
pop.data$Region <- factor(pop.data$Region, levels = c("RES", "REB", "NHB", "NWH", "SWH", "JAM", "SAN", "SEH", "NEH", "NHS", "SHS", "UNG", "FRB", "CBS", "SLE"))

# Plot all levels of k ------------------------------
gg.str.all <- Q.res %>% pivot_longer(cols =  paste0("Q", 1:6), names_to = "Group", values_to = "Q") %>% 
  #  mutate(Group = factor(Group, levels = c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6"))) %>% 
  #  dplyr::filter(ID_GQ %nin% bad.samples) %>% 
  left_join(pop.data) %>% 
  #ggplot(aes(x = reorder(ID_GQ, Longitude_echantillonnage_DD, FUN = function(x) max(as.numeric(x))), y = Q, fill = Group)) + 
  ggplot(aes(x = ID_GQ, y = Q, fill = Group)) + 
  
  geom_col() +
  #facet_grid(. ~Lieu_echantillonnage + Mois_echantillonnage, space = "free", scale = "free")  +
  facet_grid(K ~ Membership, space = "free", scale = "free") +
  scale_fill_manual(values=c("royalblue2", "deepskyblue", "orange", "springgreen4", "red3","lightblue","darkblue"))
labs(y="Membership probability") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        strip.text.y = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "white", colour  = "white"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt") )
gg.str.all


ggsave(file.path(here::here(), "/02_Results/01b_ddRAD_Bringloe/02_Admixture/Admixture_k2_to_k6.png"), plot = gg.str.all,
       width = 12, height = 6, unit = "in")

# Plot k=6 according to defined season ----------------------------
gg.str.k6 <- Q.res %>% pivot_longer(cols =  paste0("Q", 1:6), names_to = "Group", values_to = "Q") %>% 
  mutate(Group = factor(Group, levels = c("Q3", "Q6", "Q5", "Q4", "Q2", "Q1"))) %>% 
  dplyr::filter(K == 6) %>% 
  left_join(pop.data) %>% 
  mutate(Region2, factor(Region2, levels = c("SLE","CSB","FRB","UNG","SHS","NHS","NEH","EHB","BEL","JAM","WHB")),
         Season = ifelse(is.na(Month), "Unknown",
                         ifelse(Month %in% c(7,8), "Summer",
                                ifelse(Month > 3 & Month < 7, "Spring",
                                       ifelse(Month > 8 & Month < 12, "Fall",
                                              "Winter"))))) %>% 
  #ggplot(aes(x = reorder(ID_GQ, Season, FUN = function(x) max(levels(x))), y = Q, fill = Group)) + 
  ggplot(aes(x = ID_GQ, y = Q, fill = Group)) + 
  
  geom_col() +
  #facet_grid(. ~Lieu_echantillonnage + Mois_echantillonnage, space = "free", scale = "free") +
  facet_grid(. ~ Region2 + Season, space = "free", scale = "free") +
  scale_fill_brewer(palette = "Set1") +
  labs(y="Membership probability") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        strip.text.y = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "white", colour  = "white"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt") )
gg.str.k6

ggsave(file.path(here::here(), "/02_Results/01b_ddRAD_Bringloe/02_Admixture/Admixture_seasons_k6.png"), plot = gg.str.k6,
       width = 12, height = 6, unit = "in")


# Plot Fst results ------------------------------------------------------

# Plot density graphs

library(ggridges)
library(viridis)

# Read in pixy output modified with an extra column <Comparison> for pairwise comparison made
fst_distances <- read.table(file.path(here::here(), "02_Results/01b_ddRAD_Bringloe/03_FST", "MOBELS_S_20_00703_Beluga_invariant.variant_16v23_final_fst_modified.txt"), header = TRUE, sep = "\t")
                            
 gg.fst <- ggplot(fst_distances, aes(x= avg_wc_fst, y = Comparison, fill = stat(x))) + 
 scale_fill_viridis_c(name= "wc_Fst", option="D", limits=c(-0.1,0.15), oob = scales::squish) +
 geom_density_ridges_gradient(jittered_points = TRUE, quantile_lines = TRUE, scale = 0.9, alpha = 0.7, vline_size = 1, vline_color = "white",
 point_size = 0.4, point_alpha = 0.1, position = position_raincloud(adjust_vlines = FALSE)) + 
 theme_classic() + xlim(-0.1,0.6)

 gg.fst
 
 ggsave(file.path(here::here(), "/02_Results/01b_ddRAD_Bringloe/03_FST/fst_density.png"), plot = gg.fst,
 width = 6, height = 12, unit = "in")
                            
# Plot heatmaps
                            
# Plot correlation
 
# Plots coefficients of inbreeding -------------------------------------------------------------
                            
# Plot distribution maps of genetic clusters -------------------------------------------------------------
                            
library(terra)
library(ggspatial)
library(tidyterra)
                            
admin <- terra::vect(rnaturalearth::ne_countries(scale = "medium", returnclass	  = "sf")) 
#target_crs <- st_crs("+proj=stere +lat_0=90 +lat_ts=71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")
world <- ne_countries(scale = "medium", returnclass = "sf") #%>%
# sf::st_transform(crs = target_crs)
                    
pop.df <- pop.data %>% 
mutate(Region2, factor(Region2, levels = c("SLE","CSB","FRB","UNG","SHS","NHS","NEH","EHB","BEL","JAM","WHB")),
Region = ifelse(Region == "LON", "JAM", Region)) %>% 
group_by(Region2, Region) %>% summarise(LonM = mean(Lon,  na.rm = TRUE), LatM = mean(Lat,  na.rm = TRUE),
N = n()) 
pop.df
pop.df[which(pop.df$Region == "REB"), c("LonM")] <- c(-86.3) #, 66.6)
pop.df[which(pop.df$Region == "REB"), c("LatM")] <- c(66.6)
pop.df[which(pop.df$Region == "SWH"), c("LonM")] <- c(-92.24) #, 66.6)
pop.df[which(pop.df$Region == "SWH"), c("LatM")] <- c(57.2)
pop.df[which(pop.df$Region == "FRB"), c("LonM")] <- c(-67.54) #, 66.6)
pop.df[which(pop.df$Region == "FRB"), c("LatM")] <- c(62.87)
pop.df[which(pop.df$Region == "RES"), c("LonM")] <- c(-95.07) #, 66.6)
pop.df[which(pop.df$Region == "RES"), c("LatM")] <- c(74.76)
pop.df[which(pop.df$Region == "UNG"), c("LonM")] <- c(-69.52) #, 66.6)
pop.df[which(pop.df$Region == "UNG"), c("LatM")] <- c(59.59)
pop.df
                          
                            
pop <- terra::vect(pop.df,
geom = c("LonM", "LatM"),  crs = "+proj=longlat")
                
plot(crop(admin, ext(pop)))
plot(pop, add = T)
                            
gg.env <- ggplot(pop) + 
  geom_sf(data = admin, fill="gray95", size=0.5, alpha = 1) +
                              
  #geom_sf(data = pop, alpha = .8, aes(col = Region2 , size = N )) + 
  scale_size_continuous(limits=c(0,125),breaks=c(25,50,75,100,125)) + 
  scale_colour_manual(values=c("chocolate3","springgreen4","red3", "plum2","orange","salmon","palevioletred2", "lightblue","orchid3", "deepskyblue", "purple3","royalblue1","grey")) +
  # Map limits
  #scale_shape_manual(values = c(16,15,17), name = "Year") +
  coord_sf(xlim = c(-1750000, 1250000), ylim = c(350000, 3750000), crs = sf::st_crs("EPSG:6622")) +
  # Others
  #facet_grid(. ~ pred.pop) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Beluga lcWGS") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), strip.text = element_text(size=10),
  strip.background = element_rect(fill = "white"),
  legend.position = "bottom")  
gg.env
                            
ggsave(file.path(here::here(), "/02_Results/01b_ddRAD_Bringloe/04_Maps/sample_sizes_lcWGS.png"), plot = gg.env,
 width = 12, height = 6, unit = "in")
                            
# Plot worldmap for base layer depicting beluga distribution (inset of figure 1)
                            
gg.world <- ggplot() + geom_sf(data = world, fill="gray95", linewidth=0.15, alpha = 1) +
   coord_sf(xlim = c(-170, -20), ylim = c(30, 80)) + #scale_fill_brewer(palette = "Set3") +
   coord_sf(crs = "+proj=laea +lat_0=90 +lon_0=-100 +x_0=0 +y_0=0") + 
   xlab("Longitude") + ylab("Latitude") +
  theme_bw(base_size = 20) +
  theme(panel.grid = element_blank(), strip.text = element_text(size=10),
  strip.background = element_rect(fill = "white"))
gg.world
                            
ggsave(file.path(here::here(), "/02_Results/01b_ddRAD_Bringloe/04_Maps/worldmap.png"), plot = gg.world,
   width = 12, height = 6, unit = "in")
                            
# Plot maps with distribution of pie charts for genetic clusters
                            
library(scatterpie)
                            
Q.pie.membership <- read.table(file.path(here::here(), "/02_Results/01b_ddRAD_Bringloe/04_Maps/", "MOBELS-S_20_00703-Beluga-lcWGS-Membership_proportions-10v23.csv"), header = TRUE, sep = "\t")
                                                           
# Pie charts will be stretched due to coordinate system at high latitudes; we worked around this by plotting the pie charts and map seperately (geom_sf and geom_scatterpie functions), then overlaying them in a vector editor (inkscape)
# A more elegant solution to correct stretching must exist, but we did not spend the time to solve
gg.map <- ggplot() + geom_sf(data = admin, fill="gray95", size=0.5, alpha = 1) +
geom_scatterpie(aes(x=LonM, y=LatM, group = Region2, r = 1.5), 
 data = Q.pie.membership, cols = "Group",long_format = T) +
coord_sf(xlim = c(-100, -50), ylim = c(45, 75)) + #scale_fill_brewer(palette = "Set3") +
scale_fill_manual(values=c("darkblue", "royalblue2", "red3", "lightblue", "springgreen4", "orange", "deepskyblue")) + 
facet_grid(~Season) +
xlab("Longitude") + ylab("Latitude") +
theme_bw(base_size = 11) +
theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))
                                                           
ggsave(file.path(here::here(), "/02_Results/01b_ddRAD_Bringloe/04_Maps/cluster_distribution.png"), plot = gg.map,
 width = 12, height = 6, unit = "in")