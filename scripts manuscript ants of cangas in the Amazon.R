rm(list=ls(all=TRUE))

#Choose color:
cor_est <- c("#8B4513","#006400")
#These colors will be used for the stratum (Ground and Vegetation, respectively)


##__1 - Database organization----
dir()
com <- read.table("dataset_ants_canga_2024-12-24d.txt",h=T)
head(com)
names(com)

#Ant community organization
spp <- com[,10:35]
spp_simp <- data.frame(aggregate(spp, list(com$a_amb), sum, na.rm=T))
spp_com <- spp_simp[,2:27] #Community occurrence frequency
spp_com_log <- data.frame(log(spp_com+1)) #Logarithmized community occurrence frequency

#Organizing the categorical variable
nome <- data.frame(spp_simp[,1])
colnames(nome) <- "ponto"
library(tidyr)
var <- nome %>% separate(ponto, c("are","est"), "_")
cat <- data.frame(var, est_num=c(2,1,2,1,2,1,2,1), nome)

#Organizing species richness
range(spp_com)
spp_com_pa <- data.frame(ifelse(spp_com>0,1,0))
cat$rich <- rowSums(spp_com_pa)



##__2 - Richness ----
hist(cat$rich)
shapiro.test(cat$rich) #Data have a normal distribution
#We will use a paired t-test

solo <- cat[cat$est=="solo",5]
veget <- cat[cat$est=="veget",5]
t.test(solo, veget, paired = T, alternative="greater")

#Plot----
data <- data.frame(Richnes=cat$rich, Stratum=rep(c("Ground","Vegetation"),4))

solo_cor <- data[which(data$Stratum==data$Stratum[1]),]
vege_cor <- data[which(data$Stratum==data$Stratum[2]),]
data2 <- data.frame(rbind(solo_cor, vege_cor))
data2$x_cor <- c(1,1,1,1,2,2,2,2)

data3 <- data.frame(solo_cor, vege_cor)
data3$x_s <- c(1,1,1,1)
data3$x_v <- c(2,2,2,2)

library(Rmisc)
tgc <- summarySE(data, measurevar="Richnes", groupvars=c("Stratum"))

library(ggplot2)
ggplot(data = tgc, aes(x=Stratum, y=Richnes, fill=Stratum))+
  geom_col(size=.7, colour=1, alpha=0.5, width=0.7)+
  geom_errorbar(aes(ymin=Richnes-se, ymax=Richnes+se), width=.1, position=position_dodge(.7), size=.9)+
  geom_segment(data = data3, aes(x = x_s, xend = x_v, y = Richnes, yend = Richnes.1), colour = "gray50", size = .7, alpha=.6)+
  geom_point(data = data2, aes(x = x_cor, y = Richnes), colour = 1, size = 4, alpha=.7, shape=16)+
  scale_y_continuous(breaks=c(0,3,6,9,12))+
  scale_fill_manual(values=cor_est)+
  labs(x = "Stratum", y = "Ants richnes (mean±SE)", size = 3)+
  theme_classic(base_size = 15)+
  theme(legend.position = "none",
        axis.ticks = element_line(color = "black", size = .8),
        axis.line.x.bottom = element_line(color = "black", linewidth = .8),
        axis.text.x.bottom = element_text(color = "black", size = 13),
        axis.line.y.left = element_line(color = "black", linewidth = .8),
        axis.text.y.left = element_text(color = "black"))

ggsave("Fig 2 (2024-12-24d).tiff", width = 8, height = 6, dpi = 400, compression = "lzw")




##__3 - Composition----
##_3.1 - CLAM test----
library(vegan)
result_clam <- clamtest(spp_com, cat$est, alpha = 0.05)
summary(result_clam)

#Plot----
dados <- data.frame(result_clam[,c(1,4)],
                    result_clam[,c(2,3)]+1)

#Lines in figure
minval <- summary(result_clam)$minv
lineYvsgen <- minval[[1]]+1 #Line delimit specialist species Y
lineXvsgen <- minval[[2]]+1 #Line delimit specialist species X
Ymin <- minval[[1]][1,2]+1 #Line delimit rare species in Y
Xmin <- minval[[2]][1,1]+1 #Line delimit rare species in X
comp <- data.frame(x=c(Xmin, Xmin), y=c(Ymin, Ymin), z=c(1, 2)) #Complemento de linha

#Name species. Se precisar, remover alguma linha de codigo para não inserir nome
esp_s <- rbind(dados[which(dados$Classes == levels(dados$Classes)[1]), ], #Generalistas
               dados[which(dados$Classes == levels(dados$Classes)[2]), ], #Nivel 1
               dados[which(dados$Classes == levels(dados$Classes)[3]), ]) #Raras

#Colors
#cor_plot <- c("green", "blue") #Objeto com as cores utilizadas nos dois tratamentos
cor_CLAM <- c(1, cor_est, "gray75")


library(ggrepel)
ggplot(dados, aes(x=Total_veget, y=Total_solo, colour=Classes))+
  geom_curve(aes(x=Xmin, y=2, xend=2, yend=Ymin), size=1.3, colour=cor_CLAM[4], curvature = .4)+
  geom_line(data=comp, aes(x=z, y=y), colour=cor_CLAM[4], size=1.3)+ #Complemento
  geom_line(data=comp, aes(x=x, y=z), colour=cor_CLAM[4], size=1.3)+ #Complemento
  geom_line(data=lineYvsgen, aes(x=x, y=y), linetype=2, size=1.3, colour=cor_CLAM[2])+
  geom_line(data=lineXvsgen, aes(x=x, y=y), linetype=2, size=1.3, colour=cor_CLAM[3])+
  geom_point(aes(colour=Classes, shape=Classes), size=3.2)+
  geom_text_repel(data=esp_s, aes(label=gsub("_"," ",Species)), fontface="italic", colour="gray10", segment.color = '#cccccc', size = 3.6)+
  scale_shape_manual(values=c(18,19,17,15), labels = c("Generalist", "Ground","Vegetation","Too Rare"))+
  scale_colour_manual(values = cor_CLAM, labels = c("Generalist", "Ground","Vegetation","Too Rare"))+
  scale_x_continuous(trans = 'log10', breaks = c(1,3,9,30,90)) +
  scale_y_continuous(trans = 'log10', breaks = c(1,3,9,30,90))+
  labs(y="Soil (abundance + 1)", x="Vegetation (abundance + 1)")+
  theme_classic(base_size = 15)+
  theme(legend.position=c(0.08, 0.94),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.background = element_rect(color = "black", size = .5, linetype = "solid"),
        legend.key.size = unit(0, 'lines'),
        axis.title = element_text(size = 12),
        axis.line = element_line(color = "black", linewidth=.8),
        axis.ticks = element_line(color = "black", size = .8),
        axis.text = element_text(color = "black", size = 10),
        )

ggsave("Fig 4 (2024-12-24d).tiff", width = 8, height = 6, dpi = 400, compression = "lzw")



##_3.2 - PERMANOVA----
PERMANOVA <- adonis2(spp_com_log ~ cat$est, permutations=9999, method = "bray")
PERMANOVA


##_3.3 - SIMPER----
SIMPER_list <- simper(spp_com_log,cat$est)
SIMPER_tab <- summary(SIMPER_list)
SIMPER_tab2 <- data.frame(SIMPER_tab$solo_veget)

SIMPER_acum <- SIMPER_tab$solo_veget$cumsum
SIMPER_remov <- as.vector(c(0,SIMPER_acum[-26]))
SIMPER_contr <-(SIMPER_acum-SIMPER_remov)*100
SIMPER_final <- data.frame(SIMPER_tab2[,6]*100,SIMPER_contr)
rownames(SIMPER_final) <- rownames(SIMPER_tab2)
colnames(SIMPER_final) <- c("% Aumulada", "% Especie")
SIMPER_final$names <- rownames(SIMPER_tab2)
SIMPER_final #It will be important to plot in NMDS, use 70%, 10 lines...
SIMPER_final2 <- SIMPER_final[1:10,]




##_3.4 - NMDS----
comp_bray <- vegdist(spp_com_log, method="bray")
comp_bray

nmds <- metaMDS(spp_com_log, k=2, maxit=999)
stressplot(nmds)
nmds

#Plot----
summary(nmds)
nmds_data <- scores(nmds)
data.scores <- as.data.frame(nmds_data$sites)
data.scores2 <- data.frame(data.scores,cat)
summary(data.scores2)

data.scores3 <- fortify(data.scores2)

#Creating the convex hull
grp.a <- data.scores[data.scores2$est == "solo", ][chull(data.scores2[data.scores2$est == "solo", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- data.scores[data.scores2$est == "veget", ][chull(data.scores2[data.scores2$est == "veget", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
hull.data <- rbind(grp.a, grp.b)  #combine grp.a and grp.b
hull.data$grp <- as.vector(c(rep("Solo",dim(grp.a)[1]), rep("Veget",dim(grp.b)[1])))
dim(grp.a)

#Species scores
species.scores <- as.data.frame(nmds_data$species)
species.scores$species <- rownames(species.scores)
#Selecting species using SIMPER
species.scores2 <- species.scores[as.vector(SIMPER_final2$names),]


#est_cate <- as.factor(rep(c("Solo", "Vegeta??o"), 4))


##Plot
ggplot(data = data.scores3, aes(x=NMDS1, y=NMDS2)) +
  geom_vline(xintercept = 0, linetype= 2, colour = "gray50")+
  geom_hline(yintercept = 0, linetype= 2, colour = "gray50")+
  stat_ellipse(geom = "polygon", aes(fill=est), alpha=0.3)+
  geom_point(aes(shape=est, colour=est), size = 4)+
  geom_point(data = species.scores, aes(x = NMDS1, y = NMDS2), size = 2, shape=3, colour = 1, alpha=.7, position = position_jitter(width = 0.01, height = 0.01))+
  geom_text_repel(data = species.scores2, aes(x=NMDS1, y=NMDS2, label=gsub("_"," ", species)), fontface="italic", colour=1, segment.color = '#cccccc', size = 3, alpha=.7)+
  scale_color_manual("Stratum", values = cor_est, labels = c("Ground","Vegetation"))+
  scale_fill_manual("Stratum", values = cor_est, labels = c("Ground","Vegetation"))+
  scale_shape_manual("Stratum", values = c(19,17), labels = c("Ground","Vegetation"))+
  geom_text(aes(x = 0.8, y = -1.4, label = "Stress = 0.042"), size = 3.5)+
  #xlim(-1.2, 1.2)+
  #ylim(-1.4, 1.4)+
  #theme_bw() +
  #theme_light()+
  theme_classic(base_size = 15)+
  theme(legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.background = element_rect(color = "black", size = .5, linetype = "solid"),
        legend.key.size = unit(0.1, 'lines'),
        panel.background = element_rect(colour = "black", size=.5),
        panel.border = element_rect(color = "black", size=1, fill = NA),
        axis.line = element_line(color = "black", size=1),
        axis.title = element_text(size = 12),
        axis.ticks = element_line(color = "black", size = .8),
        axis.text = element_text(color = "black", size = 10))

ggsave("Fig 3 (2024-12-24d).tiff", width = 8, height = 6, dpi = 400, compression = "lzw")




##_References----
version
citation()
citation("vegan")
citation("ggplot2")

dev.off()