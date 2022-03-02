setwd("C:/Users/Charlie/Desktop/SCR Arthropod Survey")
ant_full<-read.csv("Riparian Community Final.csv")
Ant_community<-ant_full[,6:23]
ant.info <- ant_full[,1:5]
library(vegan)
sample(1:999, 1)
set.seed(257)
ant.dist<- vegdist(Ant_community, method = "bray")
ant.nmds <- metaMDS(Ant_community, k=2, model = "global", distance = "bray")
ant.nmds$stress
plot(ant.nmds)
stressplot(ant.nmds)

###NMDS plot
ant_full$Habitat.Type<-as.factor(ant_full$Habitat.Type)
summary(ant_full$Habitat.Type)
ant_full$Season<-as.factor(ant_full$Season)
pchvec<-c(15,17)
library(RColorBrewer)
brewer.pal(6, "Set1")
colvec<-c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#FFFF33")
plot(ant.nmds)
with(ant_full,
     points(ant.nmds, cex=1.0,
            col = colvec[Habitat.Type],
            pch = pchvec[Season]))
mtext(text = "Stress = 0.092", side = 3, line = -1.5, adj = 0.95, 
      cex = 1.0)
legend("bottomright", c("Arundo", "R. Channel", "R. Forest", "R. Scrub", "Replanting", "Watercress" ), title = "Microhabitat", 
       col= colvec, pch=15, cex=0.8)
legend("bottomleft", c("Summer", "Winter"), title = "Season", pch=c(15,17))
ordiellipse(ant.nmds, groups = ant_full$Habitat.Type, draw = "polygon", 
            label = TRUE, lty = 1, col = colvec)
orditorp(ant.nmds, display = "species", col="black")

plot(ant.nmds)
with(ant_full,
     points(ant.nmds, cex=1.0,
            col = colvec[Habitat.Type],
            pch = pchvec[Season]))
ordihull(ant.nmds, groups = ant_full$Habitat.Type, draw = "polygon", 
         label = TRUE, lty = 1, col = colvec) ###conne

plot(ant.nmds)
with(ant_full,
     points(ant.nmds, cex=1.0,
            col = colvec[Habitat.Type],
            pch = pchvec[Season]))
ordispider(ant.nmds, groups = ant_full$Habitat.Type, draw = "polygon", 
           label = TRUE, lty = 1, col = colvec)

###dissimilarity analysis by microhabitat
ant.ano_habitat <- anosim(ant.dist, ant.info$Habitat.Type)
summary(ant.ano_habitat)
plot(ant.ano_habitat)

###permanova

ant.dist<- vegdist(Ant_community, method = "bray")
sample(1:999, 1)
set.seed(662)
PERMANOVA<-adonis(ant.dist~ Habitat.Type, data=ant_full,
                  permutations = 10000)
PERMANOVA


###betadisper to identify groups
mod.habitat<-betadisper(ant.dist, ant_full$Habitat.Type)
anova(mod.habitat) ###significant, so PERMANOVA not ideal test
mod.habitat.HSD<-TukeyHSD(mod.habitat)
plot(mod.habitat.HSD)
mod.habitat.HSD

#environmental fit with envfit
Total.enviro<-cbind(ant_full[,3:5], ant_full[,25:81])
sample(1:999, 1)
set.seed(275)

ef<-envfit(ant.nmds, Total.enviro, permu=999, na.rm= TRUE )
ef
plot(ant.nmds)
ordiellipse(ant.nmds, groups = ant_full$Habitat.Type, draw = "polygon", 
            label = TRUE, lty = 1)
plot(ef, p.max = 0.05)


##ef just plants
plants<-ant_full[,34:81]
set.seed(275)
ef2<-envfit(ant.nmds, plants, permu=999, na.rm= TRUE )
ef2
plot(ant.nmds, display = "sites")
ordiellipse(ant.nmds, groups = ant_full$Habitat.Type, draw = "polygon", 
            label = TRUE, lty = 1)
plot(ef2, p.max = 0.05)




###db-rda (distance based redundancy analysis)
full_dbRDA<-capscale(Ant_community ~Mulch + Litter + Soil + soil.moisture +
                     percent.clay + percent.sand + canopy.cover + Erigeron.canadensis +
                       Curcubita.foetidissima + Arundo.donax + Rubus.ursinus + Salix.lasiolepis +
                       Baccharis.pilularis + Baccharis.salicifolia, Total.enviro, dist="bray", na.action = na.exclude)
summary(full_dbRDA)
ordiplot(full_dbRDA)

anova(full_dbRDA, by="terms", peru=200)
anova(full_dbRDA, by="axis")


###Variance partitioning to differentiate plants, other env factors, and locality/season?
ant.dist<- vegdist(Ant_community, method = "bray")

varpar_table_no_NA<-na.omit(ant_full)
sigplants<-cbind(varpar_table_no_NA$Arundo.donax, varpar_table_no_NA$Baccharis.pilularis, varpar_table_no_NA$Bromus.madritensis.ssp..rubens,
                 varpar_table_no_NA$Curcubita.foetidissima, varpar_table_no_NA$Heliotropium.curassavicum, varpar_table_no_NA$Hirschfeldia.incana,
                 varpar_table_no_NA$Malva.parviflora, varpar_table_no_NA$Populous.trichocarpa, varpar_table_no_NA$Rubus.ursinus,
                 varpar_table_no_NA$Salix.lasiolepis, varpar_table_no_NA$Salsola.tragus, varpar_table_no_NA$Urtica.dioica)
var_dist<-vegdist(varpar_table_no_NA[,6:23], method = "bray")
v.part<-varpart(var_dist, sigplants, varpar_table_no_NA[,25:32])
plot(v.part)
plot (v.part, digits = 2, Xnames = c('Sig. Plants', 'Other Env.'), bg = c('forestgreen', 'blue3'))

###Simper by microhabitat (dis)similarity comparison
sim_ant_habitat<-simper(Ant_community, group = ant_full$Habitat.Type)
summary(sim_ant_habitat)

###indicator species analysis
library(indicspecies)
indval = multipatt(Ant_community, ant_full$Habitat.Type, control = how(nperm=999))
summary(indval, indvalcomp=TRUE)

###chao estimators
specpool(Ant_community, pool = ant_full$Habitat.Type, smallsample = TRUE)

####ggplot NMDS plot with shape by coverage class
score.ant <-scores(ant.nmds) ###takes the NMDS values for plotting

site.scrs <- as.data.frame(scores(ant.nmds)) #save NMDS results into dataframe
site.scrs <- cbind(site.scrs, Habitat = ant.info$Habitat.Type)
site.scrs <- cbind(site.scrs, Season = ant.info$Season)
site.scrs <- cbind(site.scrs, Coverage = ant.info$Main.Cover)

NMDS_Cov = data.frame(NMDS1 = site.scrs$NMDS1, 
                      NMDS2 = site.scrs$NMDS2, group=as.factor(site.scrs$Coverage))
plot(ant.nmds$points)
ord_cov<-ordiellipse(ant.nmds, site.scrs$Coverage, display = "sites", 
                     kind = "se", conf = 0.95, label = T)

df_ell_cov <- data.frame() 
for(g in levels(NMDS_Cov$group)){
  if(is.null(ord_cov[[g]])) next
  df_ell_cov <- rbind(df_ell_cov, cbind(as.data.frame(with(NMDS_Cov[NMDS_Cov$group==g,],
                                                           vegan:::veganCovEllipse(
                                                             ord_cov[[g]]$cov,ord_cov[[g]]$center,ord_cov[[g]]$scale)))
                                        ,group=g))
}
en_coord_cont = as.data.frame(scores(ef, "vectors", p.max = 0.05)) * ordiArrowMul(ef)
NMDS_Cov = data.frame(NMDS1 = site.scrs$NMDS1, 
                      NMDS2 = site.scrs$NMDS2, group=as.factor(site.scrs$Coverage))



###getting values for ellipses
NMDS = data.frame(NMDS1 = site.scrs$NMDS1, 
                  NMDS2 = site.scrs$NMDS2, group=as.factor(site.scrs$Habitat))
plot(ant.nmds$points)
ord<-ordiellipse(ant.nmds, site.scrs$Habitat, display = "sites", 
                 kind = "se", conf = 0.95, label = T)
df_ell <- data.frame() 
for(g in levels(NMDS$group)){
  if(is.null(ord[[g]])) next
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                                                   vegan:::veganCovEllipse(
                                                     ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
}

###NMDS Plot with ordiellipses and environmental vectors
library(ggplot2) ###graphing package, will probably need to install just like tidyr. code is install.packages(ggplot2)
ggplot(site.scrs, aes(x= NMDS1, y=NMDS2))+ ###gives a plot without the environmental vectors
  geom_point(aes(color=Habitat, shape=Season),size=4) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  geom_path (data=df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1) ###+
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont)) 

###species accumulation curves
library(BiodiversityR)
ant_full$Habitat.Type<-as.factor(ant_full$Habitat.Type)
all_microhab_spcurve<-accumcomp(Ant_community, y=ant_full, factor='Habitat.Type', 
                                 method='exact', gamma = "chao", xlim=c(1,50), ylim = c(1,18) )
###use cursor to place legend after running
sp_curve_full<-specaccum(Ant_community, method = "exact", permutations = 100,
          conditioned = TRUE, gamma = "chao")
accumplot(sp_curve_full, addit = T, col="black")

###testing biogeography and ecological diversity null hypotheses
library(rareNMtests)
sample(1:999, 1)
set.seed(42)
cbecoq0<-EcoTest.sample(Ant_community, by=ant_full$Habitat.Type, MARGIN = 1, niter=200, method = "coverage", q=0)
windows();plot(cbecoq0)###sig diff, combined with non-significant bio null test means it is Cayuela et al 2015 
### scenario 2 of figure 1 (Differing community composition) for overall sample set. Now investigate each pair below

###biogeographic null model, takes longer to run
sample(1:999, 1)
set.seed(384)
bioq1<-BiogTest.sample(Ant_community, by=ant_full$Habitat.Type, MARGIN = 1, niter=50, method = "coverage", q=1)
windows();plot(bioq1)
###nmds with season
plot(ant.nmds)
with(ant_full,
     points(ant.nmds, cex=1.0,
            col = colvec[Habitat.Type],
            pch = pchvec[Season]))
mtext(text = "Stress = 0.092", side = 3, line = -1.5, adj = 0.95, 
      cex = 1.0)
legend("bottomright", c("Arundo", "R. Channel", "R. Forest", "R. Scrub", "Replanting", "Watercress" ), title = "Microhabitat", 
       col= colvec, pch=15, cex=0.8)
legend("bottomleft", c("Summer", "Winter"), title = "Season", pch=c(15,17))
ordiellipse(ant.nmds, groups = ant_full$Season, draw = "polygon", 
            label = TRUE, lty = 1,)