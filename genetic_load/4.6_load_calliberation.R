library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ggrepel)
library(grid)
library(scales)

het_mk<-read.table("mk_low_het.txt", h=T, sep="\t")

sim_targets <- data.frame(group = c("early", "middle", "recent"),mean_GL = c(0.7881, 0.6708, 0.6154))
df <- merge(het_mk, sim_targets, by = "group", all.x = TRUE)
fit_snpeff<-lm(mean_GL ~ high_all_corr + moderate_all_corr - 1, data = df)
coef(fit_snpeff)
fit_cadd<-lm(mean_GL ~ cadd_01_total_corr + cadd_05_total_corr - 1, data = df)
coef(fit_cadd)

#test for load accumulate
het_mk$realised_snpeff <- het_mk$high_hom_corr*coef(fit_snpeff)[[1]] + het_mk$moderate_hom_corr*coef(fit_snpeff)[[2]] 
het_mk$masked_snpeff <- het_mk$high_het_corr*coef(fit_snpeff)[[1]]/2 + het_mk$moderate_het_corr*coef(fit_snpeff)[[2]]/2
het_mk$total_snpeff <- het_mk$realised_snpeff + het_mk$masked_snpeff

het_mk$realised_cadd <- het_mk$cadd_01_hom_corr*coef(fit_cadd)[[1]] + het_mk$cadd_05_hom_corr*coef(fit_cadd)[[2]] 
het_mk$masked_cadd <- het_mk$cadd_01_het_corr*coef(fit_cadd)[[1]]/2 + het_mk$cadd_05_het_corr*coef(fit_cadd)[[2]]/2
het_mk$total_cadd <- het_mk$realised_cadd + het_mk$masked_cadd

effect_df <- df[which(df$mean_GL>0),]
loo_results <- lapply(1:nrow(effect_df), function(i) {
  temp_df <- effect_df[-i, ]
  fit_i <- lm(mean_GL ~ high_all_corr + moderate_all_corr - 1, data = temp_df)
  
  data.frame(
    s_HIGH = coef(fit_i)[[1]],
    s_MOD  = coef(fit_i)[[2]]
  )
})

loo_results <- bind_rows(loo_results)
summary(loo_results$s_HIGH)
summary(loo_results$s_MOD)

loo_results <- lapply(1:nrow(effect_df), function(i) {
  temp_df <- effect_df[-i, ]
  fit_i <- lm(mean_GL ~ cadd_01_total_corr + cadd_05_total_corr - 1, data = temp_df)
  
  data.frame(
    s_HIGH = coef(fit_i)[[1]],
    s_MOD  = coef(fit_i)[[2]]
  )
})

loo_results <- bind_rows(loo_results)
summary(loo_results$s_HIGH)
summary(loo_results$s_MOD)


mk_compare <- list(c("historical", "early"),c("middle", "early"), c("middle", "recent"), c("early", "recent"))
realised_snpeff <- ggplot(het_mk[which(het_mk$sample!="aMK_SB787"),], aes(x=group, y=realised_snpeff, color=group, fill=group)) +
  geom_boxplot() +
  geom_point() +
  stat_compare_means(comparisons = mk_compare, method = "wilcox.test", paired=FALSE)+
  scale_x_discrete(limit = c('historical','early','middle','recent')) +
  scale_color_manual(values = c("#562A2A","#706040","#305030","#405060"),breaks = c('early','middle','recent','historical')) +
  scale_fill_manual(values = c("#C36969","#C0AC74","#577B41","#8497B0"),breaks = c('early','middle','recent','historical')) +
  theme_classic() +
  theme(legend.position = "none")
masked_snpeff <- ggplot(het_mk[which(het_mk$sample!="aMK_SB787"),], aes(x=group, y=masked_snpeff, color=group, fill=group)) +
  geom_boxplot() +
  geom_point() +
  stat_compare_means(comparisons = mk_compare, method = "wilcox.test", paired=FALSE)+
  scale_x_discrete(limit = c('historical','early','middle','recent')) +
  scale_color_manual(values = c("#562A2A","#706040","#305030","#405060"),breaks = c('early','middle','recent','historical')) +
  scale_fill_manual(values = c("#C36969","#C0AC74","#577B41","#8497B0"),breaks = c('early','middle','recent','historical')) +
  theme_classic() +
  theme(legend.position = "none")
total_snpeff <- ggplot(het_mk[which(het_mk$sample!="aMK_SB787"),], aes(x=group, y=total_snpeff, color=group, fill=group)) +
  geom_boxplot() +
  geom_point() +
  stat_compare_means(comparisons = mk_compare, method = "wilcox.test", paired=FALSE)+
  scale_x_discrete(limit = c('historical','early','middle','recent')) +
  scale_color_manual(values = c("#562A2A","#706040","#305030","#405060"),breaks = c('early','middle','recent','historical')) +
  scale_fill_manual(values = c("#C36969","#C0AC74","#577B41","#8497B0"),breaks = c('early','middle','recent','historical')) +
  theme_classic() +
  theme(legend.position = "none")

realised_cadd <- ggplot(het_mk[which(het_mk$sample!="aMK_SB787"),], aes(x=group, y=realised_cadd, color=group, fill=group)) +
  geom_boxplot() +
  geom_point() +
  stat_compare_means(comparisons = mk_compare, method = "wilcox.test", paired=FALSE)+
  scale_x_discrete(limit = c('historical','early','middle','recent')) +
  scale_color_manual(values = c("#562A2A","#706040","#305030","#405060"),breaks = c('early','middle','recent','historical')) +
  scale_fill_manual(values = c("#C36969","#C0AC74","#577B41","#8497B0"),breaks = c('early','middle','recent','historical')) +
  theme_classic() +
  theme(legend.position = "none")
masked_cadd <- ggplot(het_mk[which(het_mk$sample!="aMK_SB787"),], aes(x=group, y=masked_cadd, color=group, fill=group)) +
  geom_boxplot() +
  geom_point() +
  stat_compare_means(comparisons = mk_compare, method = "wilcox.test", paired=FALSE)+
  scale_x_discrete(limit = c('historical','early','middle','recent')) +
  scale_color_manual(values = c("#562A2A","#706040","#305030","#405060"),breaks = c('early','middle','recent','historical')) +
  scale_fill_manual(values = c("#C36969","#C0AC74","#577B41","#8497B0"),breaks = c('early','middle','recent','historical')) +
  theme_classic() +
  theme(legend.position = "none")
total_cadd <- ggplot(het_mk[which(het_mk$sample!="aMK_SB787"),], aes(x=group, y=total_cadd, color=group, fill=group)) +
  geom_boxplot() +
  geom_point() +
  stat_compare_means(comparisons = mk_compare, method = "wilcox.test", paired=FALSE)+
  scale_x_discrete(limit = c('historical','early','middle','recent')) +
  scale_color_manual(values = c("#562A2A","#706040","#305030","#405060"),breaks = c('early','middle','recent','historical')) +
  scale_fill_manual(values = c("#C36969","#C0AC74","#577B41","#8497B0"),breaks = c('early','middle','recent','historical')) +
  theme_classic() +
  theme(legend.position = "none")

pushViewport(viewport(layout = grid.layout(2,3)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(realised_snpeff, vp = vplayout(1,1))
print(masked_snpeff, vp = vplayout(1,2))
print(total_snpeff, vp = vplayout(1,3))
print(realised_cadd, vp = vplayout(2,1))
print(masked_cadd, vp = vplayout(2,2))
print(total_cadd, vp = vplayout(2,3))
