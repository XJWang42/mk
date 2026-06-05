library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ggrepel)
library(grid)
library(scales)

het_mk<-read.table("mk_low_het.txt", h=T, sep="\t")

sim_targets <- data.frame(group = c("early", "middle", "recent"),mean_GL = c(1.8199, 1.3525, 1.2100))
df <- merge(het_mk, sim_targets, by = "group", all.x = TRUE)
fit<-lm(mean_GL ~ high_all_corr + moderate_all_corr - 1, data = df)
coef(fit)
fit<-lm(mean_GL ~ cadd_01_total_corr + cadd_05_total_corr - 1, data = df[which(df$mean_GL>1),])

het_mk$realised_snpeff <- het_mk$high_hom_corr*0.004006967 + het_mk$moderate_hom_corr*0.003411724 
het_mk$masked_snpeff <- het_mk$high_het_corr*0.004006967/2 + het_mk$moderate_het_corr*0.003411724/2
het_mk$total_snpeff <- het_mk$realised_snpeff + het_mk$masked_snpeff

het_mk$realised_cadd <- het_mk$cadd_01_hom_corr*0.037198194 + het_mk$cadd_05_hom_corr*0.002187137 
het_mk$masked_cadd <- het_mk$cadd_01_het_corr*0.037198194/2 + het_mk$cadd_05_het_corr*0.002187137/2
het_mk$total_cadd <- het_mk$realised_cadd + het_mk$masked_cadd

effect_df <- df[which(df$mean_GL>1),]
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
