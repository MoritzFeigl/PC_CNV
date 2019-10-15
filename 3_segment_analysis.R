#
# Analysis of segments
#

library(ggplot2)
library(tidyr)
library(dplyr)
library(feather)
library(tidyverse)
library(MASS)
library(gridExtra)
segments_df <- read_feather("Data/all_probands_segments_df.feather")

# subset of interesting rows
relevant_seg <- segments_df[segments_df$Length.Mbp > 0.001 &
                              segments_df$Chromosome != "chrY" &
                              segments_df$Chromosome != "chrX", ]
relevant_seg$chr_n <- as.integer(substr(relevant_seg$Chromosome, 4, 5))

# the patients with gleason 10 which have no gleson first/second data
relevant_seg[c(grep("797_", relevant_seg$Proband),
               grep("811_", relevant_seg$Proband),
               grep("863_", relevant_seg$Proband),
               grep("901_", relevant_seg$Proband),
               grep("1210_", relevant_seg$Proband),
               grep("3948_", relevant_seg$Proband),
               grep("10280_", relevant_seg$Proband),
               grep("10605_", relevant_seg$Proband),
               grep("10722_", relevant_seg$Proband),
               grep("10875_", relevant_seg$Proband),
               grep("11067_", relevant_seg$Proband)), "cc"] <- "aggressive"

relevant_seg$Group <- relevant_seg$cc
relevant_seg$Group <- factor(relevant_seg$Group, levels = c("control", "not-aggressive", "aggressive"))
names(relevant_seg)[15] <- "Group"

# Removed because of too many gains/losses: "11384"
length(grep("11384", relevant_seg$Proband))
removed_sample <- c(grep("11384", relevant_seg$Proband))
relevant_seg <- relevant_seg[-removed_sample, ]

# How many CN alterations are there in total?
# Number of segments per Group and status
table(relevant_seg$Group)
length(unique(relevant_seg$Proband))

sum(relevant_seg$status[relevant_seg$Group == "aggressive"] == "loss")
sum(relevant_seg$status[relevant_seg$Group == "aggressive"] == "gain")
sum(relevant_seg$status[relevant_seg$Group == "aggressive"] == "default")
sum(relevant_seg$status[relevant_seg$Group == "not-aggressive"] == "loss")
sum(relevant_seg$status[relevant_seg$Group == "not-aggressive"] == "gain")
sum(relevant_seg$status[relevant_seg$Group == "control"] == "loss")
sum(relevant_seg$status[relevant_seg$Group == "control"] == "gain")

# segment lenghts
boxplot(Length.Mbp ~ Group , 
        relevant_seg[relevant_seg$status == "gain", ])
boxplot(Length.Mbp ~ Group, 
        relevant_seg[relevant_seg$status == "loss", ])

seg_l <- aggregate(Length.Mbp ~ Group + status, relevant_seg, summary)
gains_sl <- seg_l[seg_l$status != "default", ][1:3, ]
losss_sl <- seg_l[seg_l$status != "default", ][4:6, ]

# relevant seg part for merging
relevant_seg_merg <- relevant_seg[, c("Proband", "Group")]
relevant_seg_merg <- relevant_seg_merg[!duplicated(relevant_seg_merg), ]

# Aggregation on proband level
gains_p <- aggregate(status ~ Proband + Chromosome, 
                     data = relevant_seg, 
                     FUN = function(x) sum(x == "gain"), drop = FALSE)
gains_p <- merge(gains_p, relevant_seg_merg, by = "Proband", all.x = TRUE)


loss_p <- aggregate(status ~ Proband + Chromosome, 
                    data = relevant_seg, 
                    FUN = function(x) sum(x == "loss"), drop = FALSE)
loss_p <- merge(loss_p, relevant_seg_merg, by = "Proband", all.x = TRUE)

cna_p <- aggregate(status ~ Proband + Chromosome, 
                   data = relevant_seg, 
                   FUN = function(x) sum(x == "loss" | x == "gain"), drop = FALSE)
cna_p <- merge(cna_p, relevant_seg_merg, by = "Proband", all.x = TRUE)

# Distribution values of segmentes for data on proband level
summary_for_type <- "loss"
type_dist <- aggregate(status ~ Proband, 
                       data = relevant_seg, 
                       FUN = function(x) sum(x == summary_for_type), drop = FALSE)
type_dist <- merge(type_dist, relevant_seg_merg, by = "Proband", all.x = TRUE)
type_dist$Group <- as.character(type_dist$Group)
summary(type_dist[type_dist$Group == "aggressive", "status"])
summary(type_dist[type_dist$Group == "not-aggressive", "status"])
summary(type_dist[type_dist$Group == "control", "status"])

# Data for Chromosomes
gains_ch <- aggregate(status ~ Group + Chromosome, 
                      data = gains_p, 
                      FUN = mean)
gains_ch$chr_n <- as.integer(substr(gains_ch$Chromosome, 4, 5))

loss_ch <- aggregate(status ~ Group + Chromosome, 
                     data = loss_p, 
                     FUN = mean)
loss_ch$chr_n <- as.integer(substr(loss_ch$Chromosome, 4, 5))

cna_ch <- aggregate(status ~ Group + Chromosome, 
                    data = cna_p, 
                    FUN = mean)
cna_ch$chr_n <- as.integer(substr(cna_ch$Chromosome, 4, 5))

# Aggregation on Chromosome level
gains_sc <- aggregate(status ~ Proband, 
                      data = gains_p, 
                      FUN = sum, drop = FALSE)
gains_p_merge <- gains_p[, c("Proband", "Group")]
gains_p_merge <- gains_p_merge[!duplicated(gains_p_merge), ]

gains_sc <- merge(gains_sc, gains_p_merge, by = "Proband", all.x = TRUE, all.y = FALSE)

loss_sc <- aggregate(status ~ Proband, 
                     data = loss_p, 
                     FUN = sum, drop = FALSE)
loss_p_merge <- loss_p[, c("Proband", "Group")]
loss_p_merge <- loss_p_merge[!duplicated(loss_p_merge), ]
loss_sc <- merge(loss_sc, loss_p_merge, by = "Proband", all.x = TRUE, all.y = FALSE)


cna_sc <- aggregate(status ~ Proband, 
                    data = cna_p, 
                    FUN = sum, drop = FALSE)
cna_p_merge <- cna_p[, c("Proband", "Group")]
cna_p_merge <- cna_p_merge[!duplicated(cna_p_merge), ]
cna_sc <- merge(cna_sc, cna_p_merge, by = "Proband", all.x = TRUE, all.y = FALSE)

# Check distribution over samples
boxplot(status~Group, loss_sc, ylim = c(0, 10))
which.max(gains_sc$status)
which.max(loss_sc$status)

# Data for Groups
gains_gr <- aggregate(status ~ Group, 
                      data = gains_sc, 
                      FUN = mean)

loss_gr <- aggregate(status ~ Group, 
                     data = loss_sc, 
                     FUN = mean)
cna_gr <- aggregate(status ~ Group, 
                    data = cna_sc, 
                    FUN = mean)

# Plots for Chromosomes
levels(gains_ch$Group) <- c("control", "non-aggressive", "aggressive")
levels(loss_ch$Group) <- c("control", "non-aggressive", "aggressive")
levels(cna_ch$Group) <- c("control", "non-aggressive", "aggressive")


p1 <- ggplot(gains_ch, aes(x = factor(chr_n), y = status, fill = Group)) + 
  geom_bar( position = "dodge", stat = "identity") +
  ggtitle("CN Gain") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6, size = 12, 
                                   margin = margin(5,0,10,0)),
        axis.text.y = element_text(size = 12, margin = margin(0,5,0,10)),
        axis.title = element_text(size = 10),
        axis.ticks.length = unit(8, "pt"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = "none") +
  xlab ("Chromosome") + ylab("mean number of segments with CN gain") +
  scale_fill_manual(values=c("#00BA38", "#619CFF", "#F8766D"))

p2 <- ggplot(loss_ch, aes(x = factor(chr_n), y = status, fill = Group)) + 
  geom_bar( position = "dodge", stat = "identity") +
  ggtitle("CN Loss") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6, size = 12, 
                                   margin = margin(5,0,10,0)),
        axis.text.y = element_text(size = 12, margin = margin(0,5,0,10)),
        axis.title = element_text(size = 10),
        axis.ticks.length = unit(8, "pt"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = "none") +
  xlab ("Chromosome") + ylab("mean number of segments with CN loss")+
  scale_fill_manual(values=c("#00BA38", "#619CFF", "#F8766D"))

p3 <- ggplot(cna_ch, aes(x = factor(chr_n), y = status, fill = Group)) + 
  geom_bar( position = "dodge", stat = "identity") +
  ggtitle("CN Alteration") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6, size = 12, 
                                   margin = margin(5,0,10,0)),
        axis.text.y = element_text(size = 12, margin = margin(0,5,0,10)),
        axis.title = element_text(size = 10),
        axis.ticks.length = unit(8, "pt"),
        legend.text = element_text(size=16, 
                                   margin = margin(t = 10, b = 10, l = 10, unit = "pt")),
        legend.title = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.key.width = unit(20, unit = "pt"),
        legend.key.height = unit(20, unit = "pt")
  ) +
  xlab ("Chromosome") + ylab("mean number of segments with CN alteration")+
  scale_fill_manual(values=c("#00BA38", "#619CFF", "#F8766D"))

legend <- cowplot::get_legend(p3)

p3 <- ggplot(cna_ch, aes(x = factor(chr_n), y = status, fill = Group)) + 
  geom_bar( position = "dodge", stat = "identity") +
  ggtitle("CN Alteration") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6, size = 12, 
                                   margin = margin(5,0,10,0)),
        axis.text.y = element_text(size = 12, margin = margin(0,5,0,10)),
        axis.title = element_text(size = 10),
        axis.ticks.length = unit(8, "pt"),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = "none") +
  xlab ("Chromosome") + ylab("mean number of segments with CN alterationloss")+
  scale_fill_manual(values=c("#00BA38", "#619CFF", "#F8766D"))
grid.arrange(p1, p2, p3, legend, nrow = 2)

levels(gain_data$Group) <- c("control", "non-aggressive", "aggressive")
levels(loss_data$Group) <- c("control", "non-aggressive", "aggressive")
levels(cna_data$Group) <- c("control", "non-aggressive", "aggressive")

# Plots for Groups
p1 <- ggplot(gain_data, aes(x = Group, y = status, fill = Group)) + 
  geom_bar( position = "dodge", stat = "identity") +
  ggtitle("CN Gain") + ylim(0, 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, 
                                   margin = margin(5,0,10,0)),
        axis.text.y = element_text(size = 12, margin = margin(0,5,0,10)),
        axis.title = element_text(size = 10),
        axis.ticks.length = unit(8, "pt"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = "none") +
  xlab ("") + ylab("mean number of segments with CN gain") +
  scale_fill_manual(values=c("#00BA38", "#619CFF", "#F8766D")) 

p2 <- ggplot(loss_data, aes(x = Group, y = status, fill = Group)) + 
  geom_bar( position = "dodge", stat = "identity") +
  ggtitle("CN Loss") + ylim(0, 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, 
                                   margin = margin(5,0,10,0)),
        axis.text.y = element_text(size = 12, margin = margin(0,5,0,10)),
        axis.title = element_text(size = 10),
        axis.ticks.length = unit(8, "pt"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = "none") +
  xlab ("") + ylab("mean number of segments with CN loss") +
  scale_fill_manual(values=c("#00BA38", "#619CFF", "#F8766D")) 
# geom_errorbar(aes(ymin = LL, ymax = UL), width=.2,
#               position=position_dodge(.9))

p3 <- ggplot(cna_data, aes(x = Group, y = status, fill = Group)) + 
  geom_bar( position = "dodge", stat = "identity") +
  ggtitle("CN Alteration") + ylim(0, 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, 
                                   margin = margin(5,0,10,0)),
        axis.text.y = element_text(size = 12, margin = margin(0,5,0,10)),
        axis.title = element_text(size = 10),
        axis.ticks.length = unit(8, "pt"),
        legend.text = element_text(size=16, 
                                   margin = margin(t = 10, b = 10, l = 10, unit = "pt")),
        legend.title = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.key.width = unit(20, unit = "pt"),
        legend.key.height = unit(20, unit = "pt")) +
  xlab ("") + ylab("mean number of segments with CN alteration") +
  scale_fill_manual(values=c("#00BA38", "#619CFF", "#F8766D")) 

legend <- cowplot::get_legend(p3)

p3 <- ggplot(cna_data, aes(x = Group, y = status, fill = Group)) + 
  geom_bar( position = "dodge", stat = "identity") +
  ggtitle("CN Alteration") + ylim(0, 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, 
                                   margin = margin(5,0,10,0)),
        axis.text.y = element_text(size = 12, margin = margin(0,5,0,10)),
        axis.title = element_text(size = 10),
        axis.ticks.length = unit(8, "pt"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = "none") +
  xlab ("") + ylab("mean number of segments with CN alteration") +
  scale_fill_manual(values=c("#00BA38", "#619CFF", "#F8766D")) 


grid.arrange(p1, p2, p3, nrow = 2)

# ROC FUNCTION 
ROC_curve <- function (logistic.model) {
  # make CV model
  ctrl <- caret::trainControl(method = "repeatedcv", number = 10, repeats = 5, 
                       savePredictions = TRUE, classProbs = TRUE)
  cvmod2 <- caret::train(aggressive ~ .,  data = logistic.model$data, 
                  method="glm", family="binomial",
                  trControl = ctrl, tuneLength = 5)
  # mean CV prediction
  mean_preds <- aggregate(aggressive ~ rowIndex, cvmod2$pred, mean)
  mean_preds <- cbind(mean_preds, "0" = NA, "1" = NA)
  for(i in mean_preds$rowIndex){
    if(logistic.model$data$aggressive[i] == "aggressive"){
      mean_preds$'0'[i] <- 0
      mean_preds$'1'[i] <- 1
    } else {
      mean_preds$'0'[i] <- 1
      mean_preds$'1'[i] <- 0
    }
  }
  Spec <- Sens <- NULL
  firsttable1 <- mean_preds[, -1]
  colnames(firsttable1)[1] <- "predicted.prob"
  firsttable <- firsttable1[, 2:3]
  rownames(firsttable) <- firsttable1$predicted.prob
  secondtable <- firsttable
  for (i in 1:length(secondtable[, 1])) {
    secondtable[i, 1] <- (sum(firsttable[, 1]) - sum(firsttable[(1:i), 
                                                                1]))/sum(firsttable[, 1])
    secondtable[i, 2] <- (sum(firsttable[, 2]) - sum(firsttable[(1:i), 
                                                                2]))/sum(firsttable[, 2])
  }
  secondtable <- as.data.frame(rbind((c(1, 1)), secondtable))
  colnames(secondtable) <- c("Spec", "Sens")
  rownames(secondtable)[1] <- "0"
  auc <- 0
  for (i in 1:(nrow(secondtable) - 1)) {
    auc <- auc + (secondtable[i, 1] - 
                    secondtable[(i + 1), 1]) * 0.5 * 
      (secondtable[i, 2] + secondtable[(i + 1), 2])
  }
  p <- ggplot( 
    data = secondtable[nrow(secondtable):1, ], aes(Spec, Sens)) +
    labs(
      x = "1-Specificity",
      y = "Sensitivity",
      title ="ROC curve",
      subtitle = "estimated with 5 times repeated 10-fold CV"
    ) + xlim(0, 1) + ylim(0, 1) + geom_step() +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20),
          plot.title =  element_text(size = 25), plot.subtitle = element_text(size = 15))
  
  p <- p + geom_segment(x = 0, y = 0, xend = 1, yend = 1, 
                        colour = "black", size = 1)
  p <- p + annotate("text", x = 0.9, y = 0.05, label = paste("AUC = ", 
                                                             formatC(auc, digits = 4)),
                    size = 10)
  p
}


# Logistic regression models -------------------------------------------------------------
model_group <- "cna" #choose either "loss", "gains" or "cna"

# Choose the data for models
type <- get(paste0(model_group, "_p")) # choose a group: gains, loss, cna
levels(type$Group) <- c("control", "non-aggressive", "aggressive")
type$aggressive <- "non.aggressive"
type$aggressive[type$Group == "aggressive"] <- "aggressive"

# aggregate all information on Proband level while keeping the chrom level infos
type_agg <- split(type, type$Proband, drop = TRUE)
cnv_aggregater <- function(x){
  chrs <- x[, c("Chromosome", "status")]
  chrs$Chromosome <- as.integer(substring(chrs$Chromosome, 4))
  df <- data.frame(Proband = x$Proband[1],
                   Group = x$Group[1],
                   aggressive = x$aggressive[1])
  for(i in 1:22){
    df <- cbind(df, chrs$status[chrs$Chromosome == i])
  }
  colnames(df)[-c(1:3)] <- paste0("Chr", 1:22)
  df$Chr_sum <- sum(df[,-c(1:3)])
  
  #df$Chr_max <- which.max(chrs$status)
  return(df)
}
type_agg <- lapply(type_agg, cnv_aggregater)
type_agg <- do.call(rbind, type_agg)
# For no control test
#type_agg <- type_agg[!(type_agg$Group %in% "control"), ]

# 0. Check correlations between chromosomes ----------------------------------------------
M <- cor(type_agg[, 4:25], method = "spearman")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
p.mat <- psych::corr.test(data.frame(type_agg[, 4:25]), adjust="none")$p
mycol <- ifelse(c(p.mat < 0.05), "black", "grey")

jpeg(filename = paste0("D:/Dropbox/CNV Literature/MANUSCRIPT/Plots/", 
                       model_group, "_correlation.jpeg"), width = 1000, height = 1000, 
     units = "px", quality = 100)
corrplot::corrplot(M, method="color", #col = col(200),  
                   type="full", order= "original", 
                   #addCoef.col = "black", 
                   tl.col="black", tl.srt = 90, addCoef.col = mycol,
                   p.mat = p.mat, insig = "blank",
                   diag=TRUE, tl.cex = 1.5, number.cex = 1.2, cl.cex = 1.8
)
dev.off()

# 1. Model: no covariates ------------------------------------------------------
model_df <- type_agg[, -c(1, 2, 26)]
model_df$aggressive <- factor(model_df$aggressive)
#model_df[, 2:23] <- apply(model_df[, 2:23], 2, log1p) 
log_model <- glm(aggressive ~ .,
                 family = binomial, 
                 data = model_df)
summary(log_model)
plot(log_model)

# 1.2 Model: no covariates and no levarge points
if(model_group %in% c("loss", "cna")){
  leverage_free_model_df <- model_df[!(type_agg$Proband %in% c("M00512_04_965_D11")), ]
  log_model_12 <- glm(aggressive ~ .,
                      family = binomial, 
                      data = leverage_free_model_df)
  summary(log_model_12)
  plot(log_model)
}

# 2. model: age as covariate ---------------------------------------------------
# Get age as covariate for model
pheno <-  read.csv("../Pheno_data/2018 10 10 Wien PRCa Update_prepared.csv", 
                   stringsAsFactors = FALSE)
# prepare Proband column in data
probands <- sapply(as.character(type_agg$Proband), function(x) {
  strs <- strsplit(x, split = "_")
  if(nchar(strs[[1]][2]) == 2) {
    return(strs[[1]][3])
  } else {
    return(strs[[1]][2])
  }
})
type_agg_age <- type_agg
type_agg_age$proband <- probands
# merge with pheno data
type_agg_age <- merge(type_agg_age, pheno[, c("proband", "Alter")], 
                      by = "proband", all.x = TRUE)

# new model df
if(model_group %in% c("loss", "cna")){
  model_df_age <- type_agg_age[!(type_agg_age$Proband %in% c("M00512_04_965_D11")), -c(1, 2, 3, 27)] # losses
} else {
  model_df_age <- type_agg_age[, -c(1, 2, 3, 27)] # for gains
} 
model_df_age$aggressive <- factor(model_df_age$aggressive)
# Model
log_model2 <- glm(aggressive ~ .,
                  family = binomial, 
                  data = model_df_age)
summary(log_model2)
plot(log_model2)

# Presented results
odds <- exp(cbind(coef(log_model2), confint(log_model2)))  
ps <- summary(log_model2)
presented_results <- data.frame(round(odds, 3), 
                                "p-value" = round(as.numeric(format(ps$coefficients[, 4], scientific = FALSE)), 4), 
                                check.names = FALSE)
colnames(presented_results)[1] <- "Odds"
write.csv2(presented_results, 
           paste0("D:/Dropbox/CNV Literature/MANUSCRIPT/Table results/", 
                  model_group, "_model_results.csv"),
           row.names = TRUE)

## ROC Curve
jpeg(filename = paste0("D:/Dropbox/CNV Literature/MANUSCRIPT/Plots/", 
                       model_group, "_ROC_curve.jpeg"), width = 800, height = 800, 
     units = "px", quality = 100)
ROC_curve(log_model2)
dev.off()









