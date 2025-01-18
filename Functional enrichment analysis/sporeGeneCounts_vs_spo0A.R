### Check and plot correlation between count of (significantly-associated) sporulation features and spo0A presence/absence

# Load libraries
library(rms)

# Make sure to load the correct feat.sel.sum in sporulationPrevalenceCorrelations.R first
feat.sel.sum.named <- cbind(feat.sel.prev$original_bin, feat.sel.sum)
#feat.sel.sum.named <- cbind(feat.sel.prev.phylum$original_bin, feat.sel.sum)
colnames(feat.sel.sum.named)[1] <- "original_bin"
plot(feat.sel.sum.named$sum, feat.sel.sum.named$spo0A)
feat_spo0A <- feat.sel.sum.named[c(match("spo0A", colnames(feat.sel.sum.named)), match("sum", colnames(feat.sel.sum.named)))]

print(ggscatter(feat_spo0A, x = "sum", y = "spo0A", add = "reg.line", size = 2, color = rgb(0.25,0.25,0.25,0.05), shape = 16, add.params = list(color = "black", linetype= "dashed")) + theme_bw() +
        stat_cor(label.x = 0.6, label.y = 0.6) + ylim(0, 1) + ggtitle("Bacillota_o_RF39_f_UBA660: sporulation gene counts") + ylab("counts") +
        stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), label.x = 0.6, label.y = 0.4) +
        rremove("grid"))

# Fit logistic regression model
model <- glm(spo0A ~ sum, data=feat_spo0A, family=binomial)
# Define new data frame that contains predictor variable
newdata <- data.frame(sum=seq(min(feat_spo0A$sum), max(feat_spo0A$sum),len=500))
# Use fitted model to predict values of response variable
newdata$spo0A = predict(model, newdata=newdata, type="response")
# Plot logistic regression curve
#pdf(file = paste0("/Users/hl636/Desktop/", "Spo0A_vs_sporeGeneCounts_logRegression_xx.pdf"), width = 5.1, height = 5) #scale 128 in AI%
plot(spo0A ~ sum, data=feat_spo0A, col=rgb(0.25,0.25,0.25,0.08), ylim=c(0,1), pch=16, cex = 1.35, xlab="Sporulation gene counts", ylab="Spo0A presence - absence", main="")
lines(spo0A ~ sum, newdata, lwd=2)
#dev.off()
# Print Nagelkerke's pseudo-R2
mod1b <- lrm(spo0A ~ sum, data=feat_spo0A)
print(mod1b)

## By phylum
phyla_sel <- phyla[c(2,3,4,5)]
phyla_sel <- phyla[c(3)]
for (p in phyla_sel) {
  print(p)
  phylum_df <- phyla_df[phyla_df$phylum == p, ]
  feat.sel.prev.phylum <- merge(phylum_df, feat.sel.prev, by = "original_bin")
  feat.sel.sum <- feat.sel.prev.phylum[,c(2:(ncol(feat.sel.prev.phylum)-1))]
  feat.sel.sum <- feat.sel.sum[,-c(1)]
  feat.sel.sum$sum <- rowSums(feat.sel.sum)
  feat.sel.sum$prevalence <-feat.sel.prev.phylum$prevalence
  feat.sel.sum.named <- cbind(feat.sel.prev.phylum$original_bin, feat.sel.sum)
  colnames(feat.sel.sum.named)[1] <- "original_bin"
  #plot(feat.sel.sum.named$sum, feat.sel.sum.named$spo0A)
  feat_spo0A <- feat.sel.sum.named[c(match("spo0A", colnames(feat.sel.sum.named)), match("sum", colnames(feat.sel.sum.named)))]
  print(ggscatter(feat_spo0A, x = "sum", y = "spo0A", add = "reg.line", size = 2, color = rgb(0.25,0.25,0.25,0.05), shape = 16, add.params = list(color = "black", linetype= "dashed")) + theme_bw() +
          stat_cor(label.x = 0.6, label.y = 0.6) + ylim(0, 1) + ggtitle(paste0(p,": Spo0A vs sporulation gene counts")) + ylab("counts") +
          stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), label.x = 0.6, label.y = 0.4) +
          rremove("grid"))
  # Fit logistic regression model
  model <- glm(spo0A ~ sum, data=feat_spo0A, family=binomial)
  # Define new data frame that contains predictor variable
  newdata <- data.frame(sum=seq(min(feat_spo0A$sum), max(feat_spo0A$sum),len=500))
  # Use fitted model to predict values of response variable
  newdata$spo0A = predict(model, newdata=newdata, type="response")
  # Plot logistic regression curve
  plot(spo0A ~ sum, data=feat_spo0A, col=rgb(0.25,0.25,0.25,0.08), ylim=c(0,1), pch=16, cex = 1.25, xlab="Sporulation gene counts", ylab="Spo0A presence - absence", main=paste0(p, ": Spo0A vs sporulation gene counts"))
  lines(spo0A ~ sum, newdata, lwd=2)
  # Print Nagelkerke's pseudo-R2
  mod1x <- lrm(spo0A ~ sum, data=feat_spo0A)
  print(mod1x)
}


# Multiple comparison adjustment on Bacillota phyla p-values
p_vals_bacillota <- c(0.0001,0.0253, NA, 0.0001) # KEGG; Bacilotta_A, Bacilotta_B, Bacilotta_C, Bacilotta ss
p.adjust(p=p_vals_bacillota, method = "fdr", n=length(p_vals_bacillota))


