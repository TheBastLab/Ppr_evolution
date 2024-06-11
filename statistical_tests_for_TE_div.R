rm(list = ls())

# Read table 
df <- read.csv("/path/to/quantification_of_TE_figures_delta_differences.csv", sep="\t", dec = ",") #or sep=";", dec = "." etc

# Test for normality (If the p-values for both tests are greater than your chosen significance level (e.g., 0.05), you can assume normality.)
shapiro.test(df$Chr.1.Alt.1)
shapiro.test(df$Chr.2.Alt.1)
shapiro.test(df$Chr.3.Alt.1)
shapiro.test(df$Chr.4.Alt.1)
shapiro.test(df$Chr.5.Alt.1)
shapiro.test(df$Chr.6.Alt.1)
shapiro.test(df$Chr.7.Alt.1)
shapiro.test(df$Chr.8.Alt.1)
shapiro.test(df$Chr.9.Alt.1)

shapiro.test(df$Chr.1.Alt.2)
shapiro.test(df$Chr.2.Alt.2)
shapiro.test(df$Chr.3.Alt.2)
shapiro.test(df$Chr.4.Alt.2)
shapiro.test(df$Chr.5.Alt.2)
shapiro.test(df$Chr.6.Alt.2)
shapiro.test(df$Chr.7.Alt.2)
shapiro.test(df$Chr.8.Alt.2)
shapiro.test(df$Chr.9.Alt.2)

#-> normality not given, no t-test possible

# Homogeneity of Variances Test (If the p-value is greater than your chosen significance level (e.g., 0.05), you can assume homogeneity of variances.)
var.test(df$Chr.1.Alt.1, df$Chr.1.Alt.2)
var.test(df$Chr.2.Alt.1, df$Chr.2.Alt.2)
var.test(df$Chr.3.Alt.1, df$Chr.3.Alt.2)
var.test(df$Chr.4.Alt.1, df$Chr.4.Alt.2)
var.test(df$Chr.5.Alt.1, df$Chr.5.Alt.2)
var.test(df$Chr.6.Alt.1, df$Chr.6.Alt.2)
var.test(df$Chr.7.Alt.1, df$Chr.7.Alt.2)
var.test(df$Chr.8.Alt.1, df$Chr.8.Alt.2)
var.test(df$Chr.9.Alt.1, df$Chr.9.Alt.2)

#-> homogeneity of variances is given

# # Perform paired t-test for the individual chromosomes (if normality and homogeneity are met)
# # For example:
# t_test_result <- t.test(df$Chr.9.Alt.1, df$Chr.9.Alt.2, paired = TRUE)
# 
# # Print the result
# print(t_test_result)

# Since t-test was not possible:
# Perform non-parametric Wilcoxon signed-rank test
wilcox.test(df$Chr.1.Alt.1, df$Chr.1.Alt.2, paired = TRUE)
wilcox.test(df$Chr.2.Alt.1, df$Chr.2.Alt.2, paired = TRUE)
wilcox.test(df$Chr.3.Alt.1, df$Chr.3.Alt.2, paired = TRUE)
wilcox.test(df$Chr.4.Alt.1, df$Chr.4.Alt.2, paired = TRUE)
wilcox.test(df$Chr.5.Alt.1, df$Chr.5.Alt.2, paired = TRUE)
wilcox.test(df$Chr.6.Alt.1, df$Chr.6.Alt.2, paired = TRUE)
wilcox.test(df$Chr.7.Alt.1, df$Chr.7.Alt.2, paired = TRUE)
wilcox.test(df$Chr.8.Alt.1, df$Chr.8.Alt.2, paired = TRUE)
wilcox.test(df$Chr.9.Alt.1, df$Chr.9.Alt.2, paired = TRUE)

# To plot the haplotype divergence:
# Take the table and reduce it to Divergence and the deviation columns of the respective chromosomes, call it df
# Reshape the dataframe from wide to long format
df_long <- tidyr::gather(df, key = "Group", value = "Value", -Divergence)

# Calculate mean values for each x-value across all subsets
mean_values <- aggregate(Value ~ Divergence, data = df_long, FUN = mean)

# Plot using ggplot2
ggplot(df_long, aes(x = Divergence, y = Value, color = Group)) +
  geom_line() +
  geom_line(data = mean_values, aes(x = Divergence, y = Value, color = "Mean"), linetype = "dashed") +  # Add mean line
  #scale_color_brewer(palette = "Set3") +  # Change palette as needed
  #scale_color_manual(values = c(RColorBrewer::brewer.pal(9, "Set3"), "black")) +  # Change to a different palette or specify colors manually
  scale_color_manual(values = c("#D2A271", "#FB7F00", "#F7E689", "#90D08B", "#C1DCDC", "#5BA8E6", "#E2B6FA", "#E27296", "#F15D59",  "black")) +
  theme_minimal() +
  labs(x = "Kimura substitution level", y = "Delta difference between haplotypes", color = "Data Subsets") +
  scale_x_continuous(breaks = seq(min(df$Divergence), max(df$Divergence), by = 1)) +  # Adjust x-axis breaks
  scale_y_continuous(breaks = seq(min(df_long$Value), max(df_long$Value), by = 0.1))  # Adjust y-axis breaks

#######################################################################################################################################################
# Testing via Anova

rm(list = ls())

# These functions are designed to perform analysis of variance (ANOVA) tests, specifically for different numbers of factors. 
# ANOVA tests are used to analyze the differences among group means in a sample.

# To decide which function to use and what arguments to pass, here are some considerations:
# Number of Factors: Determine how many independent variables (factors) you have
# Interactions: Decide whether you're interested in testing interactions between factors
# Dependent Variable: Identify the dependent variable

# anovsF1
# This function is used for ANOVA with one factor.

anovsF1 <- function(response, v1, nbrep = 5000) {
  rep <- response
  var1 <- v1
  data <- data.frame(rep = rep, var1 = var1)
  obs <- data.frame(t(anova(lm(rep ~ var1, data = data))$F))
  names(obs) <- c("Factor 1", "")
  for (i in 2:nbrep) {
    obs[i, ] <- t(anova(lm(sample(rep) ~ var1, data = data))$F)
  }
  pval <- NULL
  for (j in 1:1) {
    pval[j] <- sum(obs[, j] >= obs[1, j]) / length(obs[, j])
  }
  
  return(pval)
}


# anovsF2
# This function is used for ANOVA with two factors.
# It expects three arguments: response, v1, and v2.
# response: This is the dependent variable you want to analyze.
# v1 and v2: These are the independent variables (factors) that you want to test.
# It computes ANOVA for the main effects of each factor (Factor 1 and Factor 2) and returns the p-values for each.

anovsF2 <- function(response, v1, v2, nbrep = 5000) {
  rep <- response
  var1 <- v1
  var2 <- v2
  data <- data.frame(rep = rep, var1 = var1, var2 = var2)
  obs <- data.frame(t(anova(lm(rep ~ var1 + var2, data = data))$F))
  names(obs) <- c("Factor 1", "Factor 2", "")
  for (i in 2:nbrep) {
    obs[i, ] <- t(anova(lm(sample(rep) ~ var1 + var2, data = data))$F)
  }
  pval <- NULL
  for (j in 1:2) {
    pval[j] <- sum(obs[, j] >= obs[1, j]) / length(obs[, j])
  }
  par(mfrow = c(2, 1))
  hist(obs$"Factor 1", xlab = "F")
  abline(v = obs[1, 1], lwd = 2)
  hist(obs$"Factor 2", xlab = "F")
  abline(v = obs[1, 2], lwd = 2)
  
  return(pval)
}


# anovsF2_2
# This function is used for ANOVA with two factors and their interaction.
# It expects three arguments: response, v1, and v2.
# response: This is the dependent variable you want to analyze.
# v1 and v2: These are the independent variables (factors) that you want to test.
# It computes ANOVA for the main effects of each factor (Factor 1 and Factor 2), their interaction, and returns the p-values for each.

anovsF2_2 <- function(response, v1, v2, nbrep = 5000) {
  rep <- response
  var1 <- v1
  var2 <- v2
  data <- data.frame(rep = rep, var1 = var1, var2 = var2)
  obs <- data.frame(t(anova(lm(rep ~ var1 * var2, data = data))$F))
  names(obs) <- c("Factor 1", "Factor 2", "Interaction", "")
  for (i in 2:nbrep) {
    obs[i, ] <- t(anova(lm(sample(rep) ~ var1 * var2, data = data))$F)
  }
  pval <- NULL
  for (j in 1:3) {
    pval[j] <- sum(obs[, j] >= obs[1, j]) / length(obs[, j])
  }
  par(mfrow = c(3, 1))
  hist(obs$"Factor 1", xlab = "F")
  abline(v = obs[1, 1], lwd = 2)
  hist(obs$"Factor 2", xlab = "F")
  abline(v = obs[1, 2], lwd = 2)
  hist(obs$"Interaction", xlab = "F")
  abline(v = obs[1, 3], lwd = 2)
  
  return(pval)
}


# anovsF3
# This function is used for ANOVA with three factors and their interaction.
# It expects four arguments: response, v1, v2, and v3.
# response: This is the dependent variable you want to analyze.
# v1, v2, and v3: These are the independent variables (factors) that you want to test.
# It computes ANOVA for the main effects of each factor (Factor 1, Factor 2, and Factor 3), their interactions, and returns the p-values for each.

anovsF3 <- function(response, v1, v2, v3, nbrep = 5000) {
  rep <- response
  var1 <- v1
  var2 <- v2
  var3 <- v3
  data <- data.frame(rep = rep, var1 = var1, var2 = var2, var3 = var3)
  obs <- data.frame(t(anova(lm(rep ~ var1 + var2 * var3, data = data))$F))
  names(obs) <- c("Factor 1", "Factor 2", "Factor 3", "Interaction", "")
  for (i in 2:nbrep) {
    obs[i, ] <- t(anova(lm(sample(rep) ~ var1 + var2 * var3, data = data))$F)
  }
  pval <- NULL
  for (j in 1:4) {
    pval[j] <- sum(obs[, j] >= obs[1, j]) / length(obs[, j])
  }
  
  return(pval)
}


# Read table (note that the table contains the same information but is in long format)
data <- read.csv("/path/to/TE_div_for_anova.tsv", sep = "\t", dec = ".")

data$divergence <- as.numeric(data$divergence)

# Pass the correct arguments to the functions
anovsF2(data$divergence_value, data$chromosome, data$haplotype)
anovsF2_2(data$divergence_value, data$chromosome, data$haplotype)
#anovsF3(data$divergence_value, data$divergence, data$chromosome, data$haplotype)


# Split the sorted data into smaller sets based on subject_ID
subset_list <- split(data, data$chromosome)

# Separate into individual data frames
chr1 <- subset_list[[1]]
chr2 <- subset_list[[2]]
chr3 <- subset_list[[3]]
chr4 <- subset_list[[4]]
chr5 <- subset_list[[5]]
chr6 <- subset_list[[6]]
chr7 <- subset_list[[7]]
chr8 <- subset_list[[8]]
chr9 <- subset_list[[9]]

# Test chromosomes individually
anovsF1(chr1$divergence_value, chr1$haplotype)
anovsF1(chr2$divergence_value, chr2$haplotype)
anovsF1(chr3$divergence_value, chr3$haplotype)
anovsF1(chr4$divergence_value, chr4$haplotype)
anovsF1(chr5$divergence_value, chr5$haplotype)
anovsF1(chr6$divergence_value, chr6$haplotype)
anovsF1(chr7$divergence_value, chr7$haplotype)
anovsF1(chr8$divergence_value, chr8$haplotype)
anovsF1(chr9$divergence_value, chr9$haplotype)





