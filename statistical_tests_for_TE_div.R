rm(list = ls())

# Read table 
#df <- read.csv("/home/jens/Viki/RM_Trips_to_Plot/Ppr germany/quantification_TE_haplotypes_20.csv", sep="\t", dec = ",")
df <- read.csv("/home/jens/Viki/RM_Trips_to_Plot/Ppr germany/alt1/May results/quantification of TE figures_20.csv", sep="\t", dec = ",")

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

#-> nicht normalverteilt, kein T-Test möglich

# library(car)
# # Homogeneity of Variances Test
# leveneTest(Chr.1.Alt.1 ~ Chr.1.Alt.2, data = df)
# leveneTest(Chr.2.Alt.1 ~ Chr.2.Alt.2, data = df)
# leveneTest(Chr.3.Alt.1 ~ Chr.3.Alt.2, data = df)
# leveneTest(Chr.4.Alt.1 ~ Chr.4.Alt.2, data = df)
# leveneTest(Chr.5.Alt.1 ~ Chr.5.Alt.2, data = df)
# leveneTest(Chr.6.Alt.1 ~ Chr.6.Alt.2, data = df)
# leveneTest(Chr.7.Alt.1 ~ Chr.7.Alt.2, data = df)
# leveneTest(Chr.8.Alt.1 ~ Chr.8.Alt.2, data = df)
# leveneTest(Chr.9.Alt.1 ~ Chr.9.Alt.2, data = df)

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

# # Perform paired t-test
# t_test_result <- t.test(df$Chr.9.Alt.1, df$Chr.9.Alt.2, paired = TRUE)
# 
# # Print the result
# print(t_test_result)


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

