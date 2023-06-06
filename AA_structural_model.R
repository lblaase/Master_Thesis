#Predict Ionization efficiencies with measured ionization afficiencies logIE and selected compound parameters


# autoinstall packages
#install.packages("tidyverse","ggcorrplot", "Metrics", "ggplot2", "ggrepel", "RRF")

# include packages
#library(tidyverse, ggcorrplot, Metrics, ggplot2, ggrepel, RRF)

# working directory
setwd("directory")



#set file name
## import
Parameter_file <- "AA_PADEL_Parameter.csv"
logIE_file <- "logIE_RBE_MeOH_Glu-d3_H.csv"
# Ignore AA
ign_AA <- c("Ile", "Leu", "d3_Glu", "Arg", "His", "Lys", "Cys")
## output
experiment <- "RBE_MeOH_Glu-d3"
adduct <- "_H"

#name of output files
csv_logIE_RRF_pred = paste("logIE_RRF_structure_pred_", experiment, adduct, ".csv", sep='') # file name for exporting the predicted logIE values
csv_model_performance = paste("RRF_structure_performance_", experiment, adduct, ".csv", sep='') # file name for exporting the model evaluation table
png_logIE_correlation = paste("Plot_logIE_correlation_RRF_structure_", experiment, adduct, ".png", sep='') # file name for exporting the individual AA plots



################## prepare the parameter file ######################

#read parameter file
Parameter <- read_delim(Parameter_file, delim = ";", escape_double = FALSE, na = "empty", trim_ws = TRUE)
Parameter_red <- Parameter[ , colSums(is.na(Parameter))==0]
Parameter_red <- subset(Parameter_red, !(Name %in% ign_AA))

# Find the most common value in each column
most_common <- sapply(Parameter_red, function(x) {
  tab <- table(x)
  names(tab)[which.max(tab)]
})

# Calculate the proportion of the most common value in each column
prop_most_common <- sapply(Parameter_red, function(x) {
  tab <- table(x)
  max_val <- max(tab)
  max_prop <- max_val / sum(tab)
  return(max_prop)
})

# Save the respective column names where the most common value is in at least 95% of the rows
selected_cols <- c(names(which(prop_most_common >= 0.95)))

# Delete the respective columns
Parameter_red <- Parameter_red[, !names(Parameter_red) %in% selected_cols]

# make a correlation matrix between all parameters
corr_matrix = round(cor(Parameter_red[,-1]), 2)

# Identify the unique pairs of highly correlated columns
high_corr_pairs <- which(abs(corr_matrix) > 0.8 & upper.tri(corr_matrix, diag = FALSE), arr.ind = TRUE)

# Remove variable with higher average correlation
for(i in seq_len(nrow(high_corr_pairs))) {
  var1 <- rownames(corr_matrix)[high_corr_pairs[i, 1]]
  var2 <- colnames(corr_matrix)[high_corr_pairs[i, 2]]
  
  avg_corr_var1 <- mean(abs(corr_matrix[var1, ]))
  avg_corr_var2 <- mean(abs(corr_matrix[var2, ]))
  
  if(avg_corr_var1 > avg_corr_var2) {
    Parameter_red <- Parameter_red[, !grepl(var1, names(Parameter_red))]
  } else {
    Parameter_red <- Parameter_red[, !grepl(var2, names(Parameter_red))]
  }
}


##################### prepare the log IE file ############################


#read logIE file of measured logIE
logIE <- read_delim(logIE_file, delim = ",", escape_double = FALSE, na = "empty", trim_ws = TRUE)
logIE <- logIE[,-c(1)]
logIE[logIE == "NA"] <- NA

#deleting AA from logIE
delete_AA <- which(names(logIE) %in% ign_AA)
logIE <- logIE[,-delete_AA]

# Calculate the mean logIE of each AA
logIE_mean <- data.frame(colMeans(logIE[sapply(logIE, is.numeric)], na.rm=TRUE))
colnames(logIE_mean) <- "logIE"
logIE_mean$Name <- rownames(logIE_mean)


################## RRF ###################
set.seed(123)

# parameter file for RRF
Parameter_RRF <- merge(logIE_mean[], Parameter_red, by = 'Name', all = TRUE)
Parameter_RRF <- na.omit(Parameter_RRF)

#RRF without and with importance
RRF_all1 <- RRF(logIE ~., data = Parameter_RRF[-1])
RRF_all2 <- RRF(logIE ~., data = Parameter_RRF[-1], importance = TRUE)
RRF_pred1_all <- data.frame(RRF_all1$predicted)
RRF_pred2_all <- data.frame(RRF_all2$predicted)


low_limit <- min(logIE_mean[,1])*0.95
max_limit <- max(logIE_mean[,1])*1.05

#plot the correlation between measured and predicted logIE
x_RRF <- unlist(logIE_mean[,1])
y1_RRF <- unlist(RRF_pred1_all[,1])
y2_RRF <- unlist(RRF_pred2_all[,1])

ggplot() + 
  geom_point(data = data.frame(x_RRF, y1_RRF), aes(x = x_RRF, y = y1_RRF, color = "RRF"), size = 2) + 
  geom_point(data = data.frame(x_RRF, y2_RRF), aes(x = x_RRF, y = y2_RRF, color = "RRF importance"), size = 2) + 
  geom_abline(size = 1) +
  xlab("experimental logIE") + ylab("predicted logIE") +
  xlim(low_limit, max_limit) + ylim(low_limit, max_limit) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(size = 18), axis.text = element_text(size = 16), legend.title = element_blank(), legend.text = element_text(size = 14),
        plot.background = element_rect(fill='transparent', color=NA), panel.background = element_rect(fill='transparent'), 
        panel.border = element_rect(colour="black", fill=NA), legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#print importance plot as png
ggsave(png_logIE_correlation, dpi = 400, width = 20, height = 20, units = "cm")


#################### Data Summary #######################

 
RRF_struc_logIE_pred <- data.frame(y1_RRF, y2_RRF)
colnames(RRF_struc_logIE_pred) <- c("RRF", "RRF_imp")

RRF_struc_Eval <- data.frame(matrix(nrow = 2, ncol = 2))
colnames(RRF_struc_Eval) <- c("RRF", "RRF_imp")
rownames(RRF_struc_Eval) <- c("RMSE", "R2")

RRF_struc_Eval[1,1] <- rmse(x_RRF, y1_RRF)
RRF_struc_Eval[1,2] <- rmse(x_RRF, y2_RRF)
RRF_struc_Eval[2,1] <- summary(lm(y1_RRF ~ x_RRF))$r.squared
RRF_struc_Eval[2,2] <- summary(lm(y2_RRF ~ x_RRF))$r.squared


############### output ############

write.csv(RRF_struc_logIE_pred, csv_logIE_RRF_pred, row.names=TRUE)
write.csv(RRF_struc_Eval, csv_model_performance, row.names=TRUE)

