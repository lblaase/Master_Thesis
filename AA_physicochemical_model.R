#Predict Ionization efficiencies with measured ionization afficiencies logIE and selected compound parameters


# autoinstall packages
#install.packages("tidyverse","ggcorrplot", "Metrics", "ggplot2", "ggrepel", "RRF")

# include packages
#library(tidyverse, ggcorrplot, Metrics, ggplot2, ggrepel, RRF)

# working directory
setwd("directory")



#set file name
## import
Parameter_file <- "AA_Parameter.csv"
logIE_file <- "logIE_Ratio_MeOH_Glu-d3_QE_H.csv"
## reference compounds
ref <- "d3_Glu"
## output
experiment <- "Ratio_MeOH_Glu-d3_QE"
adduct <- "_H"


##!!! change the parameter to be deleted in line 122 after the correlation matrix


#name of output csv
csv_logIE_MLR_pred = paste("logIE_MLR_pred_", experiment, adduct, ".csv", sep='') # file name for exporting the predicted logIE values
csv_logIE_RRF_pred = paste("logIE_RRF_pred_", experiment, adduct, ".csv", sep='') # file name for exporting the predicted logIE values
csv_model_performance = paste("Model_performance_", experiment, adduct, ".csv",  sep='') # file name for exporting the model evaluation table
csv_model_coefficients = paste("Model_coefficients_", experiment, adduct, ".csv", sep='') # file name for exporting the model coefficients tables
png_logIE_MLR_correlation = paste("Plot_logIE_MLR_correlation_", experiment, adduct, ".png", sep='') # file name for exporting the correlation plot
png_logIE_RRF_correlation = paste("Plot_logIE_RRF_correlation_", experiment, adduct, ".png", sep='') # file name for exporting the correlation plot
png_importance_plot = paste("Plot_parameter_importance_", experiment, adduct, ".png", sep='') # file name for exporting the importance


#read parameter file
Parameter <- read_delim(Parameter_file, delim = ";", escape_double = FALSE, na = "empty", trim_ws = TRUE)
#deleting Ile & Leu from parameters
Parameter <- Parameter[!grepl("Ile|Leu|Arg|His|Lys|Cys", Parameter$AA),]
#read logIE file of measured logIE
logIE <- read_delim(logIE_file, delim = ",", escape_double = FALSE, na = "empty", trim_ws = TRUE)
logIE <- logIE[,-c(1)]
logIE[logIE == "NA"] <- NA

#deleting Ile & Leu from ligIE
delete_AA <- which(names(logIE)%in%c("Ile", "Leu", "d3_Glu", "Arg", "His", "Lys", "Cys"))
logIE <- logIE[,-delete_AA]

#save rownames
rownames_allAA <- Parameter$AA
All_Parameter <- Parameter[,-1]



############################## Scaled Parameter ############################


#calculate scaled parameter
##calculate mean & stdev of each parameter
### create table for mean & stdev
columns <- c("Mean", "Stdev")
Parameter_names <- unlist(colnames(Parameter[,2:10]))
Parameter_mean_stdev <- data.frame(matrix(nrow = 9, ncol = length(columns))) 
colnames(Parameter_mean_stdev) = columns
### calculate mean
Parameter_mean_stdev$Mean <- colMeans(Parameter[sapply(Parameter, is.numeric)], na.rm=TRUE)
### calulate stdev
Parameter_mean_stdev$Stdev <- sapply(Parameter[,2:10], sd, na.rm=TRUE)
### change rownames to parameter name
rownames_Parameter <- c(colnames(Parameter[2:10]))
rownames(Parameter_mean_stdev) <- rownames_Parameter
Parameter_mean_stdev <- t(Parameter_mean_stdev)
Parameter_mean_stdev <- cbind(rownames(Parameter_mean_stdev), data.frame(Parameter_mean_stdev, row.names=NULL))

## calculate scaled parameter
### create empty table
Scaled_Parameter <- Parameter
Scaled_Parameter[,2:10] <- NA
### scaled parameter
parameter = 1
for (parameter in 2:ncol(Parameter)) {
  n_row = 1
  col_name = colnames(Parameter)[parameter]
  for (n_row in 1:nrow(Parameter)) {
    Scaled_Parameter[n_row, parameter] <- (Parameter[n_row, parameter]-Parameter_mean_stdev[1, parameter])/Parameter_mean_stdev[2, parameter]
  }
}




####################### MLR #######################


#Get coefficients for the models with multiple linear regression
##calculate mean of logIE for each AA
logIE_mean <- data.frame(colMeans(logIE[sapply(logIE, is.numeric)], na.rm=TRUE))
colnames(logIE_mean) <- "logIE"
##copy the mean logIE_ values into All_Parameter by row-matching
MLR_table <- Scaled_Parameter[,-1]
rownames(MLR_table) <- rownames_allAA
MLR_table$logIE <- logIE_mean[match(rownames(MLR_table), rownames(logIE_mean)),"logIE"]
MLR_table <- cbind(Parameter[1], data.frame(MLR_table, row.names=NULL))

##get rownumber of rows with NA
row_number <- data.frame(which(is.na(MLR_table), arr.ind=TRUE))
no <- row_number$row

##delete rows with na in logIE_ that was measured
MLR_table <- na.omit(MLR_table)

##create correlation matrix between parameter
corr_matrix_parameter = round(cor(MLR_table[,2:10]), 2)
ggcorrplot(corr_matrix_parameter, hc.order = TRUE, type = "lower", lab = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), )

##when parameter have a correlation >0,8 keep only one of them
MLR_table_corr <- MLR_table[, -which(names(MLR_table) %in% c("M_1", "A", "logMV"))]

##full model
full_model <- lm(logIE ~., data = MLR_table_corr[-1])
full_MLR <- data.frame(summary(full_model)$coefficients)
full_MLR <- full_MLR[-1,1:2]

##Path3-model from Oss, M. et al. (2021): Quantitative electrospray ionization efficiency scale: 10 years after. In: Rapid communications in mass spectrometry : RCM 35 (21), e9178. DOI: 10.1002/rcm.9178.
Path3_model <- lm(logIE ~ PSA + pKa_H2O + logMV, data = MLR_table)
Path3_MLR <- data.frame(summary(Path3_model)$coefficients)
Path3_MLR <- Path3_MLR[-1,1:2]

##AA-model from Oss, M. et al. (2021): Quantitative electrospray ionization efficiency scale: 10 years after. In: Rapid communications in mass spectrometry : RCM 35 (21), e9178. DOI: 10.1002/rcm.9178.
AA_model <- lm(logIE ~ GB + M_1 + logMV, data = MLR_table)
AA_MLR <- data.frame(summary(AA_model)$coefficients)
AA_MLR <- AA_MLR[-1,1:2]

##step model - multiple linear regression model
R_model <- step(full_model, direction = "backward")
R_MLR <- data.frame(summary(R_model)$coefficients)
R_MLR <- R_MLR[-1,1:2]

##extract r2 of all models
summary(full_model)$r.squared
summary(Path3_model)$r.squared
summary(AA_model)$r.squared
summary(R_model)$r.squared


############################# Scaled logIE #########################


#Calculate the scaled parameters
##create table
columns <- c("AA","Path3", "Amino Acid", "R-Model", "Full Model")
Scaled_logIE <- data.frame(matrix(nrow = nrow(Scaled_Parameter), ncol = length(columns))) 
colnames(Scaled_logIE) = columns
Scaled_logIE[,1] <- rownames_allAA

##Path3
n_row = 1
for(n_row in 1:nrow(Scaled_Parameter)){
  Parameter_name <- rownames(Path3_MLR)
  i <- Path3_MLR[,"Estimate"]
  Scaled_logIE[n_row, "Path3"] <- sum(i*Scaled_Parameter[n_row, Parameter_name])
}

##Amino Acid
n_row = 1
for(n_row in 1:nrow(Scaled_Parameter)){
  Parameter_name <- rownames(AA_MLR)
  i <- AA_MLR[,"Estimate"]
  Scaled_logIE[n_row, "Amino Acid"] <- sum(i*Scaled_Parameter[n_row, Parameter_name])
}

##R-Model
n_row = 1
for(n_row in 1:nrow(Scaled_Parameter)){
  Parameter_name <- rownames(R_MLR)
  i <- R_MLR[,"Estimate"]
  Scaled_logIE[n_row, "R-Model"] <- sum(i*Scaled_Parameter[n_row, Parameter_name])
}

##Full-Model
n_row = 1
for(n_row in 1:nrow(Scaled_Parameter)){
  Parameter_name <- rownames(full_MLR)
  i <- full_MLR[,"Estimate"]
  Scaled_logIE[n_row, "Full Model"] <- sum(i*Scaled_Parameter[n_row, Parameter_name])
}


######################### predicted logIE ########################


#Calculate the predicted ligIE values
##create table
logIE_pred <- Scaled_logIE
logIE_pred[,2:5] <- NA

##mean and stdev over all logIE
logIE <- as.matrix(logIE)
logIE <- matrix(as.numeric(logIE),ncol = ncol(logIE))
mean_logIE <- mean(logIE, na.rm=TRUE)
stdev_logIE <- sd(logIE, na.rm=TRUE)

##calculation of the predicted logIE
n_col = 1
for(n_col in 2:ncol(Scaled_logIE)){
  n_row = 1
  for(n_row in 1:nrow(Scaled_logIE)){
    logIE_pred[n_row, n_col] <- mean_logIE + Scaled_logIE[n_row, n_col] * stdev_logIE
  }
}

logIE_pred <- logIE_pred[-no,]



######################### Reggularized Random Forest ###########################


##table for measured logIE
logIE_measured <- matrix(logIE_mean[1:nrow(logIE_mean),])
colnames(logIE_measured) <- "logIE"
logIE_measured <- cbind(logIE_pred[1], data.frame(logIE_measured, row.names=NULL))

#model with regularized random forest
##set seed for random split
set.seed(123)


##RRF without split sampleset with and without importance
RRF_all1 <- RRF(logIE ~., data = MLR_table_corr[-1])
RRF_all2 <- RRF(logIE ~., data = MLR_table_corr[-1], importance = TRUE)
RRF_pred1_all <- data.frame(RRF_all1$predicted)
RRF_pred2_all <- data.frame(RRF_all2$predicted)

##create dataframe with importance value by accuracy
imp_all = as.data.frame(importance(RRF_all2))
imp_all = cbind(vars=rownames(imp_all), imp_all)
imp_all = imp_all[,1:2]
imp_all = as.data.frame(imp_all)
imp_all = imp_all[order(imp_all[,2], decreasing = FALSE),]

  ##plot importance
ggplot(imp_all, aes(x = imp_all[,2], y = imp_all[,1], label=sprintf("%0.2f", round(imp_all[,2], digits = 2)))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  scale_y_discrete(limits = imp_all[,1]) + 
  xlab("Importance in [%]") + ylab("Parameter") +
  ggtitle("Importance full model") +
  geom_text(label = round(imp_all[,2],2), hjust = 0.7, nudge_x = -.5, colour = "black") +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = c(.05, .05)) +
  scale_fill_identity(guide = "none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(size = 18), axis.text = element_text(size = 16), legend.title = element_blank(), legend.text = element_text(size = 14),
        plot.background = element_rect(fill='transparent', color=NA), panel.background = element_rect(fill='transparent'), 
        axis.line = element_line(size = 0.4), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#print importance plot as png
ggsave(png_importance_plot, dpi = 400, width = 20, height = 20, units = "cm")


logIE_RRF_pred <- data.frame(logIE_measured, RRF_pred1_all, RRF_pred2_all)


colnames(logIE_RRF_pred) <- c("AA", "logIE measured", "RRF", "RRF importance")


########################## perfomance evaluation #######################


#Calculate the RMSE and R2 to evaluate the model performance
##table for evaluation values
Model_eval <- data.frame(matrix(ncol = 7, nrow = 4))
colnames(Model_eval) <- c("","Path3", "Amino Acid", "R-Model", "Full Model", "RRF", "RRF importance")
Model_eval[,1] <- c("RMSE", "R2", "Slope", "Intercept")


##RMSE & R2 between predicted and measured logIE
n_col = 1
for(n_col in 2:ncol(logIE_pred)){
  Model_eval[1, n_col] <- rmse(logIE_measured[,2], logIE_pred[, n_col])
  Model_eval[2, n_col] <- summary(lm(logIE_pred[, n_col] ~ logIE_measured[,2]))$r.squared
  Model_eval[3, n_col] <- summary(lm(logIE_pred[, n_col] ~ logIE_measured[,2]))$coefficients[2,1]
  Model_eval[4, n_col] <- summary(lm(logIE_pred[, n_col] ~ logIE_measured[,2]))$coefficients[1,1]
}

n_col = 1
for(n_col in 3:ncol(logIE_RRF_pred)){
  Model_name = colnames(logIE_RRF_pred[n_col])
  data <- na.omit(data.frame(x = logIE_RRF_pred[,2], y = logIE_RRF_pred[, n_col]))
  Model_eval[1, Model_name] <- rmse(data$x, data$y)
  Model_eval[2, Model_name] <- summary(lm(data = data, y ~ x))$r.squared
  Model_eval[3, Model_name] <- summary(lm(data = data, y ~ x))$coefficients[2,1]
  Model_eval[4, Model_name] <- summary(lm(data = data, y ~ x))$coefficients[1,1]
}



low_limit <- min(logIE_pred[,-1])*0.95
max_limit <- max(logIE_pred[,-1])*1.05

##correlation plot between experimental and predicted logIE for MLR
##print the correlation plots as png

x_MLR <- unlist(logIE_measured[,2])
y1_MLR <- unlist(logIE_pred[, 2])
y2_MLR <- unlist(logIE_pred[, 3])
y3_MLR <- unlist(logIE_pred[, 4])
y4_MLR <- unlist(logIE_pred[, 5])

ggplot() + 
  geom_point(data = data.frame(x_MLR, y1_MLR), aes(x = x_MLR, y = y1_MLR, color = "Path3"), size = 2) + 
  geom_point(data = data.frame(x_MLR, y2_MLR), aes(x = x_MLR, y = y2_MLR, color = "Amino Acid"), size = 2) + 
  geom_point(data = data.frame(x_MLR, y3_MLR), aes(x = x_MLR, y = y3_MLR, color = "R-Model"), size = 2) + 
  geom_point(data = data.frame(x_MLR, y4_MLR), aes(x = x_MLR, y = y4_MLR, color = "Full-Model"), size = 2) +
  geom_abline(size = 1) +
  xlab("experimental logIE") + ylab("predicted logIE") +
  xlim(low_limit, max_limit) + ylim(low_limit, max_limit) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(size = 18), axis.text = element_text(size = 16), legend.title = element_blank(), legend.text = element_text(size = 14),
        plot.background = element_rect(fill='transparent', color=NA), panel.background = element_rect(fill='transparent'), 
        panel.border = element_rect(colour="black", fill=NA), legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(png_logIE_MLR_correlation, dpi = 400, width = 20, height = 20, units = "cm")

##correlation plot between experimental and predicted logIE for RRF
##print the correlation plots as png

low_limit_RRF <- min(logIE_RRF_pred[,-1])*0.95
max_limit_RRF <- max(logIE_RRF_pred[,-1])*1.05

x_RRF <- unlist(logIE_RRF_pred[,2])
y1_RRF <- unlist(logIE_RRF_pred[, 3])
y2_RRF <- unlist(logIE_RRF_pred[, 4])

ggplot() + 
  geom_point(data = data.frame(x_RRF, y1_RRF), aes(x = x_RRF, y = y1_RRF, color = "RRF"), size = 2) + 
  geom_point(data = data.frame(x_RRF, y2_RRF), aes(x = x_RRF, y = y2_RRF, color = "RRF importance"), size = 2) + 
  geom_abline(size = 1) +
  xlab("experimental logIE") + ylab("predicted logIE") +
  xlim(low_limit_RRF, max_limit_RRF) + ylim(low_limit_RRF, max_limit_RRF) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(size = 18), axis.text = element_text(size = 16), legend.title = element_blank(), legend.text = element_text(size = 14),
        plot.background = element_rect(fill='transparent', color=NA), panel.background = element_rect(fill='transparent'), 
        panel.border = element_rect(colour="black", fill=NA), legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(png_logIE_RRF_correlation, dpi = 400, width = 20, height = 20, units = "cm")


######################## output #######################


write.csv(logIE_pred, csv_logIE_MLR_pred, row.names=TRUE)
write.csv(logIE_RRF_pred, csv_logIE_RRF_pred, row.names=TRUE)

write.csv(Model_eval, csv_model_performance, row.names=TRUE)
write.csv(rbind(Path3_MLR, AA_MLR, R_MLR, full_MLR), csv_model_coefficients, row.names=TRUE)

