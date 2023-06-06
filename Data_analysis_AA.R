# autoinstall packages
#install.packages("readr", "randomcoloR", "tidyverse", "stringr")

# include packages
#library(readr, randomcoloR, tidyverse)

# working directory
setwd("C:/Users/lenab/OneDrive - Uppsala universitet/Master Thesis/Data Analysis")
#setwd("C:/Users/LocalUser/Documents/LenaBlaase/Data Analysis")


# defines to control all variables
ign = c()# ignored rows from data set
concentration_file = "NoRatio_IS_Concentration.csv" # file name with concentrations (x-axis)
signal_file = "NoRatio_IS_signal.csv" # file name with signals (y-axis)
experiment = "NoRatio_IS" # experiment name (also used for file names)
reference = "Glu" # used to calculate RIE
delete_AA = c("Ile", "Leu") # delete certain amino acids from data sets



#plot title
main_name_H = paste("H-adducts -", gsub('[_]', ' ', experiment)) # plot title H-adducts
main_name_Na = paste("Na-adducts - ", gsub('[_]', ' ', experiment)) # plot title Na-adducts

#file names
csv_summary_file_H = paste("Summary_", experiment, "_H.csv", sep='') # file name for exporting the summary table
csv_summary_file_Na = paste("Summary_", experiment, "_Na.csv", sep='') # file name for exporting the summary table
csv_RF_file_H = paste("RF_", experiment, "_H.csv", sep='') # file name for exporting the RF table
csv_RF_file_Na = paste("RF_", experiment, "_Na.csv", sep='') # file name for exporting the RF table
csv_logIE_file_H = paste("logIE_", experiment, "_H.csv", sep='') # file name for exporting the logIE table
csv_logIE_file_Na = paste("logIE_", experiment, "_Na.csv", sep='') # file name for exporting the logIE table
pdf_full_plot_H = paste("Full_Plot_", experiment, "_H.pdf", sep='') # file name for exporting the individual AA plots
pdf_full_plot_Na = paste("Full_Plot_", experiment, "_Na.pdf", sep='') # file name for exporting the individual AA plots
pdf_indiv_plot_H = paste("Indiv_Plot_", experiment, "_H.pdf", sep='') # file name for exporting the individual AA plots
pdf_indiv_plot_Na = paste("Indiv_Plot_", experiment, "_Na.pdf", sep='') # file name for exporting the individual AA plots


# importing the calibration concentration table
## if no concentration was measured, leave cell empty
Cal_conc <- read_delim(concentration_file, delim = ";", escape_double = FALSE, na = "empty", trim_ws = TRUE)
#Cal_conc_copy<-Cal_conc[,1]
Cal_conc <- Cal_conc[, -1]
#rownames(Cal_conc) = unlist(Cal_conc_copy[,1]) # to check which rows were deleted
if(length(ign) > 0) {
	Cal_conc <- Cal_conc[-ign,]
}
#delete certain AA
delete <- which(names(Cal_conc)%in%delete_AA) # define AA to be deleted
Cal_conc <- Cal_conc[,-delete] # delete AA

# importing the signal table
## if no signal was measured, leave cell empty, delete rows from blanks
### ignore rows that are specified at the top
Cal_signal <- read_delim(signal_file, delim = ";", escape_double = FALSE, na = "empty", trim_ws = TRUE)
Cal_signal <- filter(Cal_signal, !grepl('Blank', name))
if(length(ign) > 0) {
  Cal_signal <- Cal_signal[-ign,]
}

# create dataframe for signals of protonated (H) and sodiated (Na) adducts
Cal_H_signal <- select(Cal_signal, ends_with("H"))
colnames(Cal_H_signal) <- gsub('\\s|\\s[+H]', '',colnames(Cal_H_signal))
#delete certain AA
delete <- which(names(Cal_H_signal)%in%delete_AA) # define AA to be deleted
Cal_H_signal <- Cal_H_signal[,-delete] # delete AA

Cal_Na_signal <- select(Cal_signal, ends_with("Na"))
colnames(Cal_Na_signal) <- gsub('\\N[a]', '',colnames(Cal_Na_signal))
colnames(Cal_Na_signal) <- gsub('\\s|\\s[+]', '',colnames(Cal_Na_signal))
#delete certain AA
delete <- which(names(Cal_Na_signal)%in%delete_AA) # define AA to be deleted
Cal_Na_signal <- Cal_Na_signal[,-delete] # delete AA


# reading the max concentration and max signal for the creation of the empty plot
xmax <- max(Cal_conc, na.rm = TRUE)*1.2
ymax_H <- max(Cal_H_signal, na.rm = TRUE)*1.2
ymax_Na <- max(Cal_Na_signal, na.rm = TRUE)*1.2

# crate the summary table to store R²,slope, RIE, LOD and LOQ of each curve
column_names = colnames(Cal_conc)
summary_H <- data.frame(AA = vector(length=ncol(Cal_conc)),
						            R2 = vector(mode="double", length=ncol(Cal_conc)),
						            slope = vector(mode="double", length=ncol(Cal_conc)),
						            intercept = vector(mode="double", length=ncol(Cal_conc)),
					            	RIE_slope = vector(mode="double", length=ncol(Cal_conc)),
						            logIE_slope = vector(mode="double", length=ncol(Cal_conc)),
						            LOD = vector(mode="double", length=ncol(Cal_conc)),
                        LOQ = vector(mode="double", length=ncol(Cal_conc)))


summary_Na <- data.frame(AA = vector(length=ncol(Cal_conc)),
                        R2 = vector(mode="double", length=ncol(Cal_conc)),
                        slope = vector(mode="double", length=ncol(Cal_conc)),
                        intercept = vector(mode="double", length=ncol(Cal_conc)),
                        RIE_slope = vector(mode="double", length=ncol(Cal_conc)),
                        logIE_slope = vector(mode="double", length=ncol(Cal_conc)),
                        LOD = vector(mode="double", length=ncol(Cal_conc)),
                        LOQ = vector(mode="double", length=ncol(Cal_conc)))


RF_H <- Cal_H_signal[0,]
RF_Na <- Cal_Na_signal[0,]

logIE_H <- RF_H[0,]
logIE_Na <- RF_Na[0,]

set.seed(10)
# create color list
palette <- distinctColorPalette(k = ncol(Cal_conc), altCol = FALSE, runTsne = FALSE)

#pie(rep(1, ncol(Cal_conc)), col=palette) # prints the color palette


###################### H adducts ##########################

##print the full plot as pdf
pdf(pdf_full_plot_H, width = 10, height = 10)

#create a plot with all AA and extract values for the summary table
# create an empty plot
plot_H <- plot(1, main = main_name_H, type = "n", xlab = "c(Amino Acid) in [μmol/L]", ylab = "Intensity", xlim = c(0, xmax), ylim = c(0, ymax_H), xaxs = "i", yaxs = "i")

# loop to read through the dataset and store all values / plot ines
row = 1
#for (n_aa in 1:3) { # use this instead of line 136 to plot only the 3 first columns
for (n_aa in 1:ncol(Cal_H_signal)) {
	y <- unlist(Cal_H_signal[n_aa])
	col_name <- colnames(Cal_H_signal)[n_aa]
		if (is.infinite(suppressWarnings(max(y, na.rm = TRUE))) == FALSE) {
	  x <- unlist(Cal_conc[,col_name])
	  print(col_name)
		model <- lm(y~x)
		if (format(summary(model)$r.squared) != 1) {
		  points(x, y, pch = 19, col=palette[n_aa])
			abline(model, col=palette[n_aa])
			summary_H$AA[row] = paste(col_name, sep="")
			summary_H$R2[row] = format(summary(model)$r.squared)
			summary_H$slope[row] = coef(model)["x"]
			summary_H$intercept[row] = coef(model)["(Intercept)"]
			summary_H$SE_Residuals[row] = as.numeric(format(summary(model)$sigma))
			row = row + 1
		}
	}
}

# cleanup summary Table
summary_H <- summary_H[summary_H$AA != FALSE, ]
rownames(summary_H) <- unlist(summary_H$AA)
summary_H <- summary_H[, -1]

# add a legend to the plot
legend("topleft", legend=rownames(summary_H), col=palette, cex=1.5, lty=1)

#turn off pdf plotting
dev.off()

#create individual plots
##print the individual plots as pdf
pdf(pdf_indiv_plot_H)

#loop to read through the dataset and plot graphs
row = 1
for (n_aa in 1:ncol(Cal_H_signal)) {
  y <- unlist(Cal_H_signal[n_aa])
  col_name <- colnames(Cal_H_signal)[n_aa]
  if (is.infinite(suppressWarnings(max(y, na.rm = TRUE))) == FALSE) {
    x <- unlist(Cal_conc[,col_name])
    print(col_name)
    title_H <- paste(col_name, "+ H")
    print(ggplot(data = data.frame(x, y), aes(x = x, y = y)) + geom_point(na.rm = TRUE, size = 2) + geom_smooth(method='lm', se = FALSE, fullrange = TRUE)
          + ggtitle(title_H) + xlab("c(Amino Acid) in [μmol/L]") + ylab("Intensity") 
          + scale_x_continuous(expand = c(0, 0), limits = c(0, max(x)*1.1)) + scale_y_continuous(expand = c(0, 0), limits = c(0, max(y)*1.1)) 
          + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16), axis.title = element_text(size = 14), axis.text = element_text(size = 12), 
                  plot.background = element_rect(fill='transparent', color=NA), panel.background = element_rect(fill='transparent'), 
                  panel.border = element_rect(colour="black", fill=NA)), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }
}

#turn off pdf plotting
dev.off()


# calculate RIE_slope
n_aa=1
for(n_aa in 1:nrow(summary_H)) {
	summary_H[n_aa,]$RIE_slope <- summary_H[n_aa,]$slope / summary_H[reference,]$slope
}

# calculate log(IE)_slope
n_aa=1
for(n_aa in 1:nrow(summary_H)) {
  summary_H[n_aa,]$logIE_slope <- log10(summary_H[n_aa,]$slope)
}


#calculate LOD
n_aa=1
for(n_aa in 1:nrow(summary_H)) {
  summary_H[n_aa,]$LOD <- (3.3*summary_H[n_aa,]$SE_Residuals) / summary_H[n_aa,]$slope
}

#calculate LOQ
n_aa=1
for(n_aa in 1:nrow(summary_H)) {
  summary_H[n_aa,]$LOQ <- (10*summary_H[n_aa,]$SE_Residuals) / summary_H[n_aa,]$slope
}

#calculate RF
n_aa = 1
for (n_aa in 1:ncol(Cal_H_signal)) {
  n_row = 1
  #col_name <-substring(colnames(Cal_H_signal)[n_aa], 1, 3) #only if +H or +Na is still in the name
  col_name = colnames(Cal_H_signal)[n_aa]
  #print(col_name)
  for (n_row in 1:nrow(Cal_H_signal)) {
    #print(n_row)
    RF_H[n_row, n_aa]<- Cal_H_signal[n_row, n_aa]/Cal_conc[n_row, col_name]
  }
}

#calculate logIE
n_aa=1
for (n_aa in 1:ncol(RF_H)) {
  n_row = 1
  #print(col_name)
  for (n_row in 1:nrow(RF_H)) {
    #print(n_row)
    logIE_H[n_row, n_aa] <- log10(RF_H[n_row, n_aa])
  }
}


#################### Na adducts ########################

#print the full plot as pdf
pdf(pdf_full_plot_Na, width = 10, height = 10)

#create a plot with all AA and extract values for the summary table
#create an empty plot
plot_Na <- plot(1, main = main_name_Na, type = "n", xlab = "c(Amino Acid) in [μmol/L]", ylab = "Intensity", xlim = c(0, xmax), ylim = c(0, ymax_Na), xaxs = "i", yaxs = "i")

#loop through the data and store values and make plots
row = 1
#for (n_aa in 1:3) { # use this to plot only the 3 first columns
for (n_aa in 1:ncol(Cal_Na_signal)) {
  y <- unlist(Cal_Na_signal[n_aa])
  col_name = colnames(Cal_Na_signal)[n_aa]
  if (is.infinite(suppressWarnings(max(y, na.rm = TRUE))) == FALSE) {
    x <- unlist(Cal_conc[,col_name])
    print(col_name)
    model <- lm(y~x)
    if (format(summary(model)$r.squared) != 1) {
      points(x, y, pch = 19, col=palette[n_aa])
      abline(model, col=palette[n_aa])
      summary_Na$AA[row] = paste(col_name, sep="")
      summary_Na$R2[row] = format(summary(model)$r.squared)
      summary_Na$slope[row] = coef(model)["x"]
      summary_Na$intercept[row] = coef(model)["(Intercept)"]
      summary_Na$SE_Residuals[row] = as.numeric(format(summary(model)$sigma))
      row = row + 1
    }
  }
}

# cleanup summary Table
summary_Na <- summary_Na[summary_Na$AA != FALSE, ]
rownames(summary_Na) <- unlist(summary_Na$AA)
summary_Na <- summary_Na[, -1]

# add a legend to the plot
legend("topleft", legend=rownames(summary_Na), col=palette, cex=1.5, lty=1)

#turn off pdf plotting
dev.off()

#create individual plots
##print the individual plots as pdf
pdf(file = pdf_indiv_plot_Na)

#loop to read through the dataset and plot lines
row = 1
for (n_aa in 1:ncol(Cal_Na_signal)) {
  y <- unlist(Cal_Na_signal[n_aa])
  col_name <- colnames(Cal_Na_signal)[n_aa]
  if (is.infinite(suppressWarnings(max(y, na.rm = TRUE))) == FALSE) {
    x <- unlist(Cal_conc[,col_name])
    print(col_name)
    title_Na <- paste(col_name, "+ Na")
    print(ggplot(data = data.frame(x, y), aes(x = x, y = y)) + geom_point() + geom_smooth(method='lm', se = FALSE, fullrange = TRUE)
          + ggtitle(title_Na) + xlab("c(Amino Acid) in [μmol/L]") + ylab("Intensity")
          + scale_x_continuous(expand = c(0, 0), limits = c(0, max(x)*1.1)) + scale_y_continuous(expand = c(0, 0), limits = c(0, max(y)*1.1))
          + theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(size = 14), axis.text = element_text(size = 12), 
                  plot.background = element_rect(fill='transparent', color=NA), panel.background = element_rect(fill='transparent'), panel.border = element_rect(colour="black", fill=NA)))
  }
}

#turn off pdf plotting
dev.off()


# calculate RIE_slope
n_aa=1
for(n_aa in 1:nrow(summary_Na)) {
  summary_Na[n_aa,]$RIE_slope <- summary_Na[n_aa,]$slope / summary_Na[reference,]$slope
}

# calculate log(IE)_slope
n_aa=1
for(n_aa in 1:nrow(summary_Na)) {
  summary_Na[n_aa,]$logIE_slope <- log10(summary_Na[n_aa,]$slope)
}

#calculate LOD
n_aa=1
for(n_aa in 1:nrow(summary_Na)) {
  summary_Na[n_aa,]$LOD <- (3.3*summary_Na[n_aa,]$SE_Residuals) / summary_Na[n_aa,]$slope
}

#calculate LOQ
n_aa=1
for(n_aa in 1:nrow(summary_Na)) {
  summary_Na[n_aa,]$LOQ <- (10*summary_Na[n_aa,]$SE_Residuals) / summary_Na[n_aa,]$slope
}

#calculate RF
n_aa = 1
for (n_aa in 1:ncol(Cal_Na_signal)) {
  n_row = 1
  #col_name <-substring(colnames(Cal_Na_signal)[n_aa], 1, 3) #only if +H or +Na is still in the name
  col_name = colnames(Cal_Na_signal)[n_aa]
  #print(col_name)
  for (n_row in 1:nrow(Cal_Na_signal)) {
    #print(n_row)
    RF_Na[n_row, n_aa]<- Cal_Na_signal[n_row, n_aa]/Cal_conc[n_row, col_name]
  }
}

#calculate logIE
n_aa=1
for (n_aa in 1:ncol(RF_Na)) {
  n_row = 1
  #print(col_name)
  for (n_row in 1:nrow(RF_Na)) {
    #print(n_row)
    logIE_Na[n_row, n_aa] <- log10(RF_Na[n_row, n_aa])
  }
}


###################### output ####################


# write table to csv file
write.csv(summary_H, csv_summary_file_H, row.names=TRUE)
write.csv(summary_Na, csv_summary_file_Na, row.names=TRUE)
write.csv(RF_H, csv_RF_file_H, row.names=TRUE)
write.csv(RF_Na, csv_RF_file_Na, row.names=TRUE)
write.csv(logIE_H, csv_logIE_file_H, row.names=TRUE)
write.csv(logIE_Na, csv_logIE_file_Na, row.names=TRUE)


