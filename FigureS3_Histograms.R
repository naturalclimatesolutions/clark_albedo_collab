
# install libs 
library(dplyr)
library(raster)

# setwd
setwd('C:/Users/samantha.yeo/OneDrive - The Nature Conservancy/Documents/Projects/BEF/Reforestation/Data-SampleSites/DataTables_112023/Final')

# read tabs
combined_nc <- read.csv('finalcount_NetClimateImpact.csv')

# summary & histos 
# Count values for Net Climate Imp
summary_table_verif <- albOffset %>%
  summarise(
    negative_lt0 = sum(Value < 0),
    positive_gt0 = sum(Value > 0),
    neutral_eq0 = sum(Value == 0), 
    na = sum(Value == "-nan(ind)")
  )
# View(summary_table_verif)


###########################
###### Figure S3: NCI #####
###########################

# Drop Missing IDs for albOffset 
combined_nc <- subset(combined_nc, Value != "-nan(ind)")

# Reduce space above x-axis
par(mar = c(5, 4, 4, 4)) 

## count classes - NET CLIMATE IMPACT
df_nc = combined_nc
df_nc$Value <- as.numeric(df_nc$Value)

## hist of final count data - Net Climate Impact
# Set the scientific notation for the Y-axis labels
raster::hist(df_nc$Value, main = "Net Climate Impact", xlab = "Mg CO2e ha^-1", ylab = "Pixel Count (10^4)",
             breaks = 16,  col = "#B2ABD2", border = "white", las = 1,
             cex.axis = .8, cex.lab = .8, cex.main = 1.2, cex.sub = 0.8
             , xlim = c(-500, 1000)
             , ylim = c(0, 350000) # Net Climate Impact
             , axes = FALSE
) ; abline(v = 0, col = "red", lty = "dotted")

# Add y-axis label without tick labels
axis(1, at = NULL, labels = T)

y_labels <- axTicks(2) / 1000
axis(2, at = axTicks(2), labels = y_labels)

# Add the cumulative distribution curve to the plot
# Calculate the empirical cumulative distribution function
data <- as.numeric(df_nc$Value)
ecdf_data <- ecdf(data)

# Generate x-values for plotting
x <- seq(min(data), max(data), length.out = 10000)

# Calculate y-values using the cumulative distribution function
y <- ecdf_data(x)
lines(x, y*295000, type = "l", col = "#B35806", lwd = 2) # net climate impact

# side axis !
par(new = TRUE)
plot(NULL, xlim = c(0, .5), ylim = c(0, 1.0), axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2), las = 1, cex.axis = 0.8)
mtext("Cumulative Proportion", side = 4, line = 2.5, cex = 0.8)


