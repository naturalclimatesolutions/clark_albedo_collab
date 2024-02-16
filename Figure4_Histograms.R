
# install libs 
library(dplyr)
library(raster)

# setwd
setwd('...')


# read tabs
files <- list.files(pattern = '.csv')

albOffset <- read.csv('finalcount_albedo.csv') 
combOpp <- read.csv('finalcount_CompbinedOpportunityNCI.csv') 

# summary & histos 
# Count values for combined opp 
summary_table_verif <- combOpp %>%
  summarise(
    negative_lt0 = sum(Value < 0),
    positive_gt0 = sum(Value > 0),
    neutral_eq0 = sum(Value == 0), 
    na = sum(Value == "-nan(ind)")
  )
# View(summary_table_verif)

###################################
###### Figure 4: Combined Opp #####
###################################

# Drop Missing IDs for albOffset 
albOffset <- subset(albOffset, Value != "-nan(ind)")

# join tabs & filter dat 
by = dplyr::join_by(X)
combined = full_join(albOffset, combOpp, by)
combined = subset(combined, Value.y != "-nan(ind)") %>% dplyr::select(ID.x, Value.x)

# Reduce space above x-axis
par(mar = c(5, 4, 4, 4)) 

## count classes - NET CLIMATE IMPACT
albOffset$Value <- as.numeric(albOffset$Value)
combined$Value <- as.numeric(combined$Value.x)

## histogram of final count data - Alb Offset
df_alb2 <- subset(albOffset, Value >= -100 & Value <= 100) # filter values

raster::hist(df_alb2$Value, main = "Distribution of Albedo Offset", xlab = "Percent", ylab = "Pixel Count (10^4)",
             breaks = 15,  col = "#B2ABD2", border = "white", las = .5,
             cex.axis = .8, cex.lab = .8, cex.main = 1.2, cex.sub = 0.8
             , ylim = c(0, 200000)
             , xlim = c(-50, 100) # symbols flipped
             , axes = FALSE
) 

# Add y-axis label without tick labels
axis(1, at = NULL, labels = T)

y_labels <- axTicks(2) / 1000
axis(2, at = axTicks(2), labels = y_labels)

# side axis
par(new = TRUE)
plot(NULL, xlim = c(0, .5), ylim = c(0, 1.0), axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2), las = 1, cex.axis = 1)

# histogram for combined opportunity 
par(new = TRUE)
df_comb <- subset(combined, Value >= -100 & Value <= 100)
raster::hist(df_comb$Value, main = "", xlab = "", ylab = "",
             breaks = 15,  col = "#726E86", border = "white", las = 1,
             cex.axis = 0.8, cex.lab = 0.8, cex.main = 1.2, cex.sub = 0.8
             , yaxt = 'n' , xaxt = 'n'
             , ylim = c(0, 200000)
             , xlim = c(-50, 100) # symbols flipped in table
             # , add = TRUE
) 

# Add the cumulative distribution curve to the plot
data <- as.numeric(albOffset$Value)
ecdf_data <- ecdf(data)

# Generate x-values for plotting
x <- seq(min(data), max(data), length.out = 10000)

# Calculate y-values using the cumulative distribution function
y <- ecdf_data(x)
lines(x, y*192000, type = "l", col = "#B35806", lwd = 2) # Cumulative Proportion to 87% for NCI 
mtext("Cumulative Proportion", side = 4, line = 3, cex = 0.8) 


