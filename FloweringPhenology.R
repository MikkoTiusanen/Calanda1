#FloweringPhenology
data <- read.csv(file = "PhenologyFlowers.csv", header = TRUE, sep = ",")
flowers <- names(data)[15:137]
fits <- lapply(flowers, function(x) {glm(as.formula(paste(x, "~ elevation")), data = data)})
unfit <- unlist(lapply(fits, function(x) { coef(x)[[2]] }))
# Increase margins before plotting
par(mar = c(5, 6, 4, 2))  # bottom, left, top, right
hist(unfit,
     breaks = 200,
     xlim=c(-0.05,0.1),
     main = "",#"Frequency distribution of phenology shift",
     xlab="Phenology shift (days/m)",
     cex.lab = 2,  
     cex.axis = 2,  
     cex.main = 2 
)
abline(v= 0,lwd = 2)
abline(v= 0.014147,lwd = 3, col = "red")
#(0.014147 = Average of the 123 species-specific values)
