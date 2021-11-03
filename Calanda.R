library(mgcv)
library(brms)
#GAMM Flower abundance vs elevation & DOY####
data <- read.csv(file = "FloweringInsectDensity.csv", header = TRUE, sep = ",")
data[9:177] <- lapply(data[9:177], as.numeric)

gamm(Totalflowers ~ MinOfelevation + s(DOY),random = list(Meadow = ~1,Site = ~1), data = data)

#GAMM Flowering species richness vs elevation & DOY####
gamm(Species ~ MinOfelevation + s(DOY), random = list(Meadow = ~1,Site = ~1), data = data)

#GAMM Pollinator abundance vs elevation & DOY####
gamm(InsectDensity ~ MinOfelevation + s(DOY), random = list(Meadow = ~1,Site = ~1), data = data)

#GAMM Insects per flower vs elevation & DOY####
gamm(InsectsPerFlower ~ MinOfelevation + s(DOY), random = list(Meadow = ~1,Site = ~1), data = data)

#Flowering phenology####
data <- read.csv(file = "PhenologyFlowers.csv", header = TRUE, sep = ",")
flowers <- names(data)[15:137]
fits <- lapply(flowers, function(x) {glm(as.formula(paste(x, "~ elevation")), data = data)})
unfit <- unlist(lapply(fits, function(x) { coef(x)[[2]] }))
hist(unfit,
     breaks = 200,
     xlim=c(-0.05,0.1),
     main = "Frequency distribution of phenology shift",
     xlab="Phenology shift (days/m)" 
)
abline(v= 0,lwd = 2)
abline(v= 0.014147,lwd = 3, col = "red")
#(0.014147 = Average of the 123 species-specific values)

#Niche overlaps along elevation gradient with Schoener's index####
#Within flowering species####
data <- read.csv(file = "Schoener.csv", header = TRUE, sep = ",")
brm(Schoener ~ Elevation + (1|Meadow/Site),
            data = data,
            family = 'zero_one_inflated_beta',
            iter = 1000,
            chains = 2,
            cores = 2)
#Between flowering species and insect abundance####
data <- read.csv(file = "Schoener flower overall insect.csv", header = TRUE, sep = ",")
brm(Schoener ~ Elevation + (1|Meadow/Site),
            data = data,
            family = 'beta',
            iter = 1000,
            chains = 2,
            cores = 2)
#Between flowering species and Coleoptera abundance####
data <- read.csv(file = "Schoener flower Coleoptera.csv", header = TRUE, sep = ",")
brm(Schoener ~ Elevation + (1|Meadow/Site),
            data = data,
            family = 'zero_inflated_beta',
            iter = 1000,
            chains = 2,
            cores = 2)
#Between flowering species and Diptera abundance####
data <- read.csv(file = "Schoener flower Diptera.csv", header = TRUE, sep = ",")
brm(Schoener ~ Elevation + (1|Meadow/Site),
             data = data,
             family = 'beta',
             iter = 1000,
             chains = 2,
             cores = 2)
#Between flowering species and Hymenoptera abundance####
data <- read.csv(file = "Schoener flower Hymenoptera.csv", header = TRUE, sep = ",")
brm(Schoener ~ Elevation + (1|Meadow/Site),
             data = data,
             family = 'beta',
             iter = 1000,
             chains = 2,
            cores = 2)
#Within flowering species with flowering species average phenology####
#The values of elevation and species average phenology were divided by their
#respective averages to ensure model convergence
data <- read.csv(file = "Schoener flower phenology.csv", header = TRUE, sep = ",")
brm(Schoener ~ SSpeciesAveragePhenology * SElevation + (1|Meadow/Site),
             data = data,
             family = 'zero_one_inflated_beta',
             iter = 1000,
             chains = 2,
             cores = 2)