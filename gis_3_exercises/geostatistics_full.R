# Package names
packages <- c("knitr", "sp", "gstat", "dplyr", "ggfortify", "sf", "tmap")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


# -------------------------------------

# Import the Meuse data set of the sp package including the shape of the river, the study area, and the grid
data(meuse, package = "sp")
data(meuse.riv, package = "sp")
data(meuse.area, package = "sp")
data(meuse.grid, package = "sp")

# Show the header of the data set
names(meuse)

# Display the internal structure of the meuse object
str(meuse)

# In RStudio: View the data set as a table
# View(meuse)

# Plot the data set by using ggplot
ggplot(meuse, aes(x,y)) +
  geom_point(aes(color = zinc), size = 2, alpha = 3/4) +
  ggtitle("Zinc Concentration (ppm)") + coord_equal() + theme_bw()

# Convert the Meuse data frame into an sf (simple feature) object
meuse_sf <- st_as_sf(as.data.frame(meuse), coords = c("x", "y"), crs = 28992, agr = "constant")
meuse_sf[1:3,]

# Plot data using quick tmap
qtm(meuse_sf, symbols.col = "zinc", symbols.size = 0.5)

# Designate the coordinates in the data frames as SpatialPointsDataFrame. For the grid data set, assert that it is gridded.
coordinates(meuse) <- c("x","y")
coordinates(meuse)[1:4,]
class(meuse)
summary(meuse)

# Display the data using spplot: Display the attribute "zinc" as z value
spplot(meuse, zcol="zinc", main = "zinc concentrations (ppm)")

# Create a spatial polygon object from the the river data
meuse.sr <- SpatialPolygons(list(Polygons(list(Polygon(meuse.riv)), "meuse.riv")))

# Now, put the legend to the right, set own class breaks, and add the shape of the river
spplot(meuse, zcol="zinc", col.regions=bpy.colors(),  main = "zinc concentrations (ppm)",
       cuts = c(100,200,400,700,1200,2000), key.space = "right",
       sp.layout= list("sp.polygons", meuse.sr, fill = "lightblue")
)

# Histogram of zinc concentration data
hist(meuse$zinc, breaks = 15, main = "Histogram of zinc conc. (ppm)")
hist(log(meuse$zinc), breaks = "scott", main = "Histogram of log(zinc concentration)")

# The Meuse data set provides a grid for the interpolation, prediction, and simulation. First, check its structure.
summary(meuse.grid)
class(meuse.grid)
str(meuse.grid)

# Display the prediction grid and overlay the sample data locations
p1 <- ggplot(data = as.data.frame(meuse.grid), aes(x,y)) + geom_point(size = 0.5) +
  coord_equal() + ggtitle("Grid nodes and sample data locations")
p1 +  geom_point(data = as.data.frame(meuse), aes(x,y), shape = "*", size = 5, color = "red")


# -------------------------------------

# Turn the prediction grid into a SpatialPointsDataFrame
coordinates(meuse.grid) <- c("x", "y")
class(meuse.grid)
str(meuse.grid)

# Turn the SpatialPointsDataFrame into a SpatialPixelsDataFrame
gridded(meuse.grid) <- TRUE
class(meuse.grid)

# Display the distance to the river with ggplot
ggplot(as.data.frame(meuse.grid), aes(x, y)) + geom_tile(aes(fill=dist)) +
  scale_fill_gradient(low = "red", high="yellow") + coord_equal() + theme_bw() +
  ggtitle("Distance to river Meuse")

# Display the distance to the river with spplot
spplot(meuse.grid["dist"], main="Distance from river Meuse")

# -------------------------------------

# # For the following IDW interpolation in R, we assume a linear spatial process
# idw.out <- idw(zinc ~ 1, meuse, meuse.grid, idp = 2.5)
# class(idw.out)
# spplot(idw.out, "var1.pred", main = "zinc IDW predictions")


# Now, for the following IDW interpolation, we assume a lognormal spatial process, as this is commonly used in soil science. The predicted values will be stored in the attribute var1.pred.
logzinc.idw <- idw(log(zinc)~1, meuse, meuse.grid, idp = 2.5)
class(logzinc.idw)
spplot(logzinc.idw, "var1.pred", main = "log(zinc) IDW predictions")

# Plot the histogram of the predicted values.
hist(logzinc.idw$var1.pred, xlab="IDW interpolated log(zinc)")
hist(exp(logzinc.idw$var1.pred), xlab="Back-transformed IDW zinc values")


# -------------------------------------

# Task

# The next step shows the back-transformation of the lognormally transformed predicted values. Question: Is the use of the exponential function a feasible way to back-transform the lognormally transfromed predicted values or is something missing?
zincfromlog.idw <- logzinc.idw
zincfromlog.idw$var1.pred = exp(logzinc.idw$var1.pred)
spplot(zincfromlog.idw,"var1.pred", main = "zinc from log IDW interpolations")

# Compare with IDW of zinc concentrations
zinc.idw = idw(zinc~1, meuse, meuse.grid)
spplot(zinc.idw,"var1.pred", main = "zinc IDW interpolation")

# Compare again using spplot
spplot(meuse,zcol="zinc", main = "zinc concentrations (ppm)", key.space = "right",
       sp.layout= list("sp.polygons", meuse.sr, fill = "lightblue")
       )

# Commonly, the inverse distance weighting power is 2, which means that the distance is squared. By setting the idp parameter manually, you can modify the IDW distance function.
logzinc.idw1 = idw(log(zinc)~1, meuse, meuse.grid, idp = 1)

## [inverse distance weighted interpolation]
## idp = 2 or idp = 3 -> power parameter for distance
spplot(logzinc.idw1,"var1.pred", main = "log(zinc) IDW interpolation power=1")

# Compare histograms of predicted values
hist(zincfromlog.idw$var1.pred, xlab="zinc from log IDW predictions")
hist(zinc.idw$var1.pred, xlab="zinc IDW predictions")

# Calculate the variances of the predicted data sets
c(var(zincfromlog.idw$var1.pred),var(zinc.idw$var1.pred))


# -------------------------------------

# IDW cross-validation
# create an IDW "object" in gstat
logzincIDW1.obj <- gstat(id="logzincIDW1.obj", formula = log(zinc) ~ 1, data = meuse,
                         nmax=12, set=list(idp=1)) # idp is the power exponent

# Perform cross-valiation and determine the RMSE
logzincIDW1.cv <- gstat.cv(logzincIDW1.obj, debug.level=0, random=FALSE)
str(logzincIDW1.cv) # list contents
rmse.idp1 <- sqrt(mean(logzincIDW1.cv$residual^2)) # RMS prediction error

# Repeat for idp=2
logzincIDW2.obj <- gstat(id="logzincIDW2.obj", formula = log(zinc) ~ 1, data = meuse,
                         nmax=12, set=list(idp=2)) # idp is the power exponent
# Perform cross-valiation
logzincIDW2.cv <- gstat.cv(logzincIDW2.obj, debug.level=0, random=FALSE)
rmse.idp2 <- sqrt(mean(logzincIDW2.cv$residual^2)) # RMS prediction error
c(rmse.idp1, rmse.idp2)

# IDW interpolation with local neighborhoods
logzinc.idw = idw(log(zinc)~1, meuse, meuse.grid)
logzinc.idw.loc = idw(log(zinc)~1, meuse, meuse.grid, nmax = 16, nmin = 4, maxdist = 500)

# Plot the interpolated surfaces
spplot(logzinc.idw,"var1.pred", main = "log(zinc) IDW predictions")
spplot(logzinc.idw.loc,"var1.pred", main = "Local log(zinc) IDW predictions")


# -------------------------------------

# Kriging: Introduction
# Semivariogram Estimation and Modeling









