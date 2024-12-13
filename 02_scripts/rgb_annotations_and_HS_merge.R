library(sf)
library(tidyverse)
library(ggplot2)
library(terra)
library(raster)
library(stars)
library(exactextractr)

#### Data importation ####

#Import data from Cloutier et al. (2024)
rgb_polygons <- read_sf("01_rawdata/trees_polygons_2023/trees_polygons.shp")

rgb_polygons_id <- c(1:nrow(rgb_polygons))

rgb_polygons[,5] <- rgb_polygons_id
colnames(rgb_polygons)[5] <- "id"

#visualize SBL trees
plot(rgb_polygons[1])

#sum(st_area(rgb_polygons_sf))

#group and summarise by species
#n_trees_per_sp <- 
#  rgb_polygons%>%
#  group_by(Label)%>%
#  summarise(count=n())

#Import HS CASI imagery from SBL flight lines 3 and 4
SBL_03_CASI <- terra::rast("01_rawdata/HS_SBL_L2G_aligned/220709_SBL-03-CA-2x1x1_v01-L2G_aligned.dat")
SBL_04_CASI <- terra::rast("01_rawdata/HS_SBL_L2G_aligned/220709_SBL-04-CA-2x1x1_v01-L2G_aligned.dat")

#Import HS SASI imagery from SBL flight lines 3 and 4
SBL_03_SASI <- terra::rast("01_rawdata/HS_SBL_L2G_aligned/220709_SBL-03-SA-1x1x1_v01-L2G_aligned.dat")
SBL_04_SASI <- terra::rast("01_rawdata/HS_SBL_L2G_aligned/220709_SBL-04-SA-1x1x1_v01-L2G_aligned.dat")

####flightlines extent####
#flightline SBL03

bbox_03 <- SBL_03_CASI%>%
  ext()%>%
  as.polygons(crs = crs(SBL_03_CASI))%>%
  st_as_sf()

#rgb polygons extent
bbox <- st_bbox(rgb_polygons)
print(bbox)

bbox_polygon <- st_as_sfc(x=bbox)

plot(SBL_03_CASI$`Channel 100 850.117nm +/- 2.391nm (rows 89- 90) (sc) 45 (850.117004)`)

band100_CASI_03 <- as.data.frame(SBL_03_CASI$`Channel 100 850.117nm +/- 2.391nm (rows 89- 90) (sc) 45 (850.117004)`, xy=T)

ggplot()+
  geom_sf(data=bbox_03, fill= "lightblue")+
  geom_raster(data=band100_CASI_03, aes(x = x, y = y, fill=band100_CASI_03$`Channel 100 850.117nm +/- 2.391nm (rows 89- 90) (sc) 45 (850.117004)`))+
  geom_sf(data=rgb_polygons, fill= "red")+
  theme_bw()

#flightline SBL04

bbox_04 <- SBL_04_CASI%>%
  ext()%>%
  as.polygons(crs = crs(SBL_04_CASI))%>%
  st_as_sf()

#rgb polygons extent
bbox <- st_bbox(rgb_polygons)
print(bbox)

bbox_polygon <- st_as_sfc(x=bbox)

plot(SBL_04_CASI$`Channel 100 850.117nm +/- 2.391nm (rows 89- 90) (sc) 45 (850.117004)`)

band100_CASI_04 <- as.data.frame(SBL_04_CASI$`Channel 100 850.117nm +/- 2.391nm (rows 89- 90) (sc) 45 (850.117004)`, xy=T)

ggplot()+
  geom_sf(data=bbox_04, fill= "lightblue")+
  geom_raster(data=band100_CASI_04, aes(x = x, y = y, fill=band100_CASI_04$`Channel 100 850.117nm +/- 2.391nm (rows 89- 90) (sc) 45 (850.117004)`))+
  geom_sf(data=rgb_polygons, fill= "red")+
  theme_bw()

#### Extract polygons from HS flight lines ####

#Band names extraction
#CASI
new_names <- gsub("(\\d+nm).*", "\\1",names(SBL_03_CASI))
new_names <-gsub("^Channel \\d+", "\\1", new_names)
new_names <- gsub("[^0-9.]", "", new_names)
CASI_bands <- as.numeric(new_names)

#SASI
new_names <- gsub("(\\d+nm).*", "\\1",names(SBL_03_SASI))
new_names <-gsub("^Channel \\d+", "\\1", new_names)
new_names <- gsub("[^0-9.]", "", new_names)
SASI_bands <- as.numeric(new_names)


names(SBL_03_CASI) <- paste0("B_", round(CASI_bands), "nm")
names(SBL_03_SASI) <- paste0("B_", round(SASI_bands), "nm")
names(SBL_04_CASI) <- paste0("B_", round(CASI_bands), "nm")
names(SBL_04_SASI) <- paste0("B_", round(SASI_bands), "nm")

### FLIGHT LINE 3

#flight line 3 visualization
plot(SBL_03_CASI$B_750nm)
plot(SBL_03_SASI$B_1752nm)

#Exactextract CASI FL 3
#CASI03_exact <- exact_extract(SBL_03_CASI, rgb_polygons,coverage_area=F, include_xy=T) %>%
  
#  setNames(unique(rgb_polygons$id)) %>%  # renaming each dataframe in the list
  
#  bind_rows(., .id = "ID") %>% # binding all rows, ID will be filled by names of dataframes
  
#  relocate(c(x,y,coverage_fraction), .after=ID) %>%
#  mutate_at(grep("B_", colnames(.)),funs(./10000))# Extract polygons from HS CASI image of SBL

#write_rds(CASI03_exact, "04_outputs/polygon_extraction_SBL_03_HS_CASI.rds")
CASI03_exact <- read_rds("04_outputs/polygon_extraction_SBL_03_HS_CASI.rds")
CASI03_exact[4] <- round(CASI03_exact[4],2)
colnames(CASI03_exact)[1:4] <- c("ID","X","Y","Weights")

CASI03_exact <- CASI03_exact %>%
  dplyr::filter(Weights>0.5)#filter weights<0.5

empty_rows <- apply(CASI03_exact[5:148], 1, function(x) all(x == 0))
empty_polygons <- CASI03_exact[empty_rows, ]
empty_polygons_CASI03 <- unique(empty_polygons$ID)#1680 trees out of this flight line

CASI03_exact_filt <-CASI03_exact[!empty_rows,]

#Averaging values
ID <- unique(CASI03_exact_filt$ID)
CASI03_mean <- CASI03_exact_filt %>% group_by(ID) %>% summarise_each(funs(mean)) #Average after grouping by ID
#20957 trees in this flight line


#Exactextract SASI FL 3
#SASI03_exact <- exact_extract(SBL_03_SASI, rgb_polygons, coverage_area=F, include_xy=T) %>%# check if you need coverage fraction or area for the filtering
  
#  setNames(unique(rgb_polygons$id)) %>%  # renaming each dataframe in the list
  
#  bind_rows(., .id = "ID") %>% # binding all rows, ID will be filled by names of dataframes
  
#  relocate(c(x,y,coverage_fraction), .after=ID)%>%
#  mutate_at(grep("B_", colnames(.)),funs(./10000))

#write_rds(SASI03_exact, "04_outputs/polygon_extraction_SBL_03_HS_SASI.rds")
SASI03_exact <- read_rds("04_outputs/polygon_extraction_SBL_03_HS_SASI.rds")
SASI03_exact[4] <- round(SASI03_exact[4],2)
colnames(SASI03_exact)[1:4] <- c("ID","X","Y","Weights")

SASI03_exact <- SASI03_exact %>%
  dplyr::filter(SASI03_exact$Weights>0.5)#filter weights<0.5

empty_rows <- apply(SASI03_exact[5:104], 1, function(x) all(x == 0))
empty_polygons <- SASI03_exact[empty_rows, ]
empty_polygons_SASI03 <- unique(empty_polygons$ID)#1567 trees out of this flight line

SASI03_exact_filt <-SASI03_exact[!empty_rows,]

#Averaging values
ID <- unique(SASI03_exact_filt$ID)
SASI03_mean <- SASI03_exact_filt %>% group_by(ID) %>% summarise_each(funs(mean)) #Average after grouping by ID
#20360 trees in this flight line


### FLIGHT LINE 4

#flight line 4 visualization
plot(SBL_04_CASI$B_750nm)
plot(SBL_04_SASI$B_1752nm)

#Exactextract CASI FL 4
#CASI04_exact <- exact_extract(SBL_04_CASI, rgb_polygons,coverage_area=F, include_xy=T) %>%
  
#  setNames(unique(rgb_polygons$id)) %>%  # renaming each dataframe in the list
  
#  bind_rows(., .id = "ID") %>% # binding all rows, ID will be filled by names of dataframes

#  relocate(c(x,y,coverage_fraction), .after=ID) %>%
#  mutate_at(grep("B_", colnames(.)),funs(./10000))# Extract polygons from HS CASI image of SBL

#write_rds(CASI04_exact, "04_outputs/polygon_extraction_SBL_04_HS_CASI.rds")
CASI04_exact <- read_rds("04_outputs/polygon_extraction_SBL_04_HS_CASI.rds")
CASI04_exact[4] <- round(CASI04_exact[4],2)
colnames(CASI04_exact)[1:4] <- c("ID","X","Y","Weights")

CASI04_exact <- CASI04_exact %>%
  dplyr::filter(CASI04_exact$Weights>0.5)#filter weights<0.5

empty_rows <- apply(CASI04_exact[5:148], 1, function(x) all(x == 0))
empty_polygons <- CASI04_exact[empty_rows, ]
empty_polygons_CASI04 <- unique(empty_polygons$ID)#17133 trees out of this flight line

CASI04_exact_filt <-CASI04_exact[!empty_rows,]

#Averaging values
ID <- unique(CASI04_exact_filt$ID)
CASI04_mean <- CASI04_exact_filt %>% group_by(ID) %>% summarise_each(funs(mean)) #Average after grouping by ID
#5899 trees in this flight line


#Exactextract SASI FL 3
#SASI04_exact <- exact_extract(SBL_04_SASI, rgb_polygons, coverage_area=F, include_xy=T) %>%# check if you need coverage fraction or area for the filtering
  
#  setNames(unique(rgb_polygons$id)) %>%  # renaming each dataframe in the list
  
#  bind_rows(., .id = "ID") %>% # binding all rows, ID will be filled by names of dataframes
  
#  relocate(c(x,y,coverage_fraction), .after=ID)%>%
#  mutate_at(grep("B_", colnames(.)),funs(./10000))

#write_rds(SASI04_exact, "04_outputs/polygon_extraction_SBL_04_HS_SASI.rds")
SASI04_exact <- read_rds("04_outputs/polygon_extraction_SBL_04_HS_SASI.rds")
SASI04_exact[4] <- round(SASI04_exact[4],2)
colnames(SASI04_exact)[1:4] <- c("ID","X","Y","Weights")

SASI04_exact <- SASI04_exact %>%
  dplyr::filter(SASI04_exact$Weights>0.5)#filter weights<0.5

empty_rows <- apply(SASI04_exact[5:104], 1, function(x) all(x == 0))
empty_polygons <- SASI04_exact[empty_rows, ]
empty_polygons_SASI04 <- unique(empty_polygons$ID)#16876 trees out of this flight line

SASI04_exact_filt <-SASI04_exact[!empty_rows,]

#Averaging values
ID <- unique(SASI04_exact_filt$ID)
SASI04_mean <- SASI04_exact_filt %>% group_by(ID) %>% summarise_each(funs(mean)) #Average after grouping by ID
#6027 trees in this flight line

row_data <- SASI04_mean[1, 5:104]


# Convert column names to numeric for x-axis
x_values <- SASI_bands
y_values <- as.numeric(row_data)

# Plot the data
plot(
  x = x_values, 
  y = y_values, 
  type = "l",                # Line plot
  xlab = "Band",             # x-axis label
  ylab = "Value",            # y-axis label
  main = "Spectral Data for Row 1" # Title
)

####Merge CASI_SASI####
CASI_SASI_03 <- inner_join(CASI03_mean,SASI03_mean, by= c("ID"))
CASI_SASI_03 <- CASI_SASI_03%>%
  relocate(c(X.y,Y.y, Weights.y), .after= Weights.x) %>%
  rename(X_CASI=X.x,Y_CASI=Y.x,Weights_CASI=Weights.x,X_SASI=X.y,Y_SASI=Y.y,Weights_SASI=Weights.y)
  
CASI_SASI_04 <- inner_join(CASI04_mean,SASI04_mean, by= c("ID"))
CASI_SASI_04 <- CASI_SASI_04%>%
  relocate(c(X.y,Y.y, Weights.y), .after= Weights.x) %>%
  rename(X_CASI=X.x,Y_CASI=Y.x,Weights_CASI=Weights.x,X_SASI=X.y,Y_SASI=Y.y,Weights_SASI=Weights.y)

#mean of duplicates

CASI_SASI_03 <- CASI_SASI_03 %>%
  rowwise() %>%
  mutate(B_1032nm.mean = mean(c(B_1032nm.x,B_1032nm.y), na.rm = TRUE)) %>%
  ungroup()%>%
  dplyr::select(-B_1032nm.x,-B_1032nm.y)%>%
  rename(B_1032nm=B_1032nm.mean)


CASI_SASI_04 <- CASI_SASI_04 %>%
  rowwise() %>%
  mutate(B_1032nm.mean = mean(c(B_1032nm.x,B_1032nm.y), na.rm = TRUE)) %>%
  ungroup()%>%
  dplyr::select(-B_1032nm.x,-B_1032nm.y)%>%
  rename(B_1032nm=B_1032nm.mean)



#### Convert data to spectra ####
colnames(CASI_SASI_03) <- gsub("nm", "", colnames(CASI_SASI_03)) #fix colnames
colnames(CASI_SASI_03) <- gsub("B_", "", colnames(CASI_SASI_03)) #fix colnames
colnames(CASI_SASI_04) <- gsub("nm", "", colnames(CASI_SASI_04)) #fix colnames
colnames(CASI_SASI_04) <- gsub("B_", "", colnames(CASI_SASI_04)) #fix colnames

#Remove x, y and weight

bands_03<-CASI_SASI_03[,-c(2:7)]
bands_04<-CASI_SASI_04[,-c(2:7)]

bands_03$ID <- as.numeric(bands_03$ID)
bands_04$ID <- as.numeric(bands_04$ID)

#write_rds(bands_03, "04_outputs/bands_only_03.rds")
#write_rds(bands_04,"04_outputs/bands_only_04.rds")

bands_03 <- read_rds("04_outputs/bands_only_03.rds")
bands_04 <- read_rds("04_outputs/bands_only_04.rds")
library(spectrolab)

bands_03 <- bands_03[ , order(as.numeric(colnames(bands_03)))]
which(is.na(bands_03)==T) #No NAs introduced
bands_04 <- bands_04[ , order(as.numeric(colnames(bands_04)))]
which(is.na(bands_04)==T) #No NAs introduced

#merge flight lines with HS data

bands_0304 <- bind_rows(bands_03,bands_04)
length(which(duplicated(bands_0304$ID))) #4012 duplicates

compare_lines <- bands_0304[bands_0304$ID==11258,]

#keep highest reflectance values for duplicates
bands_0304_filtered <- bands_0304 %>%
  mutate(row_mean = rowMeans(dplyr::select(., 1:243), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(ID) %>%
  dplyr::filter(row_mean>=0.01)%>%
  slice_max(order_by = row_mean, n = 1) %>%
  dplyr::select(-row_mean)%>%
  ungroup()

length(unique(bands_0304_filtered$ID))

#MUST ALSO FILTER OUT SHADOWED PIXELS AND VERY LOW REFLECTANCE PIXELS IN GENERAL

#as_spectra
spectra_sbl <- as_spectra(bands_0304_filtered, name_idx = 244)
spectral_data <- as.data.frame(spectra_sbl)


# Add the ID column to the metadata of the spectra object


# Match the reflectance data of CASI and SASI wavelengths
splice_bands_from_manufacturer <- 957.50 #first overlapping wavelengths from CASI/SASI

sbl_reflect_matched <-  match_sensors(x = spectra_sbl, 
                                        splice_at = splice_bands_from_manufacturer,
                                        fixed_sensor=2,
                                        interpolate_wvl = c(11, 1))
# check if correction worked
par(mfrow=c(1,2))
plot(spectra_sbl, lwd = 1.2, xlab="Wavelengths in nm",ylab="Reflectance",
     col = rgb(1, 0, 0))

plot_regions(spectra_sbl[4], regions = default_spec_regions(), add = TRUE)

plot(sbl_reflect_matched, lwd = 1.2, xlab="Wavelengths in nm",ylab="Reflectance",
     col = rgb(1, 0, 0))
plot_regions(sbl_reflect_matched, regions = default_spec_regions(), add = TRUE)

plot(spectra_sbl[6,])
plot(sbl_reflect_matched[6,])


#### 2) Transfer to hsdar for further processing ####
library(hsdar)
library(hyperSpec)

Spectrum_sbl <- as.matrix(sbl_reflect_matched) %>% data.frame()
Wavelengths_sbl <- gsub("X", "", Spectrum_sbl %>% colnames() %>% as.character()) %>% as.numeric()

speclib_sbl <- speclib(spectra = t(Spectrum_sbl), wavelength = Wavelengths_sbl) # load Bands into spectral library
idSpeclib(speclib_sbl) <- paste0(rownames(Spectrum_sbl)) # add IDs

speclib_sbl@transformation

### 3) Pre-processing of spectra in hsdar ####

band_mask <- data.frame(lb = c(350,920),#Removing water absorption
                        up = c(400,957))
speclib_sbl_mask <- speclib_sbl
hsdar::mask(speclib_sbl_mask) <- band_mask # mask absorption bands and noisy bands (beginnning and end)
sbl_interpolated <- interpolate.mask(speclib_sbl_mask) # interpolate masked regions
plot(speclib_sbl)
plot(sbl_interpolated)

sbl_smooth <- noiseFiltering(speclib_sbl, method = 'sgolay', n = 21, p=4) # smoothing of spectra with window of 21 bands and polynomial 4th (Wallis et al., YYYY)
plot(sbl_smooth)# p = polynomial order w = window size (must be odd) m = m-th derivative 
hsdar::spectra(sbl_smooth)[hsdar::spectra(sbl_smooth) < 0] <- 0 #set negative values to 0

#CR_sbl_bd <- transformSpeclib(sbl_smooth, out="bd", method = "sh") # calculating continuum removal: band depth (closest related to original spectrum)
#write_rds(CR_sbl_bd, "04_outputs/cont_rem_sbl_data.rds")
CR_sbl_bd <- read_rds("04_outputs/cont_rem_sbl_data.rds")
hsdar::spectra(CR_sbl_bd)[hsdar::spectra(CR_sbl_bd) < 0] <- 0 #set negative values to 0
hsdar::spectra(CR_sbl_bd)[hsdar::spectra(CR_sbl_bd) > 1] <- 1 #set values higher 1 to 1

#closest wavelengths
target_wavelengths <- c(472,601,1002,1498,2008)#472, 601, 1002,1498,2008  were removed because weird things happened
available_wavelengths <- CR_sbl_bd@wavelength
closest_wavelengths <- sapply(target_wavelengths, function(w) {
  available_wavelengths[which.min(abs(available_wavelengths - w))]
})
closest_wavelengths

featureSelection <- specfeat(CR_sbl_bd, target_wavelengths, tol = 5) #Example to isolate the absorption features around 450, 600, 1500 and 2000 nm.##Continuum removal

# Plot pre-processing of spectra ####
library(dlfUtils)

par(mfrow=c(2,2))
plot(speclib_sbl, FUN=max, new=T, col="red", main="a: Raw spectra (min/mean/max)")
plot(speclib_sbl, FUN=min, new=F, col="blue")
plot(speclib_sbl, FUN=mean, new=F)

plot(speclib_sbl_mask, FUN=max,  col="red", main="b: Atmospheric absorption regions masked")
plot(speclib_sbl_mask, FUN=mean, new=F,)
plot(speclib_sbl_mask, FUN=min, new=F, col="blue")
text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
     line2user(line=8, side=3), '\nSpectra of SBL trees', xpd=NA, cex=2, font=2)

plot(sbl_interpolated, FUN=max, col="red", main="c: Spectra with masked regions interpolated")
plot(sbl_interpolated,  new=F, FUN=mean)
plot(sbl_interpolated, FUN=min, col="blue", new=F)
plot(sbl_smooth, FUN=max, col="red", main="d: Smoothed spectra (min/mean/max)")
plot(sbl_smooth, FUN=min, col="blue",new=F)
plot(sbl_smooth, FUN=mean, col="black", new=F)

plot(CR_sbl_bd, main="e: Continuum removal - Band depth segment hull",FUN="max", col="red")
plot(CR_sbl_bd, FUN="mean", col="black", new=F)
plot(CR_sbl_bd, FUN="min",  new=F, col="blue")

plot(featureSelection, main="Absorption features",fnumber=1, col="violet")
plot(featureSelection, fnumber=2, col="orange", new=F)
plot(featureSelection, fnumber=3, col="red", new=F)
plot(featureSelection, fnumber=4, col="green", new=F)
plot(featureSelection, fnumber=5, col="blue", new=F)
# some adjustment needed for feature 3-5 
# (for single plots - could be done within the continuum removal by changing the segment hull points)
# we can use e.g. the area under the single absorption features to characterize the plots and correlate these areas against elevation

#### 4) Extract processed spectra to table with ID and species ####
SI_sbl <- (sbl_smooth@ID) # extract SI
SI_sbl <- as.data.frame(SI_sbl)
sbl_BD <- sbl_smooth@spectra@spectra_ma %>% data.frame() # extract spectral bd
colnames(sbl_BD) <- paste0("BD_Band", as.integer(sbl_smooth@wavelength), "nm") # define colnames of Bands, now Wavelenghts! not band no.
sbl_BD %>% 
  select_if(~ !is.numeric(.) || sum(.) == "NaN"|| sum(.) == 0) #no NA data
sblID_bd <- cbind(SI_sbl, sbl_BD)# dataframe of CR and Metadata for further processing
colnames(sblID_bd)[1] <- "ID"
colnames(rgb_polygons)[5] <- "ID"
rgb_polygons$ID <- as.character(rgb_polygons$ID)
sblID_bd$ID <- gsub("[^0-9.-]", "", sblID_bd$ID)

sbl_sp_ID_bd <- merge(sblID_bd,rgb_polygons,by="ID")

#write_rds(sbl_sp_ID_bd, "04_outputs/sbl_smoothed_spectra.rds")
sbl_sp_ID_bd <- read_rds("04_outputs/sbl_smoothed_spectra.rds")

write_sf(sbl_sp_ID_bd, "04_outputs/sbl_smoothed_spectra.gpkg")
write.csv(x=sbl_sp_ID_bd,file="04_outputs/sbl_smoothed_spectra.csv", row.names = F, sep = ";")

#### 5) Extract processed CR spectra to table with ID and species ####
SI_sblCR <- (CR_sbl_bd@ID) # extract SI
SI_sblCR <- as.data.frame(SI_sblCR)
sbl_BD_CR <- CR_sbl_bd@spectra@spectra_ma %>% data.frame() # extract spectral bd
colnames(sbl_BD_CR) <- paste0("BD_Band", as.integer(CR_sbl_bd@wavelength), "nm") # define colnames of Bands, now Wavelenghts! not band no.
sbl_BD_CR %>% 
  select_if(~ !is.numeric(.) || is.nan(sum(.)) || sum(.) == 0)%>%
  mutate(across(everything(), ~ ifelse(is.nan(.), 0, .)))
sblID_bd_CR <- cbind(SI_sblCR, sbl_BD_CR)# dataframe of CR and Metadata for further processing
colnames(sblID_bd_CR)[1] <- "ID"
colnames(rgb_polygons)[5] <- "ID"
sblID_bd_CR$ID <- gsub("[^0-9.-]", "", sblID_bd_CR$ID)

#sbl_sp_ID_bd_CR <- merge(sblID_bd_CR,rgb_polygons,by="ID") #matching species and area with HS data

#write_rds(sbl_sp_ID_bd_CR,"04_outputs/sbl_CR_spectra.rds")

sbl_sp_ID_bd_CR <- read_rds("04_outputs/sbl_CR_spectra.rds")

#### 6) ggplot by species with smoothed spectra####
library(reshape2)
library(colorspace)

melted_sbl <- melt(sbl_sp_ID_bd%>%dplyr::select(-geometry), id.vars = c("Label", "ID"), variable.name = "Wavenumber", value.name = "Reflectance")
melted_sbl$Wavenumber<-gsub("nm", "", melted_sbl$Wavenumber) #fix colnames
melted_sbl$Wavenumber <- gsub("BD_Band", "", melted_sbl$Wavenumber) #fix colnames
class(melted_sbl$Wavenumber)
melted_sbl$Wavenumber <- as.numeric(melted_sbl$Wavenumber)

melted_sbl_by_sp <- melted_sbl%>%
    group_by(Label,Wavenumber)%>%
    summarize(mean_reflectance=mean(Reflectance), sd_reflectance = sd(Reflectance))

# Generate 29 distinct colors
my_colors <- qualitative_hcl(length(unique(melted_sbl_by_sp$Label)), palette = "Dark3")

ggplot(melted_sbl_by_sp, aes(x=Wavenumber, y=mean_reflectance, color = Label, fill=Label)) + 
  geom_line(size=0.75) +
  #geom_ribbon(aes(ymin=mean_reflectance-sd_reflectance, ymax=mean_reflectance+sd_reflectance), alpha=0.07, color=NA)+
  scale_color_manual(values = my_colors)+
  theme_bw() + 
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  guides(color = guide_legend(nrow = 2)) + # Split legend into 2 rows
  theme(plot.margin=unit(c(1,1,0.5,1),"cm")) +
  labs(x="Longueurs d'onde (nm)", y="Reflectance")+
  ylim(0,0.62)+
  ggtitle("Figure 1 : Reflectance moyenne par espece a la SBL")+
  theme(plot.title = element_text(hjust = 0.5, size=16))

ggsave(filename="03_figures/smoothed_ref_all_SBL.jpg")

  #test avec 5 especes
melted_ACSAFAGR <- melted_sbl_by_sp %>%
  dplyr::filter(Label %in% c("ACSA","BEPA","FAGR","ABBA","THOC"))

ggplot(melted_ACSAFAGR, aes(x=Wavenumber, y=mean_reflectance, group = Label)) + 
  geom_line(size=0.75, aes(color=Label)) +
  geom_ribbon(aes(ymin=mean_reflectance-sd_reflectance, ymax=mean_reflectance+sd_reflectance, fill=Label), alpha=0.1)+
  theme_bw() + 
  theme(legend.position = "bottom") + 
  theme(plot.margin=unit(c(1,1,0.5,1),"cm")) +
  labs(x="Longueurs d'onde (nm)", y="Reflectance")+
  ylim(0,0.62)+
  ggtitle("Figure 1 : Reflectance moyenne par espece a la SBL")+
  theme(plot.title = element_text(hjust = 0.5, size=16))
ggsave(filename="03_figures/smooth_ref_5sp_SBL.jpg")


length(which(sbl_sp_ID_bd$Label=="Betula")) #only 1, which explains great variability

#### 6) ggplot by species with CR spectra####

melted_sbl_CR <- melt(sbl_sp_ID_bd_CR%>%dplyr::select(-geometry), id.vars = c("Label", "ID"), variable.name = "Wavenumber", value.name = "Reflectance")
melted_sbl_CR$Wavenumber<-gsub("nm", "", melted_sbl_CR$Wavenumber) #fix colnames
melted_sbl_CR$Wavenumber <- gsub("BD_Band", "", melted_sbl_CR$Wavenumber) #fix colnames
class(melted_sbl_CR$Wavenumber)
melted_sbl_CR$Wavenumber <- as.numeric(melted_sbl_CR$Wavenumber)

melted_sbl_CR_by_sp <- melted_sbl_CR%>%
  group_by(Label,Wavenumber)%>%
  summarize(mean_reflectance=mean(Reflectance), sd_reflectance = sd(Reflectance))

ggplot(melted_sbl_CR_by_sp, aes(x=Wavenumber, y=mean_reflectance, color = Label, fill=Label)) + 
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin=mean_reflectance-sd_reflectance, ymax=mean_reflectance+sd_reflectance), alpha=0.07, color=NA)+
  scale_color_manual(values = my_colors)+
  theme_bw() + 
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  guides(color = guide_legend(nrow = 2)) + # Split legend into 2 rows
  theme(plot.margin=unit(c(1,1,0.5,1),"cm")) +
  labs(x="Longueurs d'onde (nm)", y="Reflectance")+
  ylim(0,1)+
  ggtitle("Figure 1 : Reflectance moyenne par espece a la SBL")+
  theme(plot.title = element_text(hjust = 0.5, size=16))

ggsave(filename="03_figures/CR_ref_SBL.jpg")


#test avec 5 especes
melted_ACSAFAGR_CR <- melted_sbl_CR_by_sp %>%
  dplyr::filter(Label %in% c("ACSA","ACRU","FAGR","ABBA", "THOC"))

ggplot(melted_ACSAFAGR_CR, aes(x=Wavenumber, y=mean_reflectance, group = Label)) + 
  geom_line(size=0.75, aes(color=Label)) +
  geom_ribbon(aes(ymin=mean_reflectance-sd_reflectance, ymax=mean_reflectance+sd_reflectance, fill=Label), alpha=0.1)+
  theme_bw() + 
  theme(legend.position = "bottom") + 
  theme(plot.margin=unit(c(1,1,0.5,1),"cm")) +
  labs(x="Longueurs d'onde (nm)", y="Reflectance")+
  ylim(0,1)+
  ggtitle("Figure 1 : Mean CR reflectance by species at SBL")+
  theme(plot.title = element_text(hjust = 0.5, size=10),legend.text=element_text(size=rel(0.5)),legend.title = element_text(size=10), axis.title = element_text(size=10))

ggsave(filename="03_figures/CR_ref_5sp_SBL.jpg")

#### 8) Test PLSDA with smoothed spectra ####
library(mdatools)
library(tidyverse)
library(smotefamily)
library(sf)

set.seed(133) #For reproducibility
#import
sbl_sp_ID_bd_full <- read_rds("04_outputs/sbl_smoothed_spectra.rds")

sbl_sp_ID_bd_full$Label <- as.factor(sbl_sp_ID_bd_full$Label)

tapply(sbl_sp_ID_bd_full$Label,sbl_sp_ID_bd_full$Label, FUN=length)#length of every label class


#Global model

#Filter 2 most present species only
sbl_sp_ID_bd <- sbl_sp_ID_bd_full[,-247]%>%
  dplyr::group_by(Label) %>%  # Group by 'Label'
  dplyr::filter(Label %in% c("BEPA","ACRU")) %>%
  dplyr::ungroup() %>%   
  dplyr::filter(!is.na(Label))


#Filter for only species with more than 300 occurences
#sbl_sp_ID_bd <- sbl_sp_ID_bd_full[,-247]%>%
#  dplyr::group_by(Label) %>%  # Group by 'Label'
#  dplyr::filter(n() > 300) %>% # Keep groups with more than 300 occurrences
#  dplyr::filter(!Label %in% c("Acer","Picea")) %>%
#  dplyr::ungroup() %>%   
#  dplyr::filter(!is.na(Label))

sbl_sp_ID_bd$Label <- droplevels(sbl_sp_ID_bd$Label)#Drop unused factor levels from filtered out species
sbl_sp_ID_bd$Label <- as.factor(sbl_sp_ID_bd$Label)
levels(sbl_sp_ID_bd$Label)

table(sbl_sp_ID_bd[, 244])

nrow(sbl_sp_ID_bd)/nrow(sbl_sp_ID_bd_full)

# Apply SMOTE (maybe not for large class numbers)
#smote_sbl_total <- SMOTE(X = sbl_sp_ID_bd[,2:244], target = sbl_sp_ID_bd$Label, 
#                      K = 5, dup_size = 3)

# Combine the SMOTE result into a new data frame
smote_sbl <- sbl_sp_ID_bd[,2:244]

#smote_sbl <- data.frame(smote_sbl_total$data)
names(smote_sbl)[ncol(smote_sbl)] <- "Species"
smote_sbl$Species <- as.factor(smote_sbl$Species)

smote_sbl$Species <- droplevels(smote_sbl$Species)#Drop unused factor levels from filtered out species


# Check the distribution of the target variable after SMOTE
table(smote_sbl$Species)

cal.ind <- smote_sbl[sample(1:nrow(smote_sbl),0.80*nrow(smote_sbl), replace=F),]
length(which(cal.ind$Species=="ACRU"))

val.ind <- anti_join(smote_sbl,cal.ind)

table(val.ind$Species)

Xc <-  smote_sbl[rownames(cal.ind), 1:242]
Xv <-  smote_sbl[rownames(val.ind), 1:242]

cc.all <-  cal.ind[rownames(cal.ind), 243]
cc.all <- as.factor (cc.all$Species)
cv.all <-  val.ind[rownames(val.ind), 243]
cv.all <- as.factor(cv.all$Species)
table(cc.all)
table(cv.all)

cc.abba <- cc.all == "ABBA"
cv.abba <- cv.all == "ABBA"


m.all <-  mdatools::plsda(Xc, cc.all, ncomp = 15, cv = 1)
m.abba <- mdatools::plsda(Xc, cc.abba, ncomp=2, cv = 1, classname = "ABBA")

summary(m.all)
summary(m.abba)

getConfusionMatrix(m.all$calres)
getConfusionMatrix(m.abba$calres)

#Classification plots
par(mfrow = c(1, 2))
plotPredictions(m.all) #choose which class to show predictions for with nc= n
plotPredictions(m.abba) #choose which class to show predictions for with nc= n

#Performance plots
par(mfrow = c(3, 2))
plotMisclassified(m.all, nc = 3)
plotMisclassified(m.abba)
plotSensitivity(m.all, nc = 3)
plotSensitivity(m.abba)
plotSpecificity(m.all, nc = 3)
plotSpecificity(m.abba)

par(mfrow = c(1, 2))
plotRegcoeffs(m.all, ncomp = 3, ny = 3)
plotRegcoeffs(m.abba, ncomp = 1, show.ci = TRUE)

#Predictions for new data (validation subset)

res <- predict(m.all,Xv,cv.all)
summary(res)

plotPredictions(res)

resabba <-  predict(m.abba, Xv, cv.all)
summary(resabba)

plotPredictions(resabba)

####SUBSET WITH ABBA,BEPA and ACRU####
sbl_sp_ID_bd_sub <- subset(sbl_sp_ID_bd, Label==c("ABBA","BEPA","ACRU"))

sbl_sp_ID_bd_sub$Label <- as.factor (sbl_sp_ID_bd_sub$Label)
sbl_sp_ID_bd_sub$Label <- droplevels(sbl_sp_ID_bd_sub$Label)#Drop unused factor levels from filtered out species

sbl_sp_ID_bd_sub[, 245] <- factor(sbl_sp_ID_bd_sub[, 245])
table(sbl_sp_ID_bd_sub[, 245])

# Apply SMOTE
smote_sbl_sub_total <- SMOTE(X = sbl_sp_ID_bd_sub[,2:244], target = sbl_sp_ID_bd_sub$Label, 
                         K = 5, dup_size = 1)

# Combine the SMOTE result into a new data frame
smote_sbl_sub <- data.frame(smote_sbl_sub_total$data)
names(smote_sbl_sub)[ncol(smote_sbl_sub)] <- "Species"
smote_sbl_sub$Species <- as.factor(smote_sbl_sub$Species)

# Check the distribution of the target variable after SMOTE
table(smote_sbl_sub$Species)

cal_ind_sub <- smote_sbl_sub[sample(1:nrow(smote_sbl_sub),0.80*nrow(smote_sbl_sub), replace=F),]
length(which(cal_ind_sub$Species=="ACRU"))

val_ind_sub <- anti_join(smote_sbl_sub,cal_ind_sub)

Xc_sub = smote_sbl_sub[rownames(cal_ind_sub), 1:243]
Xv_sub = smote_sbl_sub[rownames(val_ind_sub), 1:243]

cc_all_sub = smote_sbl_sub[rownames(cal_ind_sub), 244]
cv_all_sub = smote_sbl_sub[rownames(val_ind_sub), 244]

cc_abba_sub <- as.factor(cc_all_sub=="ABBA") %>%
  droplevels()
cv_abba_sub = cv_all_sub == "ABBA"

str(cc_abba_sub)

m_all_sub = plsda(Xc_sub, cc_all_sub, ncomp = 15, cv = 1)
m_abba_sub = plsda(Xc_sub, cc_abba_sub, 5, cv = 1, classname = "ABBA")

summary(m_all_sub)
summary(m_abba_sub)

# Cross-validation to evaluate performance
set.seed(123)

# Extract classification accuracy
misclass <- m_all_sub$cvres$misclassified
accuracy <- 1-misclass

# Find the optimal number of components
optimal_comp <- which.max(accuracy)
cat("Optimal number of components:", optimal_comp, "\n")

#Classification plots
par(mfrow = c(1, 2))
plotPredictions(m_all_sub, nc=1) #choose which class to show predictions for with nc= n
plotPredictions(m_abba_sub, nc=1) #choose which class to show predictions for with nc= n

#Performance plots
par(mfrow = c(3, 2))
plotMisclassified(m_all_sub, nc = 3)
plotMisclassified(m_abba_sub)
plotSensitivity(m_all_sub, nc = 3)
plotSensitivity(m_abba_sub)
plotSpecificity(m_all_sub, nc = 3)
plotSpecificity(m_abba_sub)

par(mfrow = c(1, 2))
plotRegcoeffs(m_all_sub, ncomp = 3, ny = 3)
plotRegcoeffs(m_abba_sub, ncomp = 1, show.ci = TRUE)

#Predictions for new data (validation subset)

res_sub <- predict(m_all_sub,Xv_sub,cv_all_sub, type="prob")
summary(res_sub)

resabba_sub <-  predict(m_abba_sub, Xv_sub, cv_all_sub, )
summary(resabba_sub)


####TEST TOUTES SP IMP SBL####
set.seed(133) #For reproducibility
#import
sbl_sp_ID_bd <- read_rds("04_outputs/sbl_smoothed_spectra.rds")

sbl_sp_ID_bd$Label <- as.factor(sbl_sp_ID_bd$Label)

(unique(sbl_sp_ID_bd$Label))

#Global model
sbl_sp_ID_bd_2<- subset(sbl_sp_ID_bd, Label%in%c("ABBA","BEPA","ACRU", "ACSA", "ACPE", "BEAL", "FAGR", "LALA", "mort","PIMA", "PIST", "POGR", "THOC"))


sbl_sp_ID_bd <- sbl_sp_ID_bd[,-248]%>%
  dplyr::filter(n()>10)%>%
  dplyr::filter(!is.na(Label))

sbl_sp_ID_bd$Label <- droplevels(sbl_sp_ID_bd$Label)#Drop unused factor levels from filtered out species
sbl_sp_ID_bd$Label <- as.factor(sbl_sp_ID_bd$Label)
levels(sbl_sp_ID_bd$Label)

table(sbl_sp_ID_bd[, 245])

# Apply SMOTE (maybe not for large class numbers)
#smote_sbl_total <- SMOTE(X = sbl_sp_ID_bd[,2:244], target = sbl_sp_ID_bd$Label, 
#                      K = 5, dup_size = 3)

# Combine the SMOTE result into a new data frame
smote_sbl <- sbl_sp_ID_bd[,2:245]

#smote_sbl <- data.frame(smote_sbl_total$data)
names(smote_sbl)[ncol(smote_sbl)] <- "Species"
smote_sbl$Species <- as.factor(smote_sbl$Species)

smote_sbl$Species <- droplevels(smote_sbl$Species)#Drop unused factor levels from filtered out species


# Check the distribution of the target variable after SMOTE
table(smote_sbl$Species)

cal.ind <- smote_sbl[sample(1:nrow(smote_sbl),0.70*nrow(smote_sbl), replace=F),]
length(which(cal.ind$Species=="FAGR"))

val.ind <- anti_join(smote_sbl,cal.ind)

table(val.ind$Species)

Xc <-  smote_sbl[rownames(cal.ind), 1:243]
Xv <-  smote_sbl[rownames(val.ind), 1:243]

cc.all <-  cal.ind[rownames(cal.ind), 244]
cv.all <-  val.ind[rownames(val.ind), 244]
table(cc.all)
table(cv.all)

cc.abba <- cc.all == "ABBA"
cv.abba <- cv.all == "ABBA"

m.all <-  plsda(Xc, cc.all, 20, cv = 1)
m.abba <- plsda(Xc, cc.abba, 20, cv = 1, classname = "ABBA")

summary(m.all)
summary(m.abba)

getConfusionMatrix(m.all$calres)
getConfusionMatrix(m.abba$calres)

#Classification plots
par(mfrow = c(1, 2))
plotPredictions(m.all) #choose which class to show predictions for with nc= n
plotPredictions(m.abba) #choose which class to show predictions for with nc= n

#Performance plots
par(mfrow = c(3, 2))
plotMisclassified(m.all, nc = 3)
plotMisclassified(m.abba)
plotSensitivity(m.all, nc = 3)
plotSensitivity(m.abba)
plotSpecificity(m.all, nc = 3)
plotSpecificity(m.abba)

par(mfrow = c(1, 2))
plotRegcoeffs(m.all, ncomp = 3, ny = 3)
plotRegcoeffs(m.abba, ncomp = 1, show.ci = TRUE)

#Predictions for new data (validation subset)

res <- predict(m.all,Xv,cv.all)
summary(res)

plotPredictions(res)

resabba <-  predict(m.abba, Xv, cv.all)
summary(resabba)

plotPredictions(resabba)

