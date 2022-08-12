setwd("...")


# install.packages('aqp', dep=TRUE)
# install.packages('soilDB', dep=TRUE)
# install.packages('sharpshootR', dep=TRUE)
# install.packages("rgdal", dep = TRUE)
# install.packages("raster", dep = TRUE)
# install.packages("rgeos", dep = TRUE)
# install.packages("httr", dep=TRUE)


cords = read.csv("./selected_field_coordinates.csv")
head(cords)


library(measurements)
library(dplyr)
library(stringr)
library(aqp)
library(soilDB)
library(sp)
library(rgdal)
library(plyr)
library(raster)
library(rgeos)

cords$latitude = str_trim(gsub("-|N|W|’|’’", " ", cords$latitude))
cords$longitude = str_trim(gsub("-|N|W|’|’’", " ", cords$longitude))
head(cords)



cords$new_coords = paste0(round(as.numeric(conv_unit(cords$latitude, "deg_min_sec", "dec_deg")), 6), ",-",  round(as.numeric(conv_unit(cords$longitude, "deg_min_sec", "dec_deg")), 6))


new_cords = cords %>% 
  select(Field, new_coords)

df = new_cords

df$lats = lapply(str_split(df$new_coords, ","), "[[", 1)
df$longs = lapply(str_split(df$new_coords, ","), "[[", 2)

head(df)

q = "SELECT mukey, muname
FROM mapunit
WHERE mukey IN (
SELECT * FROM SDA_Get_Mukey_from_intersection_with_WktWgs84('point(-119.763889 45.784883))"


q <- "SELECT mukey, muname
FROM mapunit
WHERE mukey IN (
SELECT * from SDA_Get_Mukey_from_intersection_with_WktWgs84('point(-121.77100 37.368402)')
)"
res = SDA_query(q)
