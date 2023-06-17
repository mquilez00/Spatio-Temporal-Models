###ESPACIO - TEMPORAL distritos madrid
library(readxl)
data2012 <- read_excel("C:/Users/maria/OneDrive/Escritorio/DOCUMENTOS/MASTER/TFM/datos/FUNDACION MUTUALIDAD ABOGACIA/habsxcp-2012.xlsx")
data2012$year<-rep(2012,15261)
p<-sum(data2012$f_65_69,data2012$m_65_69,data2012$f_70_74,data2012$m_70_74,
       data2012$f_75_79,data2012$m_75_79,data2012$f_80_84,data2012$m_80_84,
       data2012$f_85_89,data2012$m_85_89,data2012$f_90_94,data2012$m_90_94,
       data2012$f_95_99,data2012$m_95_99,data2012$f_100,data2012$m_100)/sum(data2012$m_Total,data2012$f_Total)
data2012$Y<-data2012$f_65_69+data2012$m_65_69+data2012$f_70_74+data2012$m_70_74+
  data2012$f_75_79+data2012$m_75_79+data2012$f_80_84+data2012$m_80_84+
  data2012$f_85_89+data2012$m_85_89+data2012$f_90_94+data2012$m_90_94+
  data2012$f_95_99+data2012$m_95_99+data2012$f_100+data2012$m_100
data2012$E<-p*(data2012$m_Total+data2012$f_Total)
data2014 <- read_excel("C:/Users/maria/OneDrive/Escritorio/DOCUMENTOS/MASTER/TFM/datos/FUNDACION MUTUALIDAD ABOGACIA/habsxcp-2014.xlsx")
data2014$year<-rep(2014,15261)
p<-sum(data2014$f_65_69,data2014$m_65_69,data2014$f_70_74,data2014$m_70_74,
       data2014$f_75_79,data2014$m_75_79,data2014$f_80_84,data2014$m_80_84,
       data2014$f_85_89,data2014$m_85_89,data2014$f_90_94,data2014$m_90_94,
       data2014$f_95_99,data2014$m_95_99,data2014$f_100,data2014$m_100)/sum(data2014$m_Total,data2014$f_Total)
data2014$Y<-data2014$f_65_69+data2014$m_65_69+data2014$f_70_74+data2014$m_70_74+
  data2014$f_75_79+data2014$m_75_79+data2014$f_80_84+data2014$m_80_84+
  data2014$f_85_89+data2014$m_85_89+data2014$f_90_94+data2014$m_90_94+
  data2014$f_95_99+data2014$m_95_99+data2014$f_100+data2014$m_100
data2014$E<-p*(data2014$m_Total+data2014$f_Total)
data2016 <- read_excel("C:/Users/maria/OneDrive/Escritorio/DOCUMENTOS/MASTER/TFM/datos/FUNDACION MUTUALIDAD ABOGACIA/habsxcp-2016.xlsx")
data2016$year<-rep(2016,15261)
p<-sum(data2016$f_65_69,data2016$m_65_69,data2016$f_70_74,data2016$m_70_74,
       data2016$f_75_79,data2016$m_75_79,data2016$f_80_84,data2016$m_80_84,
       data2016$f_85_89,data2016$m_85_89,data2016$f_90_94,data2016$m_90_94,
       data2016$f_95_99,data2016$m_95_99,data2016$f_100,data2016$m_100)/sum(data2016$m_Total,data2016$f_Total)
data2016$Y<-data2016$f_65_69+data2016$m_65_69+data2016$f_70_74+data2016$m_70_74+
  data2016$f_75_79+data2016$m_75_79+data2016$f_80_84+data2016$m_80_84+
  data2016$f_85_89+data2016$m_85_89+data2016$f_90_94+data2016$m_90_94+
  data2016$f_95_99+data2016$m_95_99+data2016$f_100+data2016$m_100
data2016$E<-p*(data2016$m_Total+data2016$f_Total)
data2018 <- read_excel("C:/Users/maria/OneDrive/Escritorio/DOCUMENTOS/MASTER/TFM/datos/FUNDACION MUTUALIDAD ABOGACIA/habsxcp-2018.xlsx")
data2018$year<-rep(2018,15261)
p<-sum(data2018$f_65_69,data2018$m_65_69,data2018$f_70_74,data2018$m_70_74,
       data2018$f_75_79,data2018$m_75_79,data2018$f_80_84,data2018$m_80_84,
       data2018$f_85_89,data2018$m_85_89,data2018$f_90_94,data2018$m_90_94,
       data2018$f_95_99,data2018$m_95_99,data2018$f_100,data2018$m_100, na.rm=TRUE)/sum(data2018$m_Total,data2018$f_Total,na.rm=TRUE)
data2018$Y<-data2018$f_65_69+data2018$m_65_69+data2018$f_70_74+data2018$m_70_74+
  data2018$f_75_79+data2018$m_75_79+data2018$f_80_84+data2018$m_80_84+
  data2018$f_85_89+data2018$m_85_89+data2018$f_90_94+data2018$m_90_94+
  data2018$f_95_99+data2018$m_95_99+data2018$f_100+data2018$m_100
data2018$E<-p*(data2018$m_Total+data2018$f_Total)
data2020 <- read_excel("C:/Users/maria/OneDrive/Escritorio/DOCUMENTOS/MASTER/TFM/datos/FUNDACION MUTUALIDAD ABOGACIA/habsxcp-2020.xlsx")
data2020$year<-rep(2020,15261)
p<-sum(data2020$f_65_69,data2020$m_65_69,data2020$f_70_74,data2020$m_70_74,
       data2020$f_75_79,data2020$m_75_79,data2020$f_80_84,data2020$m_80_84,
       data2020$f_85_89,data2020$m_85_89,data2020$f_90_94,data2020$m_90_94,
       data2020$f_95_99,data2020$m_95_99,data2020$f_100,data2020$m_100)/sum(data2020$m_Total,data2020$f_Total)
data2020$Y<-data2020$f_65_69+data2020$m_65_69+data2020$f_70_74+data2020$m_70_74+
  data2020$f_75_79+data2020$m_75_79+data2020$f_80_84+data2020$m_80_84+
  data2020$f_85_89+data2020$m_85_89+data2020$f_90_94+data2020$m_90_94+
  data2020$f_95_99+data2020$m_95_99+data2020$f_100+data2020$m_100
data2020$E<-p*(data2020$m_Total+data2020$f_Total)
data2022 <- read_excel("C:/Users/maria/OneDrive/Escritorio/DOCUMENTOS/MASTER/TFM/datos/FUNDACION MUTUALIDAD ABOGACIA/habsxcp-2022.xlsx")
data2022$year<-rep(2022,15261)
p<-sum(data2022$f_65_69,data2022$m_65_69,data2022$f_70_74,data2022$m_70_74,
       data2022$f_75_79,data2022$m_75_79,data2022$f_80_84,data2022$m_80_84,
       data2022$f_85_89,data2022$m_85_89,data2022$f_90_94,data2022$m_90_94,
       data2022$f_95_99,data2022$m_95_99,data2022$f_100,data2022$m_100)/sum(data2022$m_Total,data2022$f_Total)
data2022$Y<-data2022$f_65_69+data2022$m_65_69+data2022$f_70_74+data2022$m_70_74+
  data2022$f_75_79+data2022$m_75_79+data2022$f_80_84+data2022$m_80_84+
  data2022$f_85_89+data2022$m_85_89+data2022$f_90_94+data2022$m_90_94+
  data2022$f_95_99+data2022$m_95_99+data2022$f_100+data2022$m_100
data2022$E<-p*(data2022$m_Total+data2022$f_Total)

data<-rbind(data2012,data2014,data2016,data2018,data2020,data2022)
data$SIR<-data$Y/data$E
mad<-subset(data,PROV=="Madrid")
library(dplyr)
mad <- mad %>%
  group_by(MUNICIPIO, year) %>%
  summarize(
    Y = sum(Y),
    E = sum(E),
    SIR = mean(SIR)
  )
library(tidyr)
dw <- mad %>%
  pivot_wider(
    id_cols = MUNICIPIO,
    names_from = year,
    values_from = c(Y, E, SIR)
  )

##map
library(mapSpain)
library(tidyverse)

madrid <- esp_get_munic_siane(region = "Madrid") %>%
  mutate(
    Provincia = esp_dict_translate(ine.prov.name, "es")
  )
madrid[180,6]<-"El Redegüelo"
madrid[181,6]<-"Los Baldios"
map <- merge(madrid, dw, by.x = "name", by.y = "MUNICIPIO")
library(sf)
map_sf <- st_as_sf(map)
map_sf <- gather(map, year, SIR, c("SIR_2012","SIR_2014","SIR_2016","SIR_2018",
                                   "SIR_2020","SIR_2022"))
map_sf$SIR<-ifelse(is.na(map_sf$SIR)==TRUE,mean(map_sf$SIR, na.rm=TRUE),map_sf$SIR)
map_sf$year<-ifelse(map_sf$year=="SIR_2012",2012,map_sf$year)
map_sf$year<-ifelse(map_sf$year=="SIR_2014",2014,map_sf$year)
map_sf$year<-ifelse(map_sf$year=="SIR_2016",2016,map_sf$year)
map_sf$year<-ifelse(map_sf$year=="SIR_2018",2018,map_sf$year)
map_sf$year<-ifelse(map_sf$year=="SIR_2020",2020,map_sf$year)
map_sf$year<-ifelse(map_sf$year=="SIR_2022",2022,map_sf$year)

ggplot(map_sf) + geom_sf(aes(fill = SIR),size = 0.05) +
  facet_wrap(~year, dir = "h", ncol = 3) +
  ggtitle("") + theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_gradient2(
    midpoint = 1, low = "blue", mid = "white", high = "red"
  )
g <- ggplot(map_sf, aes(x = year, y = SIR, group = name, color = name)) +
  geom_line() + geom_point(size = 2) + theme_bw() +
  guides(color = FALSE)
g

fila_max <- map_sf %>% filter(SIR == max(SIR))
library(gghighlight)
g + gghighlight(name == "Robregordo")

##INLA
nb <- poly2nb(madrid)
library(INLA)
library(spdep)
nb2INLA("map.adj", nb)
g <- inla.read.graph(filename = "map.adj")

mad$idarea <- as.numeric(as.factor(mad$MUNICIPIO))
mad$idarea1 <- mad$idarea
mad$idtime <- 1 + mad$year - min(mad$year)
mad$Y<-as.integer(mad$Y)

formula <- Y ~ f(idarea, model = "bym", graph = g) +
  f(idarea1, idtime, model = "iid") + idtime
res <- inla(formula,
            family = "poisson", data = mad, E = E,
            control.predictor = list(compute = TRUE),
            control.compute=list(dic=TRUE, waic=TRUE)
)
summary(res)

mad$RR <- res$summary.fitted.values[, "mean"]
mad$LL <- res$summary.fitted.values[, "0.025quant"]
mad$UL <- res$summary.fitted.values[, "0.975quant"]
map_sf <- merge(
  map_sf, mad,
  by.x = c("name", "year"),
  by.y = c("MUNICIPIO", "year")
)
ggplot(map_sf) + geom_sf(aes(fill = LL),size = 0.05) +
  facet_wrap(~year, dir = "h", ncol = 3) +
  ggtitle("") + theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_gradient2(
    midpoint = 1, low = "blue", mid = "white", high = "red",
    limits=c(min(map_sf$SIR.x)-0.75,max(map_sf$SIR.x)+0.75),
    name = "LL"
  )
g <- ggplot(map_sf, aes(x = year, y = LL, group = name, color = name)) +
  geom_line() + geom_point(size = 2) + theme_bw() +
  guides(color = FALSE)
g

fila_max <- map_sf %>% filter(SIR == max(SIR))
library(gghighlight)
g + gghighlight(name == "Móstoles")
#LOOCV
mad<-na.omit(mad)
n <- min(nrow(mad), 50) 
selected_rows <- sample(nrow(mad), n, replace = FALSE)  
mse <- numeric(n)

for (i in 1:n) {
  data_loo <- mad[-selected_rows[i], ]
  model_loo <- inla(formula,
                    family = "poisson", data = data_loo, E = E,
                    control.predictor = list(compute = TRUE),
                    control.compute = list(dic = TRUE, waic = TRUE))
  
  pred <- model_loo$summary.fitted.values$mean
  mse[i] <- mean((data_loo$SIR - pred)^2, na.rm = TRUE)
}

mse_inla<-mean(mse)
##RF
library(sf)
library(randomForest)
centroids <- st_centroid(madrid)

# Extract latitude and longitude values
madrid$latitude <- st_coordinates(centroids)[, "Y"]
madrid$longitude <- st_coordinates(centroids)[, "X"]

map <- merge(madrid, dw, by.x = "name", by.y = "MUNICIPIO")
map_sf <- gather(map, year, SIR, c("SIR_2012","SIR_2014","SIR_2016","SIR_2018",
                                   "SIR_2020","SIR_2022"))
data_rf<-map_sf
data_rf$Y_2012<-NULL
data_rf$Y_2014<-NULL
data_rf$Y_2016<-NULL
data_rf$Y_2018<-NULL
data_rf$Y_2020<-NULL
data_rf$Y_2022<-NULL
data_rf$E_2012<-NULL
data_rf$E_2014<-NULL
data_rf$E_2016<-NULL
data_rf$E_2018<-NULL
data_rf$E_2020<-NULL
data_rf$E_2022<-NULL
data_rf$geometry<-NULL
data_rf$year<-ifelse(data_rf$year=="SIR_2012",2012,data_rf$year)
data_rf$year<-ifelse(data_rf$year=="SIR_2014",2014,data_rf$year)
data_rf$year<-ifelse(data_rf$year=="SIR_2016",2016,data_rf$year)
data_rf$year<-ifelse(data_rf$year=="SIR_2018",2018,data_rf$year)
data_rf$year<-ifelse(data_rf$year=="SIR_2020",2020,data_rf$year)
data_rf$year<-ifelse(data_rf$year=="SIR_2022",2022,data_rf$year)
data_rf$codauto<-NULL
data_rf$name<-NULL
data_rf$ine.ccaa.name<-NULL
data_rf$ine.prov.name<-NULL
data_rf$cpro<-NULL
data_rf$cmun<-NULL
data_rf$LAU_CODE<-NULL
data_rf$Provincia<-NULL
data_rf$SIR<-as.numeric(data_rf$SIR)
data_rf<-na.omit(data_rf)
n <- nrow(data_rf)
mse <- numeric(n)
for (i in 1:n) {
  data_loo <- data_rf[-i, ]
  model_loo <- randomForest(SIR ~ ., mtry = 7,data = data_loo)
  pred <- predict(model_loo, newdata = data_loo)
  mse[i] <- mean((data_loo$SIR - pred)^2)
}
mse_rf<-mean(mse)
##GAM
library(mgcv)
data_rf<-mad
data_gam <- merge(madrid, data_rf, by.x = "name", by.y = "MUNICIPIO")
data_gam$geometry<-NULL
data_gam$name<-NULL
data_gam$codauto<-NULL
data_gam$ine.ccaa.name<-NULL
data_gam$ine.prov.name<-NULL
data_gam$cpro<-NULL
data_gam$cmun<-NULL
data_gam$LAU_CODE<-NULL
data_gam$Provincia<-NULL
library(mgcv)
data_gam$Y<-as.integer(data_gam$Y)
data_gam$year<-as.numeric(data_gam$year)
result <- gam(Y ~ ti(latitude)+ti(longitude)+ti(year)+ti(latitude,longitude)+
                ti(latitude,longitude,year), 
              offset = log(E), family = "poisson", data = data_gam)
n <- nrow(data_gam)
mse <- numeric(n)
for (i in 1:n) {
  data_loo <- data_gam[-i, ]
  model_loo <- gam(Y ~ ti(latitude)+ti(longitude)+ti(year)+ti(latitude,longitude)+
                     ti(latitude,longitude,year), 
                   offset = log(E), family = "poisson", data = data_loo)
  pred <- predict(object=model_loo, newdata = data_loo)
  mse[i] <- mean((log(data_loo$SIR) - pred)^2,na.rm=TRUE)
}
mse_gam<-mean(mse, na.rm=TRUE)
