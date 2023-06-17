###########PROVINCIAS SPATIO-TEMPORAL
library(mapSpain)

provs <- esp_get_prov_siane()

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
       data2018$f_95_99,data2018$m_95_99,data2018$f_100,data2018$m_100, na.rm = TRUE)/sum(data2018$m_Total,data2018$f_Total, na.rm = TRUE)
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
library(dplyr)
data <- data %>%
  group_by(PROV, year) %>%
  summarize(
    Y = sum(Y),
    E = sum(E),
    SIR = mean(SIR, na.rm=TRUE)
  )
data <- data[data$PROV != "Rioja", ]
data <- data[data$PROV != "Almería", ]

library(tidyr)
dw <- data %>%
  pivot_wider(
    id_cols = PROV,
    names_from = year,
    values_from = c(Y, E, SIR)
  )

provs$iso2.prov.name.es<-ifelse(provs$iso2.prov.name.es=="Álava","Alava",provs$iso2.prov.name.es)
provs$iso2.prov.name.es<-ifelse(provs$iso2.prov.name.es=="Lérida","Lleida",provs$iso2.prov.name.es)
provs$iso2.prov.name.es<-ifelse(provs$iso2.prov.name.es=="Almería","Almeria",provs$iso2.prov.name.es)
provs$iso2.prov.name.es<-ifelse(provs$iso2.prov.name.es=="Ávila","Avila",provs$iso2.prov.name.es)
provs$iso2.prov.name.es<-ifelse(provs$iso2.prov.name.es=="Cáceres","Caceres",provs$iso2.prov.name.es)
provs$iso2.prov.name.es<-ifelse(provs$iso2.prov.name.es=="Cádiz","Cadiz",provs$iso2.prov.name.es)
provs$iso2.prov.name.es<-ifelse(provs$iso2.prov.name.es=="Castellón","Castellon",provs$iso2.prov.name.es)
provs$iso2.prov.name.es<-ifelse(provs$iso2.prov.name.es=="Ciudad Real","Ciudad-real",provs$iso2.prov.name.es)
provs$iso2.prov.name.es<-ifelse(provs$iso2.prov.name.es=="Córdoba","Cordoba",provs$iso2.prov.name.es)
provs$iso2.prov.name.es<-ifelse(provs$iso2.prov.name.es=="La Coruña","Coruna",provs$iso2.prov.name.es)
provs$iso2.prov.name.es<-ifelse(provs$iso2.prov.name.es=="Jaén","Jaen",provs$iso2.prov.name.es)
provs$iso2.prov.name.es<-ifelse(provs$iso2.prov.name.es=="León","Leon",provs$iso2.prov.name.es)
provs$iso2.prov.name.es<-ifelse(provs$iso2.prov.name.es=="Málaga","Malaga",provs$iso2.prov.name.es)
provs$iso2.prov.name.es<-ifelse(provs$iso2.prov.name.es=="Las Palmas","Las-palmas",provs$iso2.prov.name.es)
provs$iso2.prov.name.es<-ifelse(provs$iso2.prov.name.es=="Santa Cruz de Tenerife","Tenerife",provs$iso2.prov.name.es)
provs$iso2.prov.name.es<-ifelse(provs$iso2.prov.name.es=="Gerona","Girona",provs$iso2.prov.name.es)
provs$iso2.prov.name.es<-ifelse(provs$iso2.prov.name.es=="Orense","Ourense",provs$iso2.prov.name.es)
provs$iso2.prov.name.es[48]<-"Gipuzkoa"
provs$iso2.prov.name.es[49]<-"Bizkaia"


map <- merge(provs,dw, by.x = "iso2.prov.name.es", by.y = "PROV")
map_sf <- gather(map, year, SIR, c("SIR_2012","SIR_2014","SIR_2016","SIR_2018",
                                   "SIR_2020","SIR_2022"))
map_sf$year<-ifelse(map_sf$year=="SIR_2012",2012,map_sf$year)
map_sf$year<-ifelse(map_sf$year=="SIR_2014",2014,map_sf$year)
map_sf$year<-ifelse(map_sf$year=="SIR_2016",2016,map_sf$year)
map_sf$year<-ifelse(map_sf$year=="SIR_2018",2018,map_sf$year)
map_sf$year<-ifelse(map_sf$year=="SIR_2020",2020,map_sf$year)
map_sf$year<-ifelse(map_sf$year=="SIR_2022",2022,map_sf$year)
ggplot(map_sf) + geom_sf(aes(fill = SIR),size = 0.01) +
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
g <- ggplot(map_sf, aes(x = year, y = SIR, group = iso2.prov.name.es, color = iso2.prov.name.es)) +
  geom_line() +
  geom_point(size = 2) +
  theme_bw() +
  labs(color = "Province") +
  theme(legend.key.size = unit(0.2, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))
g
library(gghighlight)
g + gghighlight(iso2.prov.name.es == "Segovia")

##INLA
nb <- poly2nb(map_sf)
library(INLA)
library(spdep)
nb2INLA("map.adj", nb)
g <- inla.read.graph(filename = "map.adj")
data$idarea <- as.numeric(as.factor(data$PROV))
data$idarea1 <- data$idarea
data$idtime <- 1 + data$year - min(data$year)
data$Y<-as.integer(data$Y)
formula <- Y ~ f(idarea, model = "bym", graph = g) +
  f(idarea1, idtime, model = "iid") + idtime
res <- inla(formula,
            family = "poisson", data = data, E = E,
            control.predictor = list(compute = TRUE),
            control.compute=list(dic=TRUE, waic=TRUE)
)
map_sf$RR <- res$summary.fitted.values[, "mean"]
map_sf$LL <- res$summary.fitted.values[, "0.025quant"]
map_sf$UL <- res$summary.fitted.values[, "0.975quant"]
ggplot(map_sf) + geom_sf(aes(fill = LL),size = 0.01) +
  facet_wrap(~year, dir = "h", ncol = 3) +
  ggtitle("") + theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_gradient2(
    midpoint = 1, low = "blue", mid = "white", high = "red",
    limits = c(min(map_sf$SIR)-0.8, max(map_sf$SIR)+0.8)
  )

n <- min(nrow(data), 50) 
selected_rows <- sample(nrow(data), n, replace = FALSE)  
mse <- numeric(n)

for (i in 1:n) {
  data_loo <- data[-selected_rows[i], ]
  model_loo <- inla(formula,
                    family = "poisson", data = data_loo, E = E,
                    control.predictor = list(compute = TRUE),
                    control.compute = list(dic = TRUE, waic = TRUE))
  
  pred <- model_loo$summary.fitted.values$mean
  mse[i] <- mean((data_loo$SIR - pred)^2, na.rm = TRUE)
}

mse_inla<-mean(mse,na.rm = TRUE)

# Extract latitude and longitude values
centroids <- st_centroid(provs)
dw$longitude <- st_coordinates(centroids)[, "X"]
dw$latitude <- st_coordinates(centroids)[, "Y"]

map_sf <- gather(dw, year, SIR, c("SIR_2012","SIR_2014","SIR_2016","SIR_2018",
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
data_rf$year<-ifelse(data_rf$year=="SIR_2012",2012,data_rf$year)
data_rf$year<-ifelse(data_rf$year=="SIR_2014",2014,data_rf$year)
data_rf$year<-ifelse(data_rf$year=="SIR_2016",2016,data_rf$year)
data_rf$year<-ifelse(data_rf$year=="SIR_2018",2018,data_rf$year)
data_rf$year<-ifelse(data_rf$year=="SIR_2020",2020,data_rf$year)
data_rf$year<-ifelse(data_rf$year=="SIR_2022",2022,data_rf$year)
data_rf$SIR<-as.numeric(data_rf$SIR)
data_rf$PROV<-NULL
data_rf<-na.omit(data_rf)
n <- min(nrow(data_rf), 50) 
selected_rows <- sample(nrow(data_rf), n, replace = FALSE)  
mse <- numeric(n)
library(randomForest)
for (i in 1:n) {
  data_loo <- data_rf[-selected_rows[i], ]
  model_loo <- randomForest(SIR ~ ., mtry = 7,data = data_loo)
  pred <- predict(model_loo, newdata = data_loo)
  mse[i] <- mean((data_loo$SIR - pred)^2)
}
mse_rf<-mean(mse, na.rm=TRUE)

map_sf <- gather(map, year, SIR, c("SIR_2012","SIR_2014","SIR_2016","SIR_2018",
                                   "SIR_2020","SIR_2022"))
map_sf$year<-ifelse(map_sf$year=="SIR_2012",2012,map_sf$year)
map_sf$year<-ifelse(map_sf$year=="SIR_2014",2014,map_sf$year)
map_sf$year<-ifelse(map_sf$year=="SIR_2016",2016,map_sf$year)
map_sf$year<-ifelse(map_sf$year=="SIR_2018",2018,map_sf$year)
map_sf$year<-ifelse(map_sf$year=="SIR_2020",2020,map_sf$year)
map_sf$year<-ifelse(map_sf$year=="SIR_2022",2022,map_sf$year)
model <- randomForest(SIR ~ ., mtry = 7,data = data_rf)
map_sf$SIR_rf<- predict(model, newdata = data_rf)
ggplot(map_sf) + geom_sf(aes(fill = SIR_rf),size = 0.01) +
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

##GAM
library(mgcv)
data_gam<-map_sf
data_gam_E <- gather(data_gam, year, E, c("E_2012","E_2014","E_2016","E_2018",
                                          "E_2020","E_2022"))
data_gam_Y <- gather(data_gam_E, year, Y, c("Y_2012","Y_2014","Y_2016","Y_2018",
                                            "Y_2020","Y_2022"))
data_gam<-data_gam_Y
data_gam$year<-ifelse(data_gam$year=="Y_2012",2012,data_gam$year)
data_gam$year<-ifelse(data_gam$year=="Y_2014",2014,data_gam$year)
data_gam$year<-ifelse(data_gam$year=="Y_2016",2016,data_gam$year)
data_gam$year<-ifelse(data_gam$year=="Y_2018",2018,data_gam$year)
data_gam$year<-ifelse(data_gam$year=="Y_2020",2020,data_gam$year)
data_gam$year<-ifelse(data_gam$year=="Y_2022",2022,data_gam$year)
data_gam$SIR<-as.numeric(data_gam$SIR)
data_gam<-na.omit(data_gam)
data_gam$Y<-as.integer(data_gam$Y)
data_gam$year<-as.integer(data_gam$year)

data_gam<-na.omit(data_gam)
result <- gam(Y ~ ti(latitude)+ti(longitude)+ti(year)+ti(latitude,longitude)+
                ti(latitude,longitude,year), 
              offset = log(E), family = "poisson", data = data_gam)
n <- min(nrow(data_rf), 10) 
selected_rows <- sample(nrow(data_rf), n, replace = FALSE)  
mse <- numeric(n)
for (i in 1:n) {
  data_loo <- data_gam[-selected_rows[i], ]
  model_loo <- gam(Y ~ ti(latitude)+ti(longitude)+ti(year)+ti(latitude,longitude)+
                     ti(latitude,longitude,year), 
                   offset = log(E), family = "poisson", data = data_loo)
  pred <- predict(object=model_loo, newdata = data_loo)
  log<-ifelse(log(data_loo$SIR)<(-5),-0.0000002,log(data_loo$SIR))
  mse[i] <- mean((log - pred)^2,na.rm = TRUE)
}
mse_gam<-mean(mse, na.rm=TRUE)
