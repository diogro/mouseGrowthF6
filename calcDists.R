library(geomorph)
library(dplyr)
library(ggplot2)

data = readland.tps("/home/MouseScans/Dicom/tps/834.tps", specID = "ID")
ecd = function(x, y) sqrt(sum((x - y)^2))

markers = unlist(lapply(readLines("/home/MouseScans/Dicom/tps/marker_list.txt"), function(x) strsplit(strsplit(x, "\t")[[1]][[2]], " ")[[1]][1]))

rownames(data) = markers

plot(gpagen(data))

dists = list(
  IS_PM   = list(E = c("IS", "PM-E"), D = c("IS", "PM-D")),
  IS_NSL  = list(C = c("IS", "NSL")),
  IS_PNS  = list(C = c("IS", "PNS")),	
  PM_ZS   = list(E = c("PM-E", "ZS-E"), D = c("PM-D", "ZS-D")),	
  PM_ZI   = list(E = c("PM-E", "ZI-E"), D = c("PM-D", "ZI-D")),		
  PM_MT	  = list(E = c("PM-E", "MT-E"), D = c("PM-D", "MT-D")),	
  NSL_NA  = list(C = c("NSL", "NA")),	
  NSL_ZS  = list(E = c("NSL", "ZS-E"), D = c("NSL", "ZS-D")),	
  NSL_ZI  = list(E = c("NSL", "ZI-E"), D = c("NSL", "ZI-D")),		
  NA_BR   = list(C = c("NA", "BR")),		
  NA_PNS  = list(C = c("NA", "PNS")),		
  BR_PT	  = list(E = c("BR", "PT-E"), D = c("BR", "PT-D")),
  BR_APET = list(E = c("BR", "APET-E"), D = c("BR", "APET-D")),	
  PT_APET = list(E = c("PT-E", "APET-E"), D = c("PT-D", "APET-D")),	
  PT_BA	  = list(E = c("BA", "PT-E"), D = c("BA", "PT-D")),
  PT_EAM  = list(E = c("PT-E", "EAM-E"), D = c("PT-D", "EAM-D")),	
  PT_ZYGO = list(E = c("PT-E", "ZYGO-E"), D = c("PT-D", "ZYGO-D")),	
  PT_TSP  = list(E = c("PT-E", "TSP-E"), D = c("PT-D", "TSP-D")),	
  ZS_ZI	  = list(E = c("ZS-E", "ZI-E"), D = c("ZS-D", "ZI-D")),
  ZI_MT	  = list(E = c("ZI-E", "MT-E"), D = c("ZI-D", "MT-D")),
  ZI_ZYGO = list(E = c("ZI-E", "ZYGO-E"), D = c("ZI-D", "ZYGO-D")),	
  ZI_TSP  = list(E = c("ZI-E", "TSP-E"), D = c("ZI-D", "TSP-D")),	
  MT_PNS  = list(E = c("MT-E", "PNS"), D = c("MT-D", "PNS")),	
  PNS_APET= list(E = c("APET-E", "PNS"), D = c("APET-D", "PNS")),
  APET_BA	= list(E = c("APET-E", "BA"), D = c("APET-D", "BA")),
  APET_TS	= list(E = c("APET-E", "TS-E"), D = c("APET-D", "TS-D")),
  BA_EAM	= list(E = c("BA", "EAM-E"), D = c("BA", "EAM-D")),	
  EAM_ZYGO= list(E = c("ZYGO-E", "EAM-E"), D = c("ZYGO-D", "EAM-D")),	
  ZYGO_TSP= list(E = c("ZYGO-E", "TSP-E"), D = c("ZYGO-D", "TSP-D")),		
  LD_AS   = list(E = c("LD", "AS-E"), D = c("LD", "AS-D")),		
  BR_LD	  = list(C = c("BR", "LD")),
  OPI_LD	= list(C = c("OPI", "LD")),
  PT_AS	  = list(E = c("PT-E", "AS-E"), D = c("PT-D", "AS-D")),
  JP_AS	  = list(E = c("JP-E", "AS-E"), D = c("JP-D", "AS-D")),
  BA_OPI  = list(C = c("BA", "OPI")))

tina_distances = 10*unlist(lapply(dists, function(x) mean(unlist(lapply(x, function(m) ecd(data[m[1],,], data[m[2],,]))))))

anna_distances = as.numeric(select(readRDS("./tina_anna.rds")[[2]], IS_PM:BA_OPI)[1,])

data.frame(tina = tina_distances, ms = anna_distances) %>%
  ggplot(aes(tina, ms)) + geom_point() + geom_text(aes(label = names(dists))) + geom_abline()