#setworkingdirectory

#run sample overview before and keep variables in the environment before running this code

library(ggplot2)
library(dplyr)
library(readr)
library(readxl)
library(cmocean)
library(tidyverse)
library(lubridate)
library(plotly)
library(openxlsx)
library(patchwork)
library(hrbrthemes)
library(webshot)
library(reticulate)
library(webexercises)
library(htmlwidgets)
install.packages("kaleido")

#saving the plots with names
LT_Temp_Plot <- PlotProfiles_Flexible(temperature_LT_profile, "temperature", 5, 1, "Temperature Profile", "C", RColorBrewer::brewer.pal(9, "YlOrRd"))
LT_Oxy_Plot <- PlotProfiles_Flexible(oxygen_LT_profile, "oxygen", 5, 1, "Oxygen Profile", "mg/L", rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
LT_Phs_Plot <- PlotProfiles_Flexible(phosphate_LT_profile, "phosphate", 50, 5,"Phosphate Profile", "μg/L", RColorBrewer::brewer.pal(9, "YlGnBu"))

ST_Temp_Plot <- PlotProfiles_Fixed(temperature_ST_profile, "temperature", "Temperature Profile", "°C", "YlOrRd")
ST_Oxy_Plot <- PlotProfiles_Fixed(oxygen_ST_profile, "oxygen", "Oxygen Profile", "mg/L", rev(RColorBrewer::brewer.pal(10, "RdYlBu")))
ST_Chl_Plot <- PlotProfiles_Fixed(chlorophyll_ST_profile, "chlorophyll", "Chlorophyll Profile", "μg/L", "Greens")
ST_Biotic_Plot <- biotic_plot <- plot_ly(biotic_ST_dat, x = ~date) %>%
  add_lines(y = ~bacteria, line = list(color = "#4682B4", width = 1.5), name = "Bacteria [cells/mL]") %>%
  add_trace(y = ~hnf, line = list(color = "#A0522D", width = 1.5), name = "HNF [cells/mL]", yaxis = "y2", marker = list(color = "#A0522D")) %>%
  add_markers(x = ~date, y = ~bacteria, marker = list(color = "#4682B4"), name = "Bacteria [cells/mL]") %>%
  layout(title = "Spring bloom 2021 - Epilimnion [1-8m]", xaxis = list(title = "Date")) %>%
  layout(plot_bgcolor = "white") %>%
  layout(showlegend = FALSE) %>%
  layout(yaxis = list(title = "Bacteria [million cells/mL]", titlefont = list(color = "#4682B4"), tickfont = list(color = "#4682B4"), tickvals = c(1, 2, 3, 4, 5, 6, 7, 8) * 1e6, ticktext = c("1", "2", "3", "4", "5", "6", "7", "8"), range = c(0, 8e6)),
         yaxis2 = list(title = "HNF [thousand cells/mL]", titlefont = list(color = "#A0522D"), tickfont = list(color = "#A0522D"), overlaying = "y", side = "right", tickvals = c(1, 2, 3, 4, 5, 6, 7, 8) * 1e3, ticktext = c("1", "2", "3", "4", "5", "6", "7", "8"), range = c(0, 8e3))) %>%
  layout(margin = list(l = 60, r = 80, b = 50, t = 50))

#set working directory
setwd("C:/SynologyDrive/PhD/Projects/SpringBloom2021/DataAnalysis_SpringBloom2021/Supplementary/LakeProfiles/")

#export as html

saveWidget(LT_Temp_Plot, file = "LT_Temp_Plot.html")
saveWidget(LT_Oxy_Plot, file = "LT_Oxy_Plot.html")
saveWidget(LT_Phs_Plot, file = "LT_Phs_Plot.html")
 
saveWidget(ST_Temp_Plot, file = "ST_Temp_Plot.html")
saveWidget(ST_Oxy_Plot, file = "ST_Oxy_Plot.html")
saveWidget(ST_Chl_Plot, file = "ST_Chl_Plot.html")
saveWidget(ST_Biotic_Plot, file = "ST_Biotic_Plot.html")

# orca(ST_Temp_Plot, file = "ST_Temp_Plot.png", format = "png")
# orca(ST_Chl_Plot, file = "ST_Chl_Plot.png", format = "png")
# orca(ST_Biotic_Plot, file = "ST_Biotic_Plot.png", format = "png")

export(ST_Chl_Plot, "ST_Chl_Plot.png")
export(ST_Biotic_Plot, "ST_Biotic_Plot.png")
export(ST_Temp_Plot, "ST_Temp_Plot.png")
export(ST_Oxy_Plot, "ST_Oxy_Plot.png")