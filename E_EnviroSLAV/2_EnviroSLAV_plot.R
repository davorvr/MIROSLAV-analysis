library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)

# Read the data
wd <- paste0(getwd(), "/")
df <- arrow::read_parquet(paste0(wd, "1_parsed_env/1_proto_env_resampled-15-minutes.parquet"))
df$ts <- as_datetime(df$ts, tz="UTC")
df$presence <- as.logical(df$presence)
df$temperature_na <- is.na(df$temperature)
df$humidity_na <- is.na(df$humidity)

x_breaks <- seq(from = as_datetime("2020-07-06 06:00:00"), to = max(df$ts), by = "48 hours")

###
# Illumination
###
## Variant with TRUE/FALSE lighting
df_illumination <- df %>% select(c("ts", "illumination"))
df_illumination$illumination <- if_else(df_illumination$illumination > 100, TRUE, FALSE)
p_lux <- 
  ggplot(df_illumination, aes(x = ts, y = NA, fill = illumination)) +
  geom_tile(aes(width = time_length("15 minutes"), height = Inf)) +
  scale_fill_manual(
    values = c("TRUE" = "yellow", "FALSE" = "black"), 
    labels = c("TRUE" = "Light detected"), 
    breaks = c("TRUE"),
    na.translate = TRUE,
    na.value = "pink",
    name = NULL
  ) +
  scale_x_datetime(breaks = x_breaks,
                   date_labels = "%Y-%m-%d",
                   expand = c(0,0)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank() # remove x labels as we'll use a common x axis for all plots
  )

## Variant with continuous scale
# ggplot(df, aes(x = ts, y = NA, fill = illumination)) +
#   geom_tile(aes(width = time_length("15 minutes"), height = Inf)) +
#   #scale_fill_gradient2(low="black", mid="yellow", high="yellow", midpoint=1000)+
#   scale_fill_gradientn(colors = c("black", "yellow", "white"),
#                        values = c(0, 0.7, 1),
#                        limits = c(0, 1500),
#                        na.value = "pink") +
#   scale_x_datetime(breaks = x_breaks,
#                    date_labels = "%Y-%m-%d",
#                    expand = c(0,0)) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     axis.title.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.title.y = element_blank(),
#     axis.text.y = element_blank(), # remove x labels as we'll use a common x axis for all plots
#     axis.ticks.y = element_blank()
#     #legend.title = element_text(size = 3),
#     #legend.text = element_text(size = 3)
#   )


###
# Person present in room
###
p_pres <- ggplot(df, aes(x = ts, y = NA, fill = presence)) +
  geom_tile(aes(width = time_length("15 minutes"), height = Inf)) +
  scale_fill_manual(
    values = c("TRUE" = "#7CAE00", "FALSE" = "#00000000"), # FALSE is a completely transparent colour
    labels = c("TRUE" = "Person detected"), 
    breaks = c("TRUE"),
    na.value = "pink",
    name = NULL
  ) +
  scale_x_datetime(breaks = x_breaks,
                   date_labels = "%Y-%m-%d",
                   expand = c(0,0)) +
  theme_minimal() +
  theme(
    #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1)
    #legend.title = element_text(size = 3), 
    #legend.text = element_text(size = 3)
  )

###
# Temperature
###
p_temp <- ggplot(df, aes(x = ts, y = temperature, colour = "Temperature")) +
  geom_line() +
  scale_color_manual(values = c("Temperature" = "#F8766D"), name=NULL) +
  scale_y_continuous(
    #breaks = c(21, 22, 23), 
    labels = function(x) paste0(x, "°C"),
    #labels = function(x) paste0(sprintf("%.1f", x), "°C"),
    limits = c(NA, 28.1),
  ) +
  scale_x_datetime(breaks = x_breaks,
                   date_labels = "%Y-%m-%d",
                   expand = c(0,0))+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank(),
        #panel.grid.minor.x = element_blank()
        )

###
# Humidity
###
p_hum <- ggplot(df, aes(x = ts, y = humidity, colour = "Relative Humidity")) +
  geom_line() +
  scale_color_manual(values = c("Relative Humidity" = "#00BFC4"), name=NULL) +
  scale_y_continuous(
    #breaks = c(40, 60, 80), 
    labels = function(x) paste0(x, "%"),
    #limits = c(min(df$humidity)-2.5, max(df$humidity)+2.5),
  ) +
  scale_x_datetime(breaks = x_breaks,
                   date_labels = "%Y-%m-%d",
                   expand = c(0,0))+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank(),
        #panel.grid.minor.x = element_blank(),
        )
        #axis.text.x = element_text(angle=45, hjust=1))

###
# Wrapped plots
###
p_wrapped <- wrap_plots(p_lux + theme(legend.key.size = unit(0.6, "lines")),
                        plot_spacer(),
                        p_temp,
                        plot_spacer(),
                        p_hum,
                        plot_spacer(),
                        p_pres + theme(legend.key.size = unit(0.6, "lines")),
                        ncol=1) +
  plot_layout(heights = c(1, -2.5, 10, -3.3, 10, -2.5, 1)) +
  plot_layout(guides = "keep") & theme(legend.position = "right", legend.box.just = "left", legend.justification = "left")
p_wrapped
ggsave(paste0(wd, "2_plot_env/2_env_plot.png"), plot=p_wrapped, width=7.55, height=3.81, dpi=600)
