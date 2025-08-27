library(lubridate)
library(tidyverse)
library(here)

local_wd <- here("inst","examples","dolphin","raw")
setwd(local_wd)

df <- read_csv("bottlenoseDolph.csv")
df$id <- paste0("id",1:nrow(df))
df <- df[,c(ncol(df),1:(ncol(df)-1))]

df <- df %>%
  pivot_longer(cols = `2018-06-03`:`2020-11-25`, # Select columns to pivot
               names_to = "capture_date",         # Name of the new column for old column names
               values_to = "capture")

df$capture_date <- ymd(df$capture_date)

n <- nrow(df)

occs <- data.frame(
  capture_date=seq.Date(floor_date(df$capture_date[1], "week"), ceiling_date(df$capture_date[n], "week"), length=65),
  occ_num=1:length(occ_breaks))

tojoin <- expand.grid(id=df$id, capture_date=occs$capture_date)
tojoin <- left_join(tojoin, occs)

df <- full_join(df, tojoin)
