library(rgdal)
library(tigris)
library(sp)
library(sf)
library(tidyverse)
library(maps)



setwd("/media/ryan/hdd/documents/ICES_book_chapter")
options(tibble.width = Inf)
options(tigris_use_cache = FALSE)

#dat <- read_csv("pums.csv")

#dim(dat)
# [1] 2165680     199

# variables to keep

var_keep <- c("FIPST", "SECTOR", "N07_EMPLOYER", "RG", "TABWGT", "EMPLOYMENT_NOISY",
	      "PAYROLL_NOISY", "RECEIPTS_NOISY", "PCT1", "PCT2", "PCT3", "PCT4", "ETH1",
	      "ETH2", "ETH3", "ETH4", "RACE1", "RACE2", "RACE3", "RACE4", "SEX1", "SEX2",
	      "SEX3", "SEX4", "AGE1", "AGE2", "AGE3", "AGE4", "EDUC1", "EDUC2", "EDUC3",
	      "EDUC4")

#dat2 <- dat[, var_keep]

#dim(dat2)
#[1] 2165680      32

#write_rds(dat2, "pums_small.csv")

dat2 <- read_rds("data/pums_small.csv")

table(dat2$SECTOR)
#    11     21     22     23     31     42     44     48     51     52     53 
# 16419  13815   4007 241578 113850 109581 263545 104053  46775  94862 164211 
#    54     55     56     61     62     71     72     81     99 
#294882  11361 143025  36847 163802  74079  82757 184713   1518 
# 20 different sector

#direct estimates and sampling variances by county by sector
FIPST = sort(unique(dat2$FIPST))
SECTOR = sort(unique(dat2$SECTOR))

DIRECT <- tidyr::expand_grid(FIPST = FIPST, SECTOR = SECTOR) %>%
  dplyr::mutate(EMP = NA) %>%
  dplyr::mutate(EMP_SE = NA) %>%
  dplyr::mutate(PAY = NA) %>%
  dplyr::mutate(PAY_SE = NA) %>%
  dplyr::mutate(REC = NA) %>%
  dplyr::mutate(REC_SE = NA)

#adjustment factor for variance estimation
af <- 1.992065
for (idx1 in 1:dim(DIRECT)[1]) {
  w <- which(dat2$FIPST == DIRECT$FIPST[idx1] & dat2$SECTOR == DIRECT$SECTOR[idx1])
  tmp_dat <- dat2[w, ] %>%
    dplyr::mutate(FPC = sqrt(1 - (1/TABWGT)))
  emp_rg <- rep(0, 10)
  pay_rg <- rep(0, 10)
  rec_rg <- rep(0, 10)
  for (idx2 in 1:10) {
    w_rg <- which(tmp_dat$RG == idx2)
    emp_rg[idx2] <- sum(10*tmp_dat$TABWGT[w_rg]*tmp_dat$FPC[w_rg]*tmp_dat$EMPLOYMENT_NOISY[w_rg])
    pay_rg[idx2] <- sum(10*tmp_dat$TABWGT[w_rg]*tmp_dat$FPC[w_rg]*tmp_dat$PAYROLL_NOISY[w_rg])
    rec_rg[idx2] <- sum(10*tmp_dat$TABWGT[w_rg]*tmp_dat$FPC[w_rg]*tmp_dat$RECEIPTS_NOISY[w_rg])
  }
  DIRECT$EMP[idx1] <- sum(tmp_dat$TABWGT*tmp_dat$EMPLOYMENT_NOISY)
  DIRECT$EMP_SE[idx1] <- sqrt(sum((emp_rg - mean(emp_rg))^2/(9*10))*af)
  DIRECT$PAY[idx1] <- sum(tmp_dat$TABWGT*tmp_dat$PAYROLL_NOISY)
  DIRECT$PAY_SE[idx1] <- sqrt(sum((pay_rg - mean(pay_rg))^2/(9*10))*af)
  DIRECT$REC[idx1] <- sum(tmp_dat$TABWGT*tmp_dat$RECEIPTS_NOISY)
  DIRECT$REC_SE[idx1] <- sqrt(sum((rec_rg - mean(rec_rg))^2/(9*10))*af)
  print(idx1)
}

readr::write_rds(DIRECT, file = "data/DIRECT.rds")
DIRECT <- readr::read_rds(file = "data/DIRECT.rds")

# remove us territories
dat_sf <- tigris::states(cb = TRUE) %>%
  dplyr::filter(!(STATEFP %in% c(78, 69, 66, 60, 72)))

# note: in pums some states are combined
# S1:  Alaska (02) is combined with Wyoming (56)
# S2:  Delaware (10) is combined with DC (11)
# S3:  North Dakota (38) is combined with South Dakota (46)
# S4:  Rhode Island (44) is combined with Vermont (50)

w <- which(dat_sf$STATEFP == "02")
dat_sf$STATEFP[w] <- dat_sf$GEOID[w] <- "S1"
dat_sf$STUSPS[w] <- dat_sf$NAME[w] <- "AK/WY"

w <- which(dat_sf$STATEFP == "56")
dat_sf$STATEFP[w] <- dat_sf$GEOID[w] <- "S1"
dat_sf$STUSPS[w] <- dat_sf$NAME[w] <- "AK/WY"

w <- which(dat_sf$STATEFP == "10")
dat_sf$STATEFP[w] <- dat_sf$GEOID[w] <- "S2"
dat_sf$STUSPS[w] <- dat_sf$NAME[w] <- "DE/DC"

w <- which(dat_sf$STATEFP == "11")
dat_sf$STATEFP[w] <- dat_sf$GEOID[w] <- "S2"
dat_sf$STUSPS[w] <- dat_sf$NAME[w] <- "DE/DC"

w <- which(dat_sf$STATEFP == "38")
dat_sf$STATEFP[w] <- dat_sf$GEOID[w] <- "S3"
dat_sf$STUSPS[w] <- dat_sf$NAME[w] <- "ND/SD"

w <- which(dat_sf$STATEFP == "46")
dat_sf$STATEFP[w] <- dat_sf$GEOID[w] <- "S3"
dat_sf$STUSPS[w] <- dat_sf$NAME[w] <- "ND/SD"

w <- which(dat_sf$STATEFP == "44")
dat_sf$STATEFP[w] <- dat_sf$GEOID[w] <- "S4"
dat_sf$STUSPS[w] <- dat_sf$NAME[w] <- "RI/VT"

w <- which(dat_sf$STATEFP == "50")
dat_sf$STATEFP[w] <- dat_sf$GEOID[w] <- "S4"
dat_sf$STUSPS[w] <- dat_sf$NAME[w] <- "RI/VT"

comb_sf <- dat_sf[, c("STATEFP", "NAME")]
for (idx1 in 1:4) {
  w <- which(comb_sf$STATEFP == paste0("S", idx1))
  tmp <- st_union(comb_sf[w[1], ], comb_sf[w[2], ])
  comb_sf <- comb_sf[-w, ] 
  comb_sf <- rbind(comb_sf, tmp[, 1:2])
}

comb_sf <- dplyr::left_join(comb_sf, DIRECT, by = c("STATEFP" = "FIPST"))

ggplot(subset(comb_sf, SECTOR == "21")) +
  geom_sf(aes(fill = log(EMP))) +
  coord_sf(xlim = c(-179, -60))

hist(subset(dat2$EMPLOYMENT_NOISY, dat2$SECTOR == 11))
hist(log(subset(dat2$EMPLOYMENT_NOISY, dat2$SECTOR == 11) + 1))

hist(subset(DIRECT$EMP, DIRECT$SECTOR == 21))
hist(subset(log(DIRECT$EMP), DIRECT$SECTOR == 21))

