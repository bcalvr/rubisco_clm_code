library(readxl)
library(dplyr)
library(tictoc)
library(data.table)
library(ggplot2)
library(viridis)
library(glue)
library(stringr)
theme_update(plot.title = element_text(hjust = 0.5))
library(sp)
library(gstat)
library(grid)
library(RColorBrewer)
library(here)

args <- commandArgs(TRUE)
project <- args[1]
dateStamp <- args[2]
note <- args[3]
print(glue("{project} run from {dateStamp} ({note}) data beginning"))
print(getwd())
tic(glue("{project} run from {dateStamp} ({note}) data"))
tic("Data and functions import")
setwd('..')
print(getwd())
project_dir <- getwd()
setwd("code/")
source("BCC_kriging_functions.r")
if (!file.exists(glue("{project_dir}/results/{dateStamp}-{note}-data"))) {
    dir.create(glue("{project_dir}/results/{dateStamp}-{note}-data"))
}
setwd(glue("{project_dir}/results/{dateStamp}-{note}-data"))
y_vars <- read.table("y_vars.txt")
z_vars <- read.table("z_vars.txt")
geneLength <- read.table("xLength.txt")[1, 1]
xx <- read.table("xName.txt")[1, 1]
df <- read.csv("xyz_data.csv")
toc()
for (yy in y_vars) {
    for (zz in z_vars) {
        if (yy != zz) {
            print(glue("Now starting {project}, {note} run. xvar = {xx}, y = {yy}, z = {zz}"))
            tic(glue("{project}, {note} run. xvar = {xx}, y = {yy}, z = {zz} finished"))
            filename <- glue("{dateStamp}-{xx}-{yy}-{zz}-{note}")
            export_dir <- glue("{project_dir}/results/{filename}")
            choices <- c(xx, yy, zz)
            cols <- data.frame(choices, c("x", "y", "z"))

            df1 <- df[cols[, 1]]
            colnames(df1) <- cols[, 2]
            df1$y <- as.numeric(df1$y)
            df1$z <- as.numeric(df1$z)
            df1 <- df1[complete.cases(df1), ]

            df1$y <- (df1$y - min(df1$y)) / (max(df1$y) - min(df1$y))

            dt <- data.table(df1)
            dt <- dt[, mean(z), by = c("x", "y")]
            colnames(dt) <- cols[, 2]

            xVar <- cols[1, 1]
            yVar <- cols[2, 1]
            zVar <- cols[3, 1]
            krigedata <- dt
            p_thresh <- -.5

            nMin <- 5
            nMax <- 50
            setwd(project_dir)
            maxd <- max(dist(krigedata[, 1:2]))
            coordinates(krigedata) <- c("x", "y")

            tic("Prepare for kriging")
            krigeout <- kriging_process(krigedata, notes = note, p = p_thresh, r_dir = project_dir, ex_dir = export_dir)
            pearson <- krigeout$out1
            maxd.fraction <- krigeout$out2
            vmf <- krigeout$out3
            lin_pearson <- krigeout$out4
            rm(krigeout)
            print(format(Sys.time(), "%Y%m%d-%H:%M:%S"))
            toc()

            tic("Kriging complete")
            krigemap <- do_krige(krigedata, vmf)
            print(format(Sys.time(), "%Y%m%d-%H:%M:%S"))
            toc()

            if (pearson >= p_thresh & pearson**2 > (lin_pearson**2 + 0)) {
                #  Create kriging landscapes
                tic("Landscapes")
                krigeDF <- kriging_landscapes(krigemap, krigedata, p = pearson, plots = TRUE)
                print(format(Sys.time(), "%Y%m%d-%H:%M:%S"))
                toc()

                # Weighted predictions for each residue
                tic("IVW")
                wm <- ivw(krigeDF, geneLength)
                barcode_plotter(wm, zVar, yVar, "ivw", plots = TRUE, show = TRUE)
                print(format(Sys.time(), "%Y%m%d-%H:%M:%S"))
                print(format(Sys.time(), "%Y%m%d-%H:%M:%S"))
                toc()

                # Structure mapping
                tic("Structure mapping")
                byBin <- "linear"
                ivcL <- structure_mapping(wm)
                barcode_plotter(ivcL, zVar, yVar, glue("{xx}-{yy}-{zz}-ivw-L"), bar_choice = "rgb", plots = TRUE, show = TRUE)
                byBin <- "quantile"
                ivcQ <- structure_mapping(wm)
                barcode_plotter(ivcQ, zVar, yVar, glue("{xx}-{yy}-{zz}-ivw-Q"), bar_choice = "rgb", plots = TRUE, show = TRUE)
                toc()

                print(glue("y={yVar}, z={zVar} completed"))
                gc(TRUE)
                print(format(Sys.time(), "%Y%m%d-%H:%M:%S"))
            }
            toc()
        }
    }
}
toc()
print(glue("Run complete at {Sys.time()}"))
