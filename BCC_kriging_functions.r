#' BCC Kriging Functions
#'
#' This script contains a collection of functions used for SCV kriging analysis and data output plotting. The functions are called from 'BCC_shared_kriging_run.r'.
#'
#' @section Functions:
#' - `barcode_plotter(df, title, ytitle, notes, ex_dir, bar_choice, del, plots, show)`: Plots IVW barcode sequences using ggplot2, with options for color mapping and exporting images.
#' - `ivw(krigeDF, geneLength, ex_dir, write, plots)`: Performs inverse variance weighting on kriging results, generates plots, and exports weighted means.
#' - `kriging_landscapes(krigemap, krigedata, p, ex_dir, plots)`: Creates kriging landscape visualizations and exports data.
#' - `kriging_process(kdt, notes, p, r_dir, ex_dir, xv, yv, zv, write, plots)`: Performs kriging with empirical variogram fitting, cross-validation, and result export.
#' - `do_krige(kdt, v_mf, p, ex_dir)`: Applies SCV kriging predciction to a regular grid using fitted variogram parameters and exports summary plots.
#' - `structure_mapping(wm, ex_dir)`: Maps weighted means to structure bins and generates PyMOL coloring scripts for visualization.

#' @section Dependencies:
#' - Requires packages: `ggplot2`, `data.table`, `glue`, `viridis`, and spatial analysis libraries.
#'
#' @section Usage:
#' These functions are designed to be modular and are called from the main analysis script. They support flexible plotting, data export, and integration with structural visualization tools.
#'
#' @author Ben C. Calverley
#' @date 2025

#--------------------------------------------------------
# %% Plot IVW barcode sequences
barcode_plotter <- function(df, title, ytitle, notes, ex_dir = export_dir, bar_choice = "z", del = FALSE, plots = FALSE, show = TRUE) {
    ## Plots barcode sequences using ggplot2, with options for color mapping and exporting images.
    indir <- getwd()
    if (bar_choice == "z") {
        bar <- ggplot(df, aes(x = x, y = 1, color = z))
    } else {
        rgb_cols <- rgb(df$R, df$G, df$B, maxColorValue = 255)
        bar <- qplot(1:length(rgb_cols), 1, fill = factor(1:length(rgb_cols)), geom = "tile") +
            scale_fill_manual(values = rgb_cols) +
            theme(legend.position = "none")
    }
    bar <- bar +
        geom_tile() +
        labs(x = "Sequence position", color = title, title = ytitle) +
        theme(rect = element_rect(fill = "transparent")) +
        theme(
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank()
        )
    if (del) {
        bar <- bar + scale_colour_distiller(type = "div", palette = "RdBu")
    } else {
        bar <- bar + scale_color_viridis(discrete = FALSE, option = "plasma")
    }
    if (plots) {
        setwd(file.path(ex_dir, glue('mdf-{str_replace_all(toString(maxd.fraction),fixed("."),"_")}')))
        ggsave(glue("{notes}.png"), plot = bar, width = 12, height = 3, dpi = 300, bg = "transparent")
        setwd(indir)
    }
    if (show) {
        return(bar)
    }
}

#--------------------------------------------------------
# %% Inverse variance weighting
ivw <- function(krigeDF, geneLength, ex_dir = export_dir, write = TRUE, plots = TRUE) {
    ## Performs inverse variance weighting on a kriging data frame. It calculates the weighted mean of the predicted values, using the inverse of the variance as the weight. It also calculates the standard deviation of the mean and generates a plot of the weighted z values. Finally, it saves the results to a file. The function takes several parameters, such as the input data frame (krigeDF), the gene length, and a flag to specify whether to write the results to a file and plot the data.
    tic("Inverse variance weighting")
    indir <- getwd()
    krigeDT <- data.table(krigeDF)
    meanVar <- krigeDT[, mean(var1.var), by = x]
    meanSd <- sqrt(meanVar)
    wm <- krigeDT[, weighted.mean(var1.pred, 1 / var1.var, na.rm = TRUE), by = x]
    wm$x <- 1 + (geneLength - 1) * wm$x
    zeroVars <- (wm$x[is.na(wm$V1)] - 1) / (geneLength - 1)
    for (xx in zeroVars) {
        wm$V1[wm$x == (1 + xx * (geneLength - 1))] <- mean(krigeDT$var1.pred[krigeDT$var1.var == 0 & mapply(function(x) {
            isTRUE(all.equal(x, unname(xx)))
        }, krigeDT$x)])
    }
    krigeDT$var1.pred[mapply(function(x) {
        isTRUE(all.equal(x, unname(xx)))
    }, krigeDT$x)]
    wm$low <- wm$V1 - meanSd$V1
    wm$low[wm$low < 0] <- 0
    wm$high <- wm$V1 + meanSd$V1
    wm$z <- wm$V1

    dir1 <- file.path(ex_dir, glue('mdf-{str_replace_all(toString(maxd.fraction),fixed("."),"_")}'))
    if (!file.exists(dir1)) {
        dir.create(ex_dir)
        dir.create(dir1)
    }
    setwd(dir1)
    if (plots) {
        weighted.z <- ggplot(wm, aes(x, V1)) +
            geom_point(size = .1) +
            ylim(0.95 * min(wm$low), 1.05 * max(wm$high)) +
            geom_ribbon(data = wm, aes(ymin = low, ymax = high), alpha = .3)
        ggsave("inv_var_weighted_z.png", plot = weighted.z)
    }
    if (write) {
        write.table(data.frame(wm$x, wm$V1, meanVar$V1), file = "weighted_z.csv", sep = ",", append = FALSE, quote = FALSE, col.names = TRUE, row.names = FALSE)
    }
    setwd(indir)
    toc()
    print(format(Sys.time(), "%Y%m%d-%H:%M:%S"))
    return(wm)
}

#--------------------------------------------------------
# %% Create kriging landscapes
kriging_landscapes <- function(krigemap, krigedata, p = pearson, ex_dir = export_dir, plots = TRUE) {
    ## Creates kriging landscape visualizations and exports data. It takes a kriged map and the original kriging data as input, along with parameters for exporting and plotting. The function generates a ggplot2 visualization of the kriged landscape, overlays the original data points, and saves the plot and data to files if specified. It returns the kriging data frame for further analysis.
    indir <- getwd()
    krigeDF <- as(krigemap, "data.frame")
    krigepdf <- as(krigedata, "data.frame")
    write.table(krigeDF, glue("{filename}.csv"), sep = ",", row.names = FALSE)
    if (plots) {
        setwd(ex_dir)
        krigeplot <- try(ggplot() +
            geom_tile(data = krigeDF, aes(x, y, color = var1.pred)) +
            theme(rect = element_rect(fill = "transparent")) +
            geom_point(data = krigepdf, aes(x, y), size = 1, color = "white") +
            coord_equal() +
            scale_color_viridis(discrete = FALSE, option = "plasma") +
            labs(color = zVar, fill = "Variance") +
            scale_fill_gradient(low = "black", high = "white") +
            theme(aspect.ratio = 1) +
            theme_bw() +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank()))
        if (!(class(krigeplot) %in% "try-error")) {
            varPlot <- ggplot(krigeDF, aes(x, y)) +
                geom_tile(aes(fill = var1.var)) +
                geom_point(data = krigepdf[, 1:2], size = .1) +
                coord_equal() +
                scale_fill_distiller(palette = "Reds", direction = 1) +
                ggtitle(glue("{xVar} OK variance r = {signif(p,3)}")) +
                xlab("Sequence position") +
                ylab(yVar) +
                labs(fill = "Variance") +
                theme(aspect.ratio = 1)
            ggsave("OKplot.png", plot = krigeplot)
            ggsave("Varplot.png", plot = varPlot)
            setwd(indir)
        } else {
            print("ERROR HERE LOOK AT ME!")
        }
    }
    return(krigeDF)
}

#--------------------------------------------------------
# %% Perform kriging
kriging_process <- function(kdt, notes = note, p = p_thresh, r_dir = parent_dir, ex_dir = export_dir, xv = xVar, yv = yVar, zv = zVar, write = TRUE, plots = TRUE) {
    ## Performs kriging on a given dataset using a range of distance cutoffs, and returns the Pearson correlation coefficient between the predicted and observed values, the range fraction that yielded the highest correlation, and the fitted variogram model.
    indir <- getwd()
    # Calculate empirical variogram and fit variogram model
    tic("Variogram and fit for kriging")
    pearson <- 0
    for (maxi in seq(.05, 1, .05)) {
        maxd1 <- maxi
        v1 <- variogram(z ~ 1, data = kdt, cutoff = maxd1 * maxd, width = maxd1 * maxd / 50)
        if (class(v1) == "NULL") {
            next
        }
        vmf1 <- fit.variogram(v1, vgm(c("Nug", "Wav", "Exp", "Sph", "Gau", "Exc", "Mat", "Ste", "Cir", "Lin", "Pen", "Per", "Hol", "Log", "Bes", "Spl")), debug.level = 1) # ,'Pow'
        vLine1 <- try(variogramLine(vmf1, maxdist = maxd1 * maxd, n = 1000))
        if (class(vLine1) %in% "try-error") {
            next
        }

        #  Krige cross validation and test accuracy of predictions
        nMin <- 5
        nMax <- 50
        kcv1 <- krige.cv(z ~ 1, locations = kdt, model = vmf1, nmin = nMin, nmax = nMax)
        pearson1 <- try(cor.test(kcv1$var1.pred, kcv1$observed, method = c("pearson")))
        if (class(pearson1) %in% "try-error") {
            next
        }

        if (is.na(pearson1$estimate)) {
            next
        }
        if ((pearson1$estimate)**2 > pearson**2) {
            maxd.fraction <- maxd1
            vmf <- vmf1
            kcv <- kcv1
            pearson <- pearson1$estimate
            vLine <- vLine1
            v <- v1
        }
        lin_pearson <- cor(kdt$y, kdt$z)
    }
    if (write) {
        setwd(r_dir)
        if (file.exists(glue("{dateStamp}-pearson_R_values_v2.csv"))) {
            r_csv <- read.csv(glue("{dateStamp}-pearson_R_values_v2.csv"), header = FALSE, numerals = "no.loss")
            if (length(r_csv[r_csv[2] == notes & r_csv[4] == yv & r_csv[5] == zv & mapply(function(x) {
                isTRUE(all.equal(x, unname(pearson)))
            }, r_csv[, 6]) & r_csv[7] == maxd.fraction & r_csv[8] == lin_pearson]) == 0) {
                write.table(data.frame(format(Sys.time(), "%Y%m%d%H%M%S"), notes, xv, yv, zv, pearson, maxd.fraction, lin_pearson), file = glue("{dateStamp}-pearson_R_values_v2.csv"), sep = ",", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)
            }
        } else {
            write.table(data.frame(format(Sys.time(), "%Y%m%d%H%M%S"), notes, xv, yv, zv, pearson, maxd.fraction, lin_pearson), file = glue("{dateStamp}-pearson_R_values_v2.csv"), sep = ",", append = TRUE, quote = FALSE, col.names = c("Date & time", "Note", "Protein (x)", "y", "z", "R", "range fraction", "linear R"), row.names = FALSE)
        }
    }
    if (pearson >= p) {
        if (plots) {
            if (!file.exists(ex_dir)) {
                dir.create(ex_dir)
            }
            dir.create(file.path(ex_dir, glue('mdf-{str_replace_all(toString(maxd.fraction),fixed("."),"_")}')))
            setwd(file.path(ex_dir, glue('mdf-{str_replace_all(toString(maxd.fraction),fixed("."),"_")}')))
            print(getwd())
            ggsave("variogram_scatter.png", plot = ggplot(v, aes(dist, gamma)) +
                geom_point(size = .5))
            resPlot <- ggplot(as(kcv, "data.frame"), aes(var1.pred, residual)) +
                geom_point(size = .1) +
                geom_hline(yintercept = 0, linetype = 2) +
                ggtitle(paste(xVar, zVar, "OK prediction residuals", sep = " "))
            ggsave("predVresidual.png", plot = resPlot)
            vPlot <- ggplot(v, aes(dist, gamma)) +
                geom_point(size = .5) +
                geom_line(data = vLine, colour = "blue") +
                theme_bw(base_size = 36) +
                theme(axis.title.x = element_blank(), axis.title.y = element_blank())
            ggsave("variogram.png", plot = vPlot)
        }
    }
    setwd(indir)
    toc()
    print(format(Sys.time(), "%Y%m%d-%H:%M:%S"))
    return(list(out1 = pearson, out2 = maxd.fraction, out3 = vmf, out4 = lin_pearson))
}

#--------------------------------------------------------
# %% Perform kriging
do_krige <- function(kdt, v_mf, p = p_thresh, ex_dir = export_dir) {
    ## Uses the parameters obtained from kriging_process() to perform kriging on a regular grid of points, and returns a kriged map.
    indir <- getwd()
    tic("Do krige")
    if (pearson >= p) {
        krigegrid <- expand.grid(x = seq(0, 1, by = 1 / (geneLength - 1)), y = seq(0, 1, by = 0.01))
        coordinates(krigegrid) <- ~ x + y
        gridded(krigegrid) <- TRUE
        krigemap <- krige(z ~ 1, locations = kdt, newdata = krigegrid, model = v_mf, nmin = nMin, nmax = nMax)
        summary_export(krigemap, kdt, export_dir)
    } else {
        krigemap <- c()
        print("NULL krigemap")
    }
    setwd(indir)
    toc()
    print(format(Sys.time(), "%Y%m%d-%H:%M:%S"))
    return(krigemap)
}

#--------------------------------------------------------
# %% Map to structure
structure_mapping <- function(wm, ex_dir = export_dir) {
    ## Maps weighted means to structure bins and generates PyMOL coloring scripts for visualization.
    indir <- getwd()
    if (byBin == "quantile") {
        n <- 20 # length(wm$V1)/20
        colours <- try(data.frame(table(cut_number(wm$V1, n)))[1])
        if (class(colours) %in% "try-error") {
            next
        }
        names(colours)[1] <- "bin"
        colours[c("R", "G", "B")] <- transpose(data.frame(col2rgb(plasma(n))))
        wm$bin <- cut_number(wm$V1, n)
    } else if (byBin == "linear") {
        n <- 100
        colours <- try(data.frame(table(cut_interval(wm$V1, n)))[1])
        if (class(colours) %in% "try-error") {
            next
        }
        names(colours)[1] <- "bin"
        colours[c("R", "G", "B")] <- transpose(data.frame(col2rgb(plasma(n))))
        wm$bin <- cut_interval(wm$V1, n)
    }
    wm$bin <- as.character(wm$bin)
    invVarcolour <- merge(wm[, c("x", "V1", "bin")], colours, by.x = "bin", by.y = "bin", sort = F)
    invVarcolour <- invVarcolour[invVarcolour$x != 0, ]

    setwd(file.path(ex_dir, glue('mdf-{str_replace_all(toString(maxd.fraction),fixed("."),"_")}')))
    unlink(glue("pymol-{byBin}.txt"))
    for (i in 1:nrow(invVarcolour)) {
        cat("select ", i, ", resi ", as.character(invVarcolour[i, 2]), " and chain A", "\n",
            "set_color color", i, "=[", as.character(invVarcolour[i, 4]), ",", as.character(invVarcolour[i, 5]), ",", as.character(invVarcolour[i, 6]), "]", "\n",
            "color color", i, ",", i, "\n",
            file = glue("pymol-{byBin}.txt"), sep = "", append = TRUE
        )
        cat("select ", i, ", resi ", as.character(invVarcolour[i, 2]), " and chain B", "\n",
            "set_color color", i, "=[", as.character(invVarcolour[i, 4]), ",", as.character(invVarcolour[i, 5]), ",", as.character(invVarcolour[i, 6]), "]", "\n",
            "color color", i, ",", i, "\n",
            file = glue("pymol-{byBin}.txt"), sep = "", append = TRUE
        )
    }
    setwd(indir)
    return(invVarcolour)
}
