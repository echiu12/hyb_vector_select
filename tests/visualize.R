libs=c("ggplot2", "RColorBrewer", "grid", "cowplot")
lapply(libs, function(x) if (!require(x, character.only=TRUE)) install.packages(x, dependencies=TRUE))
lapply(libs, library, character.only=TRUE)

options(width=200)

# https://github.com/wilkelab/cowplot/issues/202
get_legend_35 <- function(plot) {
  # return all legend candidates
  legends <- get_plot_component(plot, "guide-box", return_all = TRUE)
  # find non-zero legends
  nonzero <- vapply(legends, \(x) !inherits(x, "zeroGrob"), TRUE)
  idx <- which(nonzero)
  # return first non-zero legend if exists, and otherwise first element (which will be a zeroGrob) 
  if (length(idx) > 0) {
    return(legends[[idx[1]]])
  } else {
    return(legends[[1]])
  }
}

# Parse arguments.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 & length(args) != 2) {
    stop("Requires one or two positional arguments.\n If one positional argument is provided, it is the directory in which all .out files will be visualized.\n If it is two positional arguments, it is interpreted as in_file out_file")
}

if (length(args) == 1) {
    in_list <- list.files(path=args[1], pattern="bench.*\\.out", full.names=TRUE)
    out_list <- gsub(".out", ".pdf", in_list)
    f_list<-data.frame(in_list,out_list)
} else {
    in_list <- c(args[1])
    out_list <- c(args[2])
    f_list<-data.frame(in_list,out_list)
}

print(in_list)

for(i in 1:nrow(f_list)) {
    print(f_list[i, ]$in_list)

    in_file <- f_list[i, ]$in_list
    out_file <- f_list[i, ]$out_list

    # Construct dataframe
    data <- read.table(in_file, header=TRUE)

    # Additional features for plotting
    data$bit_vector_group <- unlist(lapply(data$index_type, function(x) unlist(strsplit(x, "<"))[1]))

    # Groups
    data$bit_vector_group <- gsub("sparse_hyb_vector", "HYB", data$bit_vector_group)
    data$bit_vector_group <- gsub("sparse_sd_vector", "SD", data$bit_vector_group)

    data$bit_vector_group <- gsub("hyb_vector", "HYB", data$bit_vector_group)
    data$bit_vector_group <- gsub("rrr_vector", "RRR", data$bit_vector_group)
    data$bit_vector_group <- gsub("sd_vector", "SD", data$bit_vector_group)
    data$bit_vector_group <- gsub("bit_vector", "BV", data$bit_vector_group)
    data$bit_vector_group <- gsub("uncompressed", "BV", data$bit_vector_group)

    data$bit_vector_group <- gsub("huff", "HUFF", data$bit_vector_group)
    data$bit_vector_group <- gsub("blcd", "BLCD", data$bit_vector_group)

    data$bit_vector_group <- gsub("^HYB$", "PLCP_HYB", data$bit_vector_group)
    data$bit_vector_group <- gsub("^RRR$", "PLCP_RRR", data$bit_vector_group)
    data$bit_vector_group <- gsub("^SD$", "PLCP_SD", data$bit_vector_group)
    data$bit_vector_group <- gsub("^BV$", "PLCP_BV", data$bit_vector_group)

    data$bit_vector_type <- factor(
        unlist(lapply(data$bit_vector_group, function(x) tail(unlist(strsplit(x, "_")), n=1))),
        levels=c(
            "HYB",
            "RRR",
            "SD",
            "BV"
        ),
        ordered=TRUE
    )
    
    data$bit_vector_factor <- factor(
        data$bit_vector_group, 
        levels=c(
            "PLCP_HYB", 
            "HUFF_HYB", 
            "BLCD_HYB", 
            "PLCP_RRR", 
            "HUFF_RRR", 
            "BLCD_RRR", 
            "PLCP_SD", 
            "HUFF_SD", 
            "BLCD_SD", 
            "PLCP_BV",
            "HUFF_BV", 
            "BLCD_BV"
        ), 
        ordered=TRUE
    )
    print(data$bit_vector_group)
    print(data$bit_vector_factor)
    
    # Plot
    pdf(out_file, width=4, height=2)
    
    fnames <- unique(data$text_name)
    
    for (f in fnames) {
        data_subset <- data[data$text_name==f,]
    
        # data_subset$order <- unlist(lapply(data_subset$bit_vector_group, function(x) if(grepl("HYB", x, fixed = TRUE)) 1 else 0))
        # data_subset$order <- data_subset$order + unlist(lapply(data_subset$bit_vector_group, function(x) if(grepl("RRR", x, fixed = TRUE)) 2 else 0))
        # data_subset$order <- data_subset$order + unlist(lapply(data_subset$bit_vector_group, function(x) if(grepl("SD", x, fixed = TRUE)) 3 else 0))
        # data_subset$order <- data_subset$order + unlist(lapply(data_subset$bit_vector_group, function(x) if(grepl("BV", x, fixed = TRUE)) 4 else 0))
        # print(data_subset$order)
    
        # data_subset <- data_subset[order(data_subset$order),]
        # print(data_subset$order)
    
        p <- ggplot(
                data_subset,
                aes(
                    x=space_pct, 
                    y=time_us, 
                    color=bit_vector_factor,
                    shape=bit_vector_factor,
                    group=bit_vector_factor
                )
            )
        p <- p + theme_bw()
        p <- p + theme(
            text = element_text(family = "serif"),
            legend.title=element_blank(),
            legend.position = "none",
        )
        if (any(grepl('BLCD', data_subset$bit_vector_group))) {
        # if (sum(str_detect(data_subset$bit_vector_group, 'BLCD')) > 0) {
            # p <- p + scale_color_brewer(palette = "Paired")
            # p <- p + scale_shape_manual(values=c(15, 15, 16, 16, 17, 17, 18, 18))
            # p <- p + scale_size_manual(values = c(3, 3, 3, 3, 3, 3, 4, 4))
    
            colors <- brewer.pal(12, "Paired")
            p <- p + scale_color_manual(values = c(colors[1], colors[2], colors[3], colors[4], colors[5], colors[6], colors[7], colors[8]))
            p <- p + scale_shape_manual(values = c(15, 15, 16, 16, 17, 17, 18, 18))
            p <- p + scale_size_manual(values = c(3, 3, 3, 3, 3, 3, 4, 4))
        }
        else {
            colors <- brewer.pal(12, "Paired")
            p <- p + scale_color_manual(values = colors[seq_along(colors) %% 2 == 0])
            p <- p + scale_shape_manual(values=c(15, 16, 17, 18))
            p <- p + scale_size_manual(values = c(3, 3, 3, 4))
        }
        p <- p + geom_point(aes(size=bit_vector_factor))
        p <- p + geom_path()
        p <- p + labs(
            title=f, 
            x="Space (%)", 
            y=expression(paste("Time (", mu, "s)")), 
            color="", 
            shape="",
            size=""
        )
        p <- p + expand_limits(x = 0, y = 0)
        # p <- p + scale_x_continuous(expand = c(0, 0), limits = c(0, NA))
        # p <- p + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
        print(p)
    }
    
    p <- p + theme(
        text = element_text(family = "serif"),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.key.size = unit(0.0, "cm"),
    )
    grid.newpage()
    grid.draw(get_legend_35(p))
    
    dev.off()

}
