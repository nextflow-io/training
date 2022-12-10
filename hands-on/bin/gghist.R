#!/usr/bin/env Rscript


options(stringsAsFactors=FALSE)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"), default="stdin",
	help="File or stdin [default=%default]"),

make_option(c("-o", "--output"), default="gghist.out.pdf",
	help="Output file name [default=%default]"),

make_option(c("-x", "--x_axis"), default=1,
	help="Index of the column with values, or labels if you already have counts [default=%default]"),

make_option(c("-y", "--y_axis"), default=NULL, type="integer",
	help="Index of the column with values, in case x provides counts. This will plot identity. Leave empty for default histogram [default=%default]"),

make_option(c("--header"), action="store_true", default=FALSE,
	help="Use this if the input has a header [default=%default]"),

make_option(c("--position"), default='dodge',
	help="Position for histogram [default=%default]"),

make_option(c("--scale_x_log10"), action="store_true", default=FALSE,
	help="log10-transform x scale [default=%default]"),

make_option(c("--scale_y_log10"), action="store_true", default=FALSE,
	help="log10-transform y scale [default=%default]"),

make_option(c("--y_title"), type="character", default="count",
	help="Title for the y axis [default=%default]"),

make_option(c("--x_title"), type="character", default=NULL,
	help="Title for the x axis [default=%default]"),

make_option(c("--title"), type="character", default=NULL,
	help="Title for the plot [default=%default]"),

make_option(c("-f", "--fill"), default="aquamarine",
	help="choose the color which you want to fill the histogram with"),

make_option(c("-c", "--color"), default="white",
	help="choose the color which you want to contour the histogram with"),

make_option(c("-F", "--fill_by"), type='numeric',
	help="the column index with the factor to fill by. Leave empty for no factor."),

make_option(c("-C", "--color_by"), type='numeric',
	help="the column index with the factor to color by. Leave empty for no factor."),

make_option(c("-A", "--alpha_by"), type='numeric',
	help="the column index with the factor to fill by. Leave empty for no factor."),

make_option(c("-P", "--palette"), 
    help='File with colors for the lines. Leave empty to use even color spacing'),

make_option(c("--sort"), action="store_true", default=FALSE,
	help="Sort the columns in decreasing order [default=%default]"),

make_option(c("--facet_by"), type='numeric',
	help="the column index with the factor to facet by. Leave empty for no factor."),

make_option(c("--facet_scale"), type='character', default="fixed",
	help="the scale of faceting: <fixed|free|free_y|free_x> [default=%default]"),

make_option(c("--facet_nrow"), type="numeric", 
	help="Number of row for faceting. Leave empty for auto [default=%default]"),

make_option(c("-W", "--width"), default=7,
	help="width of the plot in inches. [default=%default]"),

make_option(c("-H", "--height"), default=5,
	help="height of the plot in inches. [default=%default]"),

make_option(c("-B", "--base_size"), default=20,
	help="BAse size. [default=%default]"),

make_option(c("-b", "--binwidth"), type="double", 
	help="Specify binwidth. Leave empty for default"),

make_option(c("--flip"), action="store_true", default=FALSE,
	help="Flip coordinates [default=%default]"),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="if you want more output [default=%default]")

)

parser <- OptionParser(
	usage = "%prog [options] file", 
	option_list=option_list,
	description = "Reads the values on the first column and outputs a histogram"
)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}

#------------
# LIBRARIES
#------------ 

if (opt$verbose) {cat("Loading libraries... ")}
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(plyr))
if (opt$verbose) {cat("DONE\n\n")}


# ======== #
# BEGIN    #
# ======== #


# Read data
if (opt$input == "stdin") {input=file("stdin")} else {input=opt$input}
m = read.table(input, sep="\t", h=opt$header, quote=NULL) 

df = m


# Read facet
if (!is.null(opt$facet_by)) {facet_formula = as.formula(sprintf("~%s", colnames(df)[opt$facet_by]))}
# Read columns
x_col = colnames(df)[opt$x_axis]
if (!is.null(opt$y_axis)) {y_col = colnames(df)[opt$y_axis]}
if (!is.null(opt$fill_by)) {F_col = colnames(df)[opt$fill_by]}
if (!is.null(opt$color_by)) {C_col = colnames(df)[opt$color_by]}
if (!is.null(opt$alpha_by)) {A_col = colnames(df)[opt$alpha_by]}

if (!is.null(opt$fill_by)) {
	if (F_col == x_col) {
		df[paste(x_col, "fill", sep=".")] = df[,F_col]
		F_col = paste(x_col, "fill", sep=".")
	}
}

# Read palette
if (!is.null(opt$palette)) {
	palette = as.character(read.table(opt$palette, h=F, comment.char="%")$V1)
}

# Correct newlines if column is character
if (is.character(df[,x_col])) {
	df[,x_col] <- gsub("\\\\n", "\n", df[,x_col])
}

#================
# GGPLOT
#================

theme_set(theme_bw(base_size=opt$base_size))
theme_update(
	axis.text.x=element_text(angle=45, hjust=1, vjust=1),
	legend.key = element_rect(color='white'),
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank()
)




# Sort bars by abundance
if (opt$sort) {
	oldLev = levels(factor(df[,x_col]))
	if (is.null(opt$y_axis)) {
		lev = levels(as.factor(df[,x_col]))[order(table(df[,x_col]), decreasing=TRUE)]
		df[x_col] <- factor(df[,x_col], levels=lev)
	}
	if (!is.null(opt$y_axis)) {
		lev = df[order(df[,y_col], decreasing=TRUE), x_col]   # have to remove duplicated x
		df[x_col] <- factor(df[,x_col], levels=lev)
	}
}

# Params
geom_params = list()
#geom_params$color = opt$color
#geom_params$closed = c("right", "left")
#geom_params$include.lowest=TRUE

# specify binwidth
geom_params$binwidth = opt$binwidth

stat = "bin"

if (is.factor(df[,x_col]) | is.character(df[,x_col])) {
	stat = "count"
}

# Stat parameters 
stat = ifelse(is.null(opt$y_axis), stat, "identity")


stat_params = list(
	right=TRUE, 
	include.lowest=TRUE
)

mapping = list()

mapping <- modifyList(mapping, aes_string(x=x_col))

if (!is.null(opt$y_axis)) {
	mapping <- modifyList(mapping, aes_string(y=y_col))
}

# specify fill column
if (!is.null(opt$fill_by)) {
	mapping <- modifyList(mapping, aes_string(fill=F_col, order=rev(F_col)))
} else {
	geom_params$fill = opt$fill
}


# specify color column
if (!is.null(opt$color_by)) {
	mapping <- modifyList(mapping, aes_string(color=F_col, order=rev(F_col)))
} else {
	geom_params$color = opt$color
}

# specify alpha column
if (!is.null(opt$alpha_by)) {
	mapping <- modifyList(mapping, aes_string(alpha=A_col))
} 


class(mapping) <- "uneval"

# define histogram layer 
histLayer <- layer(
    geom = "bar",
    params = geom_params,
	position = opt$position,
	mapping = mapping,
    stat = stat
)


# start the plot
gp = ggplot(df) + histLayer

if (!is.character(df[,x_col]) & !is.factor(df[,x_col])) {
	avg = mean(df[,x_col], na.rm=TRUE)
	med = median(df[,x_col], na.rm=TRUE)
	gp = gp + geom_point(aes(x=avg, y=0), size=2)
	gp = gp + geom_vline(xintercept=med, linetype=2)
}

# Fill scale
if (!is.null(opt$fill_by)) {
	if (!is.null(opt$palette)) {
		gp = gp + scale_fill_manual(values=palette)
	} else {
		gp = gp + scale_fill_hue()
	}
}

# Color scale
if (!is.null(opt$color_by)) {
	if (!is.null(opt$palette)) {
		gp = gp + scale_color_manual(values=palette)
	} else {
		gp = gp + scale_color_hue()
	}
}

if (!is.null(opt$alpha_by)) {
	gp = gp + scale_alpha_discrete(range=rev(c(0.4,1)))
}

if (!is.null(opt$facet_by)) {
	gp = gp + facet_wrap(facet_formula, scales=opt$facet_scale, nrow=opt$facet_nrow)
}

if (opt$scale_x_log10) {gp = gp + scale_x_log10()}
if (opt$scale_y_log10) {gp = gp + scale_y_log10()}

if (!is.null(opt$x_title)) {gp = gp + labs(x=opt$x_title)}
if (!is.null(opt$title)) {gp = gp + ggtitle(opt$title)}

gp = gp + labs(y=opt$y_title)

gp = gp + coord_cartesian()

if (opt$flip) {
	gp = gp + coord_flip()
}

#gp = gp + geom_density(aes_string(x=x_col))

#gp = gp + scale_y_continuous(limits=c(0,20000))

ggsave(opt$output, h=opt$height, w=opt$width, title=opt$output)

# EXIT
quit(save='no')
