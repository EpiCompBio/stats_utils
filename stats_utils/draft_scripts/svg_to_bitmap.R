######################
# Import SVGs into R and process into multi-plot figure
# Antonio Berlanga-Taylor
# July 2018
######################


######################
# source('http://bioconductor.org/biocLite.R')
# biocLite('grImport2')

# For image rendering from SVGs:
library(rsvg)
# https://cran.r-project.org/web/packages/rsvg/index.html
# https://www.opencpu.org/posts/svg-release/
# https://github.com/jeroen/rsvg

# Instead of rsvg maybe use grImport2 as bitmap array needs convert to R graphics
# object in order to create multi-plots and work interactively:
# https://www.stat.auckland.ac.nz/~paul/R/grImport2/grImport2.pdf
# https://stackoverflow.com/questions/50325139/plot-graph-in-device
library(grImport2)

# image creation but not needed for rendering:
library(svglite) # prefer over base R svg() but not always! Faster but not better quality
# https://blog.rstudio.com/2015/12/10/svglite-1-0-0/
library(ggplot2)

# For processing multi-plot figures, labels, etc.:
library(cowplot)
# https://cran.r-project.org/web/packages/cowplot/vignettes/plot_grid.html

# Extensive editing in R might be done with:
library(magick)
# https://cran.r-project.org/web/packages/magick/vignettes/intro.html
# https://stackoverflow.com/questions/9917049/inserting-an-image-to-ggplot2# See:
######################

######################
# setwd('imputation/')

# # Get file names:
# input_name_1 <- 'missingness_pattern_nhanes.svg'
# input_name_2 <- 'missingness_pattern_nhanes.svg'

# Read in files:
# plot_1 <- rsvg(input_name)
# plot_2 <- rsvg(input_name)

# Adjust with cowplot
# fig_x <- plot_grid(plot_1, plot_2, labels = c('A', 'B'))

# Render and save with rsvg:
# rsvg_pdf(fig_x, "fig_x.pdf")

# Next steps might be:
# Pull in rendered (eg PDFs) back in to an Rmd file
# Put together with legends
# Put together with other figures into one document
######################

######################
# Examples
# Create an svg image and save it:
svglite("plot.svg") #, width = 10, height = 7)
qplot(mpg, wt, data = mtcars, colour = factor(cyl))
dev.off()

# Create a second image:
svglite("plot2.svg") #, width = 10, height = 7)
qplot(mpg, wt, data = mtcars, colour = factor(cyl))
dev.off()

# Read in plot and render it into a bitmap array:
bitmap_1 <- rsvg("plot.svg")
dim(bitmap)
class(bitmap)
str(bitmap)
bitmap_2 <- rsvg("plot2.svg")

# Check if magick is better:
plot_1_mag <- image_read_svg('plot.svg') # requires the rsvg package, gives better quality
class(plot_1_mag)
print(plot_1_mag)
plot_2_mag <- image_read_svg('plot2.svg') # requires the rsvg package, gives better quality

# Check if grImport2 is better
# And at the same time convert bitmap array to graphics object for further
# processing
# plot_1_gr <- readPicture("plot.svg") # errors
# Try:
# https://stackoverflow.com/questions/50325139/plot-graph-in-device

library(magrittr)
infile <- 'plot.svg'
plot_1_gr <- system(paste("cat", infile), intern = TRUE) %>%
	paste0(., collapse = "") %>%
	charToRaw(.) %>%
	rsvg::rsvg_svg(NULL, file = NULL) %>%
	rawToChar(.) %>%
	grImport2::readPicture(.)

grImport2::grid.picture(plot_1_gr)
class(plot_1_gr)

infile <- 'plot2.svg'
plot_2_gr <- system(paste("cat", infile), intern = TRUE) %>%
	paste0(., collapse = "") %>%
	charToRaw(.) %>%
	rsvg::rsvg_svg(NULL, file = NULL) %>%
	rawToChar(.) %>%
	grImport2::readPicture(.)

grImport2::grid.picture(plot_2_gr)

# Convert bitmap arrays into class "ggplot", "gtable", "grob", "recordedplot",
# or a function that plots to an R graphicsdevice when called
# https://cran.r-project.org/web/packages/magick/vignettes/intro.html#raster_images
plot_1 <- plot_1_gr
plot_2 <- plot_2_gr

# Adjust with cowplot:
plot_grid(bitmap_1, bitmap_2, labels = c('A', 'B')) # errors, wrong object
plot_grid(plot_1, plot_2, labels = c('A', 'B')) # errors, wrong object

# TO DO: fix above and below!
# These don't communicate, different devices:
plot_grid(NULL, NULL, labels = c('A', 'B'))
plot_grid(grImport2::grid.picture(plot_1_gr),
					grImport2::grid.picture(plot_2_gr),
					# labels = c('A', 'B')
					)
dev.off()



# Write to format with other tools:
# png::writePNG(bitmap, "bitmap.png")
# jpeg::writeJPEG(bitmap, "bitmap.jpg", quality = 1)
# webp::write_webp(bitmap, "bitmap.webp", quality = 100)

# Write to format with rsvg:
rsvg_svg(fig_x, "fig_x.svg")

# You can also read in plot and render it directly into output format with rsvg:
rsvg_pdf("plot.svg", "out.pdf")
rsvg_svg("plot.svg", "out.svg")

# Cowplot examples
# https://cran.r-project.org/web/packages/cowplot/vignettes/plot_grid.html

theme_set(theme_cowplot(font_size=12)) # reduce default font size
# Plot 1:
plot.mpg <- ggplot(mpg, aes(x = cty, y = hwy, colour = factor(cyl))) +
	geom_point(size=2.5)

# Plot 2:
plot.diamonds <- ggplot(diamonds, aes(clarity, fill = cut)) + geom_bar() +
	theme(axis.text.x = element_text(angle=70, vjust=0.5))

# Adjust with cowplot:
plot_grid(plot.mpg, plot.diamonds, labels = c('A', 'B'))

# Render and save with rsvg:
# rsvg_svg(fig_x, "fig_x.svg")
######################