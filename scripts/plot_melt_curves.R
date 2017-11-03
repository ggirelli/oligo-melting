#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: plot single-oligo melting curves.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages({
	if (!require('argparser')) {
		install.packages("argparser", repos = "http://cran.uk.r-project.org")
		if (!require('argparser')) {
			stop(paste0('The R package argparser was not found. ',
				'Tried to automatically install the package but failed.',
				'Install it manually with install.packages(argparser)'))
		}
	}
	if (!require('ggplot2')) {
		install.packages("ggplot2", repos = "http://cran.uk.r-project.org")
		if (!require('ggplot2')) {
			stop(paste0('The R package ggplot2 was not found. ',
				'Tried to automatically install the package but failed.',
				'Install it manually with install.packages(ggplot2)'))
		}
	}
})

# INPUT ========================================================================

# Create arguent parser
parser = arg_parser('Plot single-oligo melting curves.',
	name = 'melt_curve_plot.R')

# Define mandatory arguments
parser = add_argument(parser, arg = 'input_tsv',
	help = paste0('Path to input tsv file. Three columns: oligo name, ',
	'temperature and dissociated fraction.'))
parser = add_argument(parser, arg = 'output_pdf',
	help = 'Path to output pdf file.')

# Define elective arguments
parser = add_argument(parser, arg = '--probe-name', short = '-n',
	type = class(''), help = 'Probe name.', default = "", nargs = 1)

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

# RUN ==========================================================================

# Read input
t = read.delim(input_tsv, as.is = T, header = F)

# Point to pdf
pdf(output_pdf, height = 10, width = 10)

# Structure plot
p = ggplot(t, aes(x = V2, y = V3, color = V1)) + geom_line() + guides(color = F)
p = p + xlab("Temperature [degC]") + ylab("Dissociated fraction")
p = p + ggtitle(probe_name)
print(p)

# Produce plot
graphics.off()

# END --------------------------------------------------------------------------

################################################################################
