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

suppressMessages(library(argparser))
library(ggplot2)

# INPUT ========================================================================

# Create arguent parser
parser = arg_parser(
	'Plot single-oligo hybridization / secondary structure melting curves.',
	name = 'melt_curve_plot.R')

# Define mandatory arguments
parser = add_argument(parser, arg = 'target_tsv',
	help = paste0('Path to input target tsv file. Three columns: oligo name, ',
	'temperature and dissociated fraction.'))
parser = add_argument(parser, arg = 'second_tsv',
	help = paste0('Path to input secondary structure tsv file. ',
		'Three columns: oligo name, temperature and dissociated fraction.'))
parser = add_argument(parser, arg = 'output_pdf',
	help = 'Path to output pdf file.')

# Define elective arguments
parser = add_argument(parser, arg = '--probe-name', short = '-n',
	type = class(''), help = 'Probe name.', default = "", nargs = 1)
parser = add_argument(parser, arg = '--thybr', short = '-t', type = class(0.),
	help = 'Hybridization temperature.', default = NA, nargs = 1)
parser = add_argument(parser, arg = "--trange", type = class(0),
	help = 'Temperature range for plot limits.', default = 40, nargs = 1)

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

# Manipulations
xrange = c(thybr - trange / 2, thybr + trange / 2)

# RUN ==========================================================================

# Read input
t = read.delim(target_tsv, as.is = T, header = F)
s = read.delim(second_tsv, as.is = T, header = F)

# Point to pdf
pdf(output_pdf, height = 10, width = 10)

for (name in unique(t$V1)) {

	subt = t[t$V1 == name,]
	subt$type = "target"
	subs = s[do.call(rbind, strsplit(s$V1, ','))[,1] == name,]
	subs$type = "secondary structure"

	merged = rbind(subt, subs)

	# Structure plot
	p = ggplot(merged, aes(x = V2, y = V3, color = type))
	p = p + geom_line()# + guides(color = F)
	p = p + xlab("Temperature [degC]") + ylab("Dissociated/Unfolded fraction")
	if ( !is.na(thybr) ) p = p + geom_vline(
		aes(xintercept = thybr, color = "Thybrid"), linetype = 2)
	p = p + geom_hline(aes(yintercept = 0.5,
		color = "Melting point"), linetype = 2)
	p = p + ggtitle(name) + xlim(xrange)
	options(warn=-1); print(p); options(warn=0)

}

# Produce plot
graphics.off()

# END --------------------------------------------------------------------------

################################################################################
