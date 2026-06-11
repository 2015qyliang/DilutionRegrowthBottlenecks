require(tidyverse)
require(scales)
require(RColorBrewer)

mytheme <- theme_bw() + theme(
  legend.title  = element_text( size=17),
  #  legend.position = "bottom",
  #	legend.direction = "horizontal",
  legend.key = element_blank(),
  legend.text  = element_text( size=17),
  panel.background = element_rect(fill=NA),
  panel.grid = element_blank(),
  text = element_text( family="Helvetica", size=19),
  panel.border = element_rect( colour = "black", size=2),
  axis.ticks = element_line(size = 1.),
  legend.background = element_rect(fill=NA)
)

mytheme_main <- theme_bw() + theme(
  legend.title  = element_text(family="Helvetica", size=17, color = "#222222"),
  legend.key = element_blank(),
  legend.text  = element_text(family="Helvetica", size=17, color = "#222222"),
  panel.background = element_rect(fill="transparent"),
  #plot.background = element_rect(fill="transparent", colour = NA),
  panel.grid = element_blank(),
  text = element_text( family="Helvetica", size=17, color = "#222222"),
  panel.border = element_blank(),
  axis.title = element_text( family="Helvetica", size=10, color = "#222222"),
  axis.text = element_text( family="Helvetica", size=8, color = "#222222"),
  axis.line = element_line(size = 1., color = "#222222"),
  axis.ticks = element_line(size = 1.,color = "#222222"),
  legend.background = element_rect(fill="transparent", colour = NA)
)


fancy_linear <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = FALSE)
  # return this as an expression
  parse(text=l)
}



fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # e+00 becomes 1
  l <- gsub("e\\+00", "", l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove prefactor 1
  l <- gsub("'1'e", "10^", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # remove plus
  l <- gsub("\\+", "", l)
  # return this as an expression
  parse(text=l)
}


fancy_scientificb <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove prefactor 1
  l <- gsub("'1'e", "10^", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # remove plus
  l <- gsub("\\+", "", l)
  # return this as an expression
  parse(text=l)
}


scientific_10_exp_labels <- trans_format("log10", math_format(10^.x) )
scientific_10_exp_breaks <- trans_format("log10", function(x) 10^x )

# scalecols <- scale_colour_manual(values = c(
#   "D0" = "#00d667", "D1" = "#0000FF", 
#   "D3" = "#9933FF", "D5" = "#6699FF",
#   "D6" = "#FF0000", "D7" = "#FF6600" 
#   ) )

scalecols <- scale_colour_manual(values = pal_nejm()(6),
                                 breaks = c('D0', 'D1', 'D3', 'D5', 'D6', 'D7') )

# scaleshapes = scale_shape_manual( values = c( "D0" = 2, "D1" = 3, "D3" = 5, "D5" = 9, "D6" = 0, "D7" = 7))
scaleshapes = scale_shape_manual(values = rep(16, 6),
                                 breaks = c('D0', 'D1', 'D3', 'D5', 'D6', 'D7') )
