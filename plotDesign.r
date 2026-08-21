######## 
# COLORS
################################################################


### Colors 
colors <- list()

# base
colors$base$mild <- "lightsteelblue4"
colors$base$moderate <- "#F39F07"
colors$base$severe <- "#EB5D12"
colors$base$PBC <- "#7BC086"
colors$base$AIH <- "#EB5D12"
# mHAI Plot
colors$mHAI_tri <- c(mild = colors$base$mild, moderate = colors$base$moderate, severe = colors$base$severe)
colors$mHAI_tri_with_PBC <- c(PBC = colors$base$PBC, colors$mHAI_tri)
# ALT/AST
colors$ast_alt <- c(astULN = colors$base$moderate, altULN = colors$base$mild)


# palette
## Heatmap Scaling
colors$palette$DarkBlue_Red$basefunc <- colorRampPalette(c("#3F517C", "whitesmoke", "#EB5D12"))
colors$palette$DarkBlue_Red$ggcolor <- function(...)  scale_color_gradient2(low = "#3F517C", 
                                                                                     high = "#EB5D12", 
                                                                                     mid = "whitesmoke", 
                                                                                     ...)
colors$palette$DarkBlue_Red$ggfill <- function(...)  scale_fill_gradient2(low = "#3F517C", 
                                                                                     high = "#EB5D12", 
                                                                                     mid = "whitesmoke", 
																				      ...)

## PBC vs AIH Scaling
colors$palette$Green_Red$basefunc <- colorRampPalette(c("#7BC086", "whitesmoke", "#EB5D12"))
colors$palette$Green_Red$ggcolor <- function(...)  scale_color_gradient2(low = "#7BC086", 
                                                                                     high = "#EB5D12", 
                                                                                     mid = "whitesmoke", 
                                                                                     midpoint = midpoint)
colors$palette$Green_Red$ggfill <- function(...)  scale_fill_gradient2(low = "#7BC086", 
                                                                                     high = "#EB5D12", 
                                                                                     mid = "whitesmoke", 
                                                                                     midpoint = midpoint)


######## 
# ggTheme
################################################################
theme_paper <- 
function(bordertype = "closed", # closed (surrounded by border), open (x an y axis), blank (no axis line / border)
         aspect.ratio = 1/1, # ratio of plot -> y/x
         base_size = 12, # font size
         base_family = "Helvetica", # font familiy
         legend_position = "right", 
         title_size= base_size){
    
    theme_closed <- 
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        
    theme(
        aspect.ratio = aspect.ratio, 
        # Panel
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "#000000",
                                    linewidth= unit(units = "npc", x = 0.5)), 
        
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        # Axis
        axis.line = element_blank(), axis.ticks = element_line(colour = "#000000"),
        axis.text = element_text(size = base_size), plot.title=element_text(size = title_size),
        axis.title = element_text(size = base_size), legend.text=element_text(size = base_size),
        # Legend
        legend.key = element_blank(),
        legend.title=element_text(size = base_size), strip.text=element_text(size = base_size),
        legend.position = legend_position, 
        legend.background = element_blank()
    )
    
    if(bordertype == "open")
        return(theme_closed %+replace% theme(axis.line = element_line(colour = "#000000", linewidth= unit(units = "npc", x = 0.5)), panel.border = element_blank()))
    else if(bordertype == "blank")
        return(theme_closed %+replace% theme(axis.line = element_blank(), panel.border = element_blank()))
    
    return(theme_closed)
    
}




# Example Application

# library(ggplot2)
# data <- data.frame(data1 = rnorm(n = 20, mean = 0, sd = 10), 
#                   data2 = rnorm(n = 20, mean = 0, sd = 10), 
#                   label = sample(LETTERS, 20))
#
#ggplot(data,
#      aes(data1, data2, color = data1)) + geom_point() + colors$palette$Green_Red$ggcolor() +
#theme_paper(bordertype =  "open")
#
#
#ggplot(data,
#      aes(y=label, x= "", fill = data1)) + geom_tile() + 
#colors$palette$DarkBlue_Red$ggfill() +
#theme_paper(bordertype =  "blank")
#
#heatmap(x = as.matrix(data[,1:2]), col=colors$palette$DarkBlue_Red$basefunc(50))