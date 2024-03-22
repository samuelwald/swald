#' Plot cluster proportions with highlighted clusters

#' @export 

dist_plot_highlight <- function(
    cluster, 
    cols, 
    percent = 50, 
    data,
    split_by,
    fill_by
){
    t_col <- function(
    color, 
    percent
    ){

    ## Get RGB values for named color
    rgb.val <- col2rgb(color)

    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100 - percent) * 255 / 100)

    ## Save the color
    return(t.col)
}
    `%ni%` <- Negate(`%in%`)
    for (i in names(cols)){
        if (i %ni% cluster){
            cols[[i]] <- t_col(color = cols[[i]], percent = percent)
        }
    }
    plot <- ggplot(data=data, aes(split_by)) +
                geom_bar(aes(fill=fill_by), position="fill") +
                theme_cowplot() + 
                theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
                theme(legend.title = element_blank()) +
                theme(aspect.ratio = 1.5) +
                scale_fill_manual(values= cols)
    return(plot)
}