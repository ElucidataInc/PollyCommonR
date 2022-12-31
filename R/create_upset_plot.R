#' create_upset_plot
#'
#' It creates an UpSet plot from list of vectors elements. For more information type help(UpSetR::upset)
#'
#' @param sets_list A list of sets to be compared.
#' @param nsets Number of sets to look at
#' @param nintersects Number of intersections to plot. If set to NA, all intersections will be plotted.
#' @param sets Specific sets to look at (Include as combinations. Ex: c("Name1", "Name2"))
#' @examples 
#' create_upset_plot(sets_list, nsets = 5, nintersects = 10, sets = NULL)
#' @import UpSetR
#' @export
create_upset_plot <- function(sets_list, nsets = 20, nintersects = 100, sets = NULL,
                              keep.order = F, set.metadata = NULL, intersections = NULL,
                              matrix.color = "gray23", main.bar.color = "gray23",
                              mainbar.y.label = "Intersection Size", mainbar.y.max = NULL,
                              sets.bar.color = "gray23", sets.x.label = "Set Size",
                              point.size = 3.0, line.size = 1.8, mb.ratio = c(0.7, 0.3),
                              expression = NULL, att.pos = NULL, att.color = main.bar.color,
                              order.by = c("freq", "degree"), decreasing = c(T, F),
                              show.numbers = "yes", number.angles = 0, group.by = "degree",
                              cutoff = NULL, queries = NULL, query.legend = "none",
                              shade.color = "gray88", shade.alpha = 0.25, matrix.dot.alpha = 0.5,
                              empty.intersections = NULL, color.pal = 1, boxplot.summary = NULL,
                              attribute.plots = NULL, scale.intersections = "identity",
                              scale.sets = "identity", text.scale = 1.1, set_size.angles = 0,
                              set_size.show = TRUE, set_size.numbers_size = NULL,
                              set_size.scale_max = NULL){
  message("Create UpSet Plot Started...")
  
  comb_elements <- unique(unlist(sets_list))
  upsetr_input <- data.frame(id = comb_elements)
  sets_names <- names(sets_list)
  sets_names <- sets_names[!is.na(sets_names)]
  
  for (set_name in sets_names){
    elements_bool <- as.numeric(comb_elements %in% sets_list[[set_name]])
    upsetr_input[[set_name]] <- elements_bool
  }
  
  if (identical(set_size.scale_max, NULL)){
    set_size.scale_max <- length(upsetr_input$id) + round(length(upsetr_input$id)*0.05, 0)
  }
  
  p <- UpSetR::upset(upsetr_input, nsets = nsets, nintersects = nintersects, sets = sets,
                     keep.order = keep.order, set.metadata = set.metadata, intersections = intersections,
                     matrix.color = matrix.color, main.bar.color = main.bar.color,
                     mainbar.y.label = mainbar.y.label, mainbar.y.max = mainbar.y.max,
                     sets.bar.color = sets.bar.color, sets.x.label = sets.x.label,
                     point.size = point.size, line.size = line.size, mb.ratio = mb.ratio,
                     expression = expression, att.pos = att.pos, att.color = att.color,
                     order.by = order.by, decreasing = decreasing,
                     show.numbers = show.numbers, number.angles = number.angles, group.by = group.by,
                     cutoff = cutoff, queries = queries, query.legend = query.legend,
                     shade.color = shade.color, shade.alpha = shade.alpha, matrix.dot.alpha = matrix.dot.alpha,
                     empty.intersections = empty.intersections, color.pal = color.pal, boxplot.summary = boxplot.summary,
                     attribute.plots = attribute.plots, scale.intersections = scale.intersections,
                     scale.sets = scale.sets, text.scale = text.scale, set_size.angles = set_size.angles,
                     set_size.show = set_size.show, set_size.numbers_size = set_size.numbers_size,
                     set_size.scale_max = set_size.scale_max
  )
  
  message("Create UpSet Plot Completed...")
  
  return (p)
}