#' create_pairwise_posthoc_boxplot
#'
#' Creates a pairwise post-hoc boxplot for selected metabolite and comparison.
#'
#' @param posthoc_data_long_df A dataframe containing long format post-hoc data.
#' @param intensity_data A dataframe containing intensity data.
#' @param selected_metabolite The selected metabolite for analysis.
#' @param selected_comparison The selected comparison for analysis.
#' @return A plotly object representing the pairwise post-hoc boxplot.
#' @import dplyr
#' @import plotly
#' @export
create_pairwise_posthoc_boxplot <- function(posthoc_data_long_df, intensity_data, selected_metabolite, selected_comparison) {
  if (!is.null(selected_metabolite) && !is.null(selected_comparison)) {
    message("Plot Posthoc Started...")
    require(dplyr)
    require(plotly)
    
    # Filter the posthoc_result dataframe based on the selected metabolite and interaction
    filtered_posthoc <- dplyr::filter(posthoc_data_long_df, Comparison == selected_comparison & ID == selected_metabolite)
    
    filtered_posthoc <- filtered_posthoc %>%
      dplyr::mutate(group_1 = sub(":.*", "", Comparison),
             group_2 = sub(".*:", "", Comparison))
    
    selected_interaction <- unique(filtered_posthoc$Interaction)
    
    # Filter the intensity_data dataframe based on the selected metabolite
    filtered_intensity_data <- dplyr::filter(intensity_data, id == selected_metabolite)
    
    # Subset filtered_intensity_data based on group1 and group2
    filtered_intensity_data <- filtered_intensity_data %>%
      dplyr::filter(get(selected_interaction) %in% c(filtered_posthoc$group_1, filtered_posthoc$group_2))
    
    p <- plotly::plot_ly(
      data = filtered_intensity_data,
      x = ~.data[[selected_interaction]],
      y = ~value,
      type = "box",
      boxpoints = "all",
      jitter = 0.3,
      pointpos = 0,
      marker = list(color = ~.data[[selected_interaction]]),
      boxmean = TRUE,
      fillcolor = ~.data[[selected_interaction]],
      text = ~paste("Sample:", .data[["variable"]], "<br>Intensity:", .data[["value"]]), # Tooltip text
      hoverinfo = "text" # Show tooltip text
    ) %>%
      layout(
        xaxis = list(
          title = paste(selected_metabolite, "\n", selected_interaction),
          titlefont = list(size = 16),  # Adjust x-axis title font size
          showticklabels = FALSE  # Hide tick labels
        ),
        yaxis = list(
          title = "Normalized Intensity",
          titlefont = list(size = 16)  # Adjust y-axis title font size
        ),
        paper_bgcolor = "white",
        plot_bgcolor = "white"
      )

    # Format adjusted p-value and add text annotation
    adj_p_formatted <- sprintf("%.6f", filtered_posthoc$adj.p)

    # Add text annotation for adjusted p-value
    p <- p %>%
      add_annotations(
        x = 0.5,
        y = max(filtered_intensity_data$value),
        text = paste("Adj. p =", adj_p_formatted),
        showarrow = FALSE,
        font = list(size = 12)
      )

    p <- update_plotly_config(p)
    
    message("Plot Posthoc Completed...")
    return(p)
  }
}