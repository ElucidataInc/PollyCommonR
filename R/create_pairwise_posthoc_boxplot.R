#' create_pairwise_posthoc_boxplot
#'
#' Creates a pairwise post-hoc boxplot for selected metabolite and interaction.
#'
#' @param posthoc_data_long_df A dataframe containing long format post-hoc data.
#' @param intensity_data A dataframe containing intensity data.
#' @param selected_metabolite The selected metabolite for analysis.
#' @param selected_interaction The selected interaction for analysis.
#' @return A plotly object representing the pairwise post-hoc boxplot.
#' @import dplyr
#' @import plotly
#' @import ggpubr
#' @export
create_pairwise_posthoc_boxplot <- function(posthoc_data_long_df = NULL, intensity_data = NULL, selected_metabolite = NULL, selected_interaction = NULL) {
  if (selected_metabolite != "" && selected_interaction != "" && !is.null(posthoc_data_long_df) && !is.null(intensity_data)) {
    message("Plot Posthoc Started...")
    require(dplyr)
    require(plotly)
    require(ggpubr)

    # Filter the posthoc_result dataframe based on the selected metabolite and interaction
    filtered_posthoc <- dplyr::filter(posthoc_data_long_df, Interaction == selected_interaction & ID == selected_metabolite)

    if (nrow(filtered_posthoc) == 0) {
      message("No posthoc data found for the selected metabolite and interaction.")
    }

    if (nrow(intensity_data) == 0) {
      message("No intensity data found for the selected metabolite and interaction.")
    }

    filtered_posthoc <- filtered_posthoc %>%
      dplyr::mutate(group1 = sub(":.*", "", Comparison),
                    group2 = sub(".*:", "", Comparison))

    selected_interaction <- unique(filtered_posthoc$Interaction)

    # Filter the intensity_data dataframe based on the selected metabolite
    filtered_intensity_data <- dplyr::filter(intensity_data, id == selected_metabolite)

    # Combine desired columns with selected_interaction
    column_to_keep <- c("id", "value", selected_interaction)
    filtered_intensity_data <- filtered_intensity_data[, column_to_keep]

    # Adjust interaction column to factor
    filtered_intensity_data[[selected_interaction]] <- factor(filtered_intensity_data[[selected_interaction]])

    # Format the p-values to 5 decimal places and color them
    filtered_posthoc <- filtered_posthoc %>%
      mutate(
        label = paste0("adj.p = ", formatC(adj.p, format = "f", digits = 5)),
        color = ifelse(adj.p < 0.05, "red", "black")  # Set color based on adj.p value
      )

    # Set y-position to the maximum value present in the value column
    max_value <- max(filtered_intensity_data$value, na.rm = TRUE) * 1.05

    # Create the box plot and add p-values
    p <- ggboxplot(filtered_intensity_data, x = selected_interaction, y = "value", fill = selected_interaction) +
      stat_pvalue_manual(
        filtered_posthoc,
        y.position = max_value, step.increase = 0.2,
        label = "label",
        color = "color"
      ) +
      scale_color_identity() +
      labs(
        x = paste(selected_metabolite, "\n", selected_interaction),
        y = "Normalized Intensity") +
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14), # Change x-axis tick label font size
        axis.text.y = element_text(size = 14), # Change y-axis tick label font size
        legend.position = "none", # Remove the legend
        plot.title = element_blank() # Remove the main title
      )

    message("Plot Posthoc Completed...")
    return(p)
  } else {
    message("Invalid input: Please ensure selected_metabolite, selected_interaction, posthoc data, and long intensity data are provided.")
  }
}