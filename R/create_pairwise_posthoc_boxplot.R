#' create_pairwise_posthoc_boxplot
#'
#' Creates a pairwise post-hoc boxplot for selected metabolite and interaction.
#'
#' @param posthoc_data_long_df A dataframe containing long format post-hoc data.
#' @param intensity_data A dataframe containing intensity data.
#' @param selected_metabolite The selected metabolite for analysis.
#' @param selected_interaction The selected interaction for analysis.
#' @param filter_pvalues Logical indicating whether to filter comparisons with adj. p-value >= 0.05.
#' @param plot_axis_rotation Logical indicating whether to rotate axis labels and wrap text.
#' @return A ggplot2 object representing the pairwise post-hoc boxplot.
#' @import dplyr
#' @import plotly
#' @import ggpubr
#' @import stringr
#' @export
create_pairwise_posthoc_boxplot <- function(posthoc_data_long_df = NULL, intensity_data = NULL, selected_metabolite = NULL, selected_interaction = NULL, filter_pvalues = FALSE, plot_axis_rotation = FALSE) {
  if (selected_metabolite != "" && selected_interaction != "" && !is.null(posthoc_data_long_df) && !is.null(intensity_data)) {
    message("Plot Posthoc Started...")
    require(dplyr)
    require(plotly)
    require(ggpubr)
    require(stringr)

    # Filter the posthoc_result dataframe based on the selected metabolite and interaction
    filtered_posthoc <- dplyr::filter(posthoc_data_long_df, ID == selected_metabolite & Interaction %in% selected_interaction)

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

    # Reshape the data to long format
    long_intensity_data <- tidyr::pivot_longer(filtered_intensity_data, 
                                              cols = all_of(selected_interaction), 
                                              names_to = "Interaction", 
                                              values_to = "Groups")

    # Create a custom sorting order based on selected_interactions
    interaction_order <- factor(selected_interaction, levels = unique(selected_interaction))

    # Reorder the levels of the Interaction column
    long_intensity_data$Interaction <- factor(long_intensity_data$Interaction, levels = interaction_order)

    # Select relevant columns and arrange the data
    filtered_intensity_data <- long_intensity_data %>% 
      dplyr::select(id, value, Interaction, Groups) %>%
      arrange(Interaction, Groups) # Arrange by Interaction & Groups

    # Filter out p-values >= 0.05 if filter_pvalues is TRUE
    if (filter_pvalues) {
      filtered_posthoc <- filtered_posthoc %>% dplyr::filter(adj.p < 0.05)
    }

    # Format the p-values to 5 decimal places and color them
    filtered_posthoc <- filtered_posthoc %>%
      mutate(
        label = paste0("adj.p = ", adj.p),
        color = ifelse(adj.p < 0.05, "red", "black")  # Set color based on adj.p value
      )

    # Set y-position to the maximum value present in the value column
    max_value <- max(filtered_intensity_data$value, na.rm = TRUE) * 1.05

    # Adjust axis labels and rotation based on plot_axis_rotation
    if (plot_axis_rotation) {
      angle_value <- 45
      hjust <- 1
      vjust <- 1
      labels_fun <- function(x) stringr::str_wrap(x, width = 5)
    } else {
      angle_value <- 0
      hjust <- 0.5
      vjust <- 0.5
      labels_fun <- waiver()
    }

    # Create the box plot and add p-values
    p <- ggboxplot(filtered_intensity_data, x = "Groups", y = "value", fill = "Groups") +
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
      scale_x_discrete(labels = labels_fun) + # Apply wrapping function if rotation is enabled
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 10, angle = angle_value, hjust = hjust, vjust = vjust),  # Rotate x-axis labels, Change x-axis tick label font size
        axis.text.y = element_text(size = 10), # Change y-axis tick label font size
        legend.position = "none", # Remove the legend
        plot.title = element_blank() # Remove the main title
      )

    message("Plot Posthoc Completed...")
    return(p)
  } else {
    message("Invalid input: Please ensure selected_metabolite, selected_interaction, posthoc data, and long intensity data are provided.")
  }
}