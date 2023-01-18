# PollyCommonR

## Guide to PollyCommonR


```R
library(PollyCommonR)
```

### Read Raw Data


```R
raw_data <- read.csv("data/demo_raw_intensity.csv", stringsAsFactors = FALSE, check.names = FALSE)
```

### Read Metadata


```R
metadata <- read.csv("data/demo_metadata.csv", stringsAsFactors = FALSE)
```

### Create Sample Matrix


```R
raw_data$id <- paste0(raw_data$groupId, "_", raw_data$compound)
sample_raw_mat <- sample_intensity_matrix(raw_intensity_df = raw_data, metadata_df = metadata,rownames_col = "id")

# Shifting data by 1
sample_raw_mat <- sample_raw_mat + 1
```

### COV calculation


```R
cov_cal_df <- calculate_cohortwise_cov(raw_matrix = sample_raw_mat, metadata = metadata, cohort_col = "Cohort")
```

```R
p <- create_cohortwise_cov_boxplot(calculated_cov_df = cov_cal_df, interactive = FALSE)
```

```R
p <- create_cohortwise_cov_barplot(calculated_cov_df = cov_cal_df,id_order = '1087_Std-L-Phenylalanine', id_col = 'id')
```
 
### Pre Normalization

#### Boxplot on Raw Data


```R
p <- create_boxplot_on_matrix(sample_raw_mat = sample_raw_mat, x_label = "Sample",y_label = "Raw Intensity",title_label = "Boxplot on Raw Data")
``` 

#### Density Plot on Raw Data


```R
p <- create_densityplot_on_matrix(sample_raw_mat)
```

#### PCA on Raw Data


```R
pca_compute <- compute_pca(sample_raw_mat = sample_raw_mat)
```

```R
p <- plot_proportion_of_variance(PCAObj_Summary = pca_compute)
```

```R
p <- plot_pca(pca_compute, metadata = metadata, condition = 'Cohort', title_label = "PCA on Raw Data", interactive = FALSE)
```

```R
p <- plot_pca3d(pca_compute, metadata = metadata, condition = 'Cohort', pc_x = 1, pc_y = 2, pc_z = 3, title_label = "PCA on Raw Data")
```

### Post Normalization

#### Normalization by Internal Standard


```R
norm_agent <- t(sample_raw_mat["1087_Std-L-Phenylalanine",])
```


```R
norm_mat <- normalize_by_scaling_factor(sample_raw_mat, normalization_agent = norm_agent, scaling_factor_col = 1)
```


#### Log2 Transformation


```R
log2_norm_mat <- log2(norm_mat)
log2_norm_mat_shift <- max(abs(log2_norm_mat)) + log2_norm_mat
```

#### Boxplot on Normalized Data


```R
p <- create_boxplot_on_matrix(log2_norm_mat_shift, x_label = "Sample",y_label = "Normalized Intensity",title_label = "Boxplot on Normalized Data")
```


#### Density Plot on Normalized Data


```R
p <- create_densityplot_on_matrix(log2_norm_mat_shift)
```


#### PCA on Normalized Data


```R
norm_pca_compute <- compute_pca(log2_norm_mat_shift)
```
 

```R
p <- plot_proportion_of_variance(PCAObj_Summary = norm_pca_compute)
```


```R
p <- plot_pca(norm_pca_compute,metadata = metadata, condition = 'Cohort', title_label = "PCA on Normalized Data", interactive = FALSE)
```

```R
p <- plot_pca3d(norm_pca_compute, metadata = metadata, condition = 'Cohort', pc_x = 1, pc_y = 2, pc_z = 3, title_label = "PCA on Normalized Data")
```


#### PCA on Filtered Normalized Data


```R
filtered_metadata <- filter_metadata_by_cohorts(metadata, condition = "Cohort", selected_cohorts = c("Cohort_2", "Cohort_5", "Cohort_3"))
filtered_log2_norm_mat_shift <- log2_norm_mat_shift[,filtered_metadata[,1]]
```
  
```R
filtered_pca_compute <- compute_pca(filtered_log2_norm_mat_shift)
``` 

```R
p <- plot_pca(filtered_pca_compute, metadata = metadata, condition = 'Cohort', title_label = "PCA on Normalized Data", interactive = FALSE)
```


#### Samplewise bar plot for single metabolite


```R
p <- create_samplewise_barplot(log2_norm_mat_shift, metadata, id_name = "2_GMP", cohort_col = "Cohort",
                               x_label = "Sample", y_label = "Normalized Intensity",  title_label = "GMP")
```


### Differential Expression Analysis


```R
diff_exp <- compute_differential_expression(norm_data = log2_norm_mat_shift, metadata = metadata, cohort_col = 'Cohort', cohort_a = "Cohort_2", cohort_b = "Cohort_1", algo = "limma")


```


```R
p <- plot_volcano_from_limma(diff_exp = diff_exp, log2fc_range = 0, p_val_cutoff = 0.05, interactive = FALSE)
```
 
### Perform Anova Test


```R
log2_norm_mat_shift_anova <- compute_anova(norm_data = log2_norm_mat_shift, metadata = metadata, cohort_col = "Cohort")
```
    