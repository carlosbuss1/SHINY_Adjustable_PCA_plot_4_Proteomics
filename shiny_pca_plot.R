# Adjustable PCA Analysis with Shiny

## Load Required Libraries
library(shiny)
library(ggplot2)
library(ggrepel)
library(limma)
library(edgeR)

# Load Data
file_path <- "proteomics_counts.csv"  # Update file path if necessary
data <- read.csv(file_path, header = TRUE, sep = "\t")

# Preprocess Data
sample_mapping <- c(
  'KO1' = 'KO1', 'KO2' = 'KO2', 'KO3' = 'KO3', 'KO4' = 'KO4', 'KO5' = 'KO5',
  'WT1' = 'WT1', 'WT2' = 'WT2', 'WT3' = 'WT3', 'WT4' = 'WT4', 'WT5' = 'WT5'
)
colnames(data)[2:ncol(data)] <- sample_mapping[colnames(data)[2:ncol(data)]]

# Create Groups
samples <- c('KO1', 'KO2', 'KO3', 'KO4', 'KO5', 'WT1', 'WT2', 'WT3', 'WT4', 'WT5')
group <- factor(c(rep('KO', 5), rep('WT', 5)), levels = c('WT', 'KO'))
numeric_data <- as.matrix(data[, samples])
rownames(numeric_data) <- data$Gene.Symbol

# Replace Missing Values
for (i in 1:ncol(numeric_data)) {
  min_val <- min(numeric_data[, i], na.rm = TRUE)
  numeric_data[is.na(numeric_data[, i]), i] <- min_val
}

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Define UI
ui <- fluidPage(
  titlePanel("PCA Analysis for Normalized Counts"),
  sidebarLayout(
    sidebarPanel(
      actionButton("run_pca", "Run PCA"),
      selectInput("normalization", "Normalization Method", choices = c("Voom", "Log CPM")),
      numericInput("min_expression", "Minimum Expression Threshold", value = 0.5, step = 0.1),
      checkboxInput("filter_low", "Filter Low Variance Genes", value = TRUE),
      sliderInput("low_variation_threshold", 
                  "Low Variation Filter (Quantile):", 
                  min = 0, max = 1, value = 0.25, step = 0.05)
    ),
    mainPanel(
      plotOutput("pcaPlot"),
      verbatimTextOutput("summaryText")
    )
  )
)

# Define Server Logic
server <- function(input, output, session) {
  # Reactive Preprocessing
  processed_counts <- reactive({
    req(input$run_pca)
    dge <- DGEList(counts = numeric_data)
    dge <- calcNormFactors(dge)

    # Normalize Data
    if (input$normalization == "Voom") {
      v <- voom(dge, design, plot = FALSE)
      counts <- v$E
    } else if (input$normalization == "Log CPM") {
      counts <- cpm(dge, log = TRUE)
    } else {
      counts <- numeric_data
    }

    # Filter by Minimum Expression
    counts <- counts[rowMeans(counts) >= input$min_expression, ]

    # Optionally Filter Low-Variance Genes
    if (input$filter_low) {
      variances <- apply(counts, 1, var)
      quantile_threshold <- input$low_variation_threshold
      counts <- counts[variances > quantile(variances, quantile_threshold), ]
    }

    return(counts)
  })

  # PCA Analysis
  pca_results <- reactive({
    counts <- processed_counts()
    pca <- prcomp(t(counts), scale. = TRUE)
    explained_variance <- round((pca$sdev^2 / sum(pca$sdev^2)) * 100, 2)
    list(pca = pca, explained_variance = explained_variance)
  })

  # Render PCA Plot
  output$pcaPlot <- renderPlot({
    req(pca_results())
    pca <- pca_results()$pca
    explained_variance <- pca_results()$explained_variance
    pca_data <- as.data.frame(pca$x)
    pca_data$Group <- group

    ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
      geom_point(size = 4, alpha = 0.8) +
      geom_polygon(data = pca_data, aes(x = PC1, y = PC2, fill = Group, group = Group), alpha = 0.2) +
      geom_text_repel(aes(label = rownames(pca_data)), size = 3, color = "black") +
      scale_color_manual(values = c("WT" = "royalblue", "KO" = "red3")) +
      scale_fill_manual(values = c("WT" = "lightblue", "KO" = "pink")) +
      labs(
        title = "PCA of Processed Counts",
        x = paste0("PC1: ", explained_variance[1], "% Variance"),
        y = paste0("PC2: ", explained_variance[2], "% Variance")
      ) +
      theme_minimal(base_size = 15) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
      )
  })

  # PCA Summary
  output$summaryText <- renderPrint({
    req(pca_results())
    summary(pca_results()$pca)
  })
}

# Run Shiny Application
shinyApp(ui = ui, server = server)

