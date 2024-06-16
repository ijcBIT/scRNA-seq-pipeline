

library(shiny)
library(xfun)
library(dplyr)
library(hdf5r)
library(Seurat)
library(ggplot2)
library(ggraph)
library(RColorBrewer)
library(data.table)
library(SingleR)
library(celldex)
library(shinythemes)
library(shinyWidgets)

options(shiny.maxRequestSize = 50 * 1024^2)

ui <- fluidPage(
  theme = shinytheme("flatly"),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
    tags$style(HTML("
      .scrollable-sidebar {
        height: 90vh;
        max-height: 90vh;
        overflow-y: auto;
        padding: 10px 15px;
        background-color: #f7f7f7;
        border-right: 1px solid #ddd;
      }
    "))
  ),
  titlePanel("scRNAseq Analysis Pipeline"),
  sidebarLayout(
    sidebarPanel(
      class = "scrollable-sidebar",
      width = 3,
      div(class = "sidebar-header", "Control Panel"),
      fileInput("dataFile", "Choose the .h5 file", accept = c(".h5")),
      actionButton("loadData", "Load Data", icon = icon("upload")),
      hr(),
      h4("QC and Filtering"),
      sliderInput("nUMI", "Minimum nUMI:", min = 0, max = 10000, value = 500),
      sliderInput("nGene", "Minimum nGene:", min = 0, max = 10000, value = 250),
      sliderInput("log10GenesPerUMI", "Minimum log10GenesPerUMI:", min = 0, max = 5, value = 0.8),
      sliderInput("mitoRatio", "Maximum Mito Ratio:", min = 0, max = 1, value = 0.2),
      actionButton("filterData", "Filter Data", icon = icon("filter")),
      hr(),
      h4("Normalization"),
      actionButton("normalizeData", "Normalize Data", icon = icon("balance-scale")),
      hr(),
      h4("Cell Cycle Scoring"),
      actionButton("cellCycleScoring", "Score Cell Cycle", icon = icon("sync")),
      hr(),
      h4("Variable Features"),
      actionButton("variableFeatures", "Find Variable Features", icon = icon("search")),
      hr(),
      h4("Scaling the Data"),
      actionButton("scaleData", "Scale Data", icon = icon("arrows-alt")),
      hr(),
      h4("Dimensional Reduction"),
      actionButton("runPCA", "Run PCA", icon = icon("chart-bar")),
      actionButton("runUMAP", "Run UMAP", icon = icon("project-diagram")),
      actionButton("runTSNE", "Run tSNE", icon = icon("braille")),
      hr(),
      h4("Clustering"),
      actionButton("findNeighbors", "Find Neighbors", icon = icon("sitemap")),
      actionButton("findClusters", "Find Clusters", icon = icon("object-group")),
      hr(),
      h4("Marker Identification"),
      actionButton("findMarkers", "Find Markers", icon = icon("map-marker-alt")),
      hr(),
      h4("Cell Type Annotation"),
      actionButton("cellTypeAnnotation", "Annotate Cell Types", icon = icon("tags"))
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Data Summary", verbatimTextOutput("summary")),
        tabPanel("QC Plots", plotOutput("qcCorrelationPlot"), plotOutput("qcViolinPlot")),
        tabPanel("PCA", plotOutput("pcaPlot"), plotOutput("pcaLoadingsPlot")),
        tabPanel("UMAP", plotOutput("umapPlot")),
        tabPanel("tSNE", plotOutput("tsnePlot")),
        tabPanel("Markers", verbatimTextOutput("markers"), plotOutput("markersPlot")),
        tabPanel("Scaling Data", verbatimTextOutput("scalingDescription"), plotOutput("variableFeaturesPlot"), plotOutput("labelPointsPlot")),
        tabPanel("Annotation", plotOutput("annotationPlot")),
        tabPanel("PCA Heatmaps", plotOutput("pcaHeatmap1"), plotOutput("pcaHeatmap2"))
      )
    )
  )
)

server <- function(input, output, session) {
  pbmc.data <- eventReactive(input$loadData, {
    req(input$dataFile)
    Read10X_h5(input$dataFile$datapath)
  })
  
  pbmc <- eventReactive(input$loadData, {
    req(pbmc.data())
    CreateSeuratObject(counts = pbmc.data(), project = "scRNA_seq", min.cells = 3, min.features = 100)
  })
  
  output$summary <- renderPrint({
    req(pbmc())
    pbmc()
  })
  
  observeEvent(input$filterData, {
    req(pbmc())
    seurat_obj <- pbmc()
    seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
    seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-") / 100
    seurat_obj@meta.data <- seurat_obj@meta.data %>%
      dplyr::rename(seq_folder = orig.ident, nUMI = nCount_RNA, nGene = nFeature_RNA)
    
    cells_before <- ncol(seurat_obj)
    
    filtered_obj <- tryCatch(
      subset(seurat_obj, subset = (nUMI >= input$nUMI) & 
               (nGene >= input$nGene) & 
               (log10GenesPerUMI > input$log10GenesPerUMI) & 
               (mitoRatio < input$mitoRatio)),
      error = function(e) {
        showNotification("No cells found after filtering. Please adjust the filtering criteria.", type = "error")
        return(NULL)
      }
    )
    
    if (!is.null(filtered_obj)) {
      cells_after <- ncol(filtered_obj)
      pbmc <<- reactive(filtered_obj)
      output$summary <- renderPrint({ pbmc() })
      showNotification(paste("Filtered from", cells_before, "to", cells_after, "cells."), type = "message")
    }
  })
  
  observeEvent(input$normalizeData, {
    req(pbmc())
    seurat_obj <- NormalizeData(pbmc(), normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
    pbmc <<- reactive(seurat_obj)
    output$summary <- renderPrint({ pbmc() })
  })
  
  observeEvent(input$cellCycleScoring, {
    req(pbmc())
    s.genes <- cc.genes.updated.2019$s.genes
    g2m.genes <- cc.genes.updated.2019$g2m.genes
    seurat_obj <- CellCycleScoring(pbmc(), s.features = s.genes, g2m.features = g2m.genes)
    pbmc <<- reactive(seurat_obj)
    output$summary <- renderPrint({ table(pbmc()$Phase) })
  })
  
  observeEvent(input$variableFeatures, {
    req(pbmc())
    seurat_obj <- FindVariableFeatures(pbmc(), selection.method = "vst", nfeatures = 2000)
    pbmc <<- reactive(seurat_obj)
    output$summary <- renderPrint({ pbmc() })
  })
  
  observeEvent(input$scaleData, {
    req(pbmc())
    seurat_obj <- ScaleData(pbmc())
    pbmc <<- reactive(seurat_obj)
    
    output$scalingDescription <- renderText({
      "Scaling the data involves removing unwanted sources of variation; e.g. technical variation, batch effects, biological variations due to cell cycle state, mitochondrial contamination.
      
      Linear transformation ('scaling'), a standard pre-processing step prior to dimensional reduction techniques like PCA. The 'ScaleData' function:
      · Shifts the expression of each gene, so that the mean expression across cells is 0.
      · Scales the expression of each gene, so that the variance across cells is 1. This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate."
    })
    
    output$variableFeaturesPlot <- renderPlot({
      var_features <- VariableFeaturePlot(pbmc())
      var_features
    })
    
    output$labelPointsPlot <- renderPlot({
      top10 <- head(VariableFeatures(pbmc()), 10)
      var_features <- VariableFeaturePlot(pbmc())
      lab_points <- LabelPoints(plot = var_features, points = top10, repel = TRUE)
      lab_points
    })
  })
  
  observeEvent(input$runPCA, {
    req(pbmc())
    seurat_obj <- RunPCA(pbmc(), features = VariableFeatures(object = pbmc()))
    pbmc <<- reactive(seurat_obj)
    output$pcaPlot <- renderPlot({ DimPlot(pbmc(), reduction = "pca", group.by = "Phase", split.by = "Phase") })
    output$pcaLoadingsPlot <- renderPlot({
      VizDimLoadings(pbmc(), dims = 1:2, reduction = "pca")
    })
    output$pcaHeatmap1 <- renderPlot({
      DimHeatmap(pbmc(), dims = 1, cells = 500, balanced = TRUE)
    })
    output$pcaHeatmap2 <- renderPlot({
      DimHeatmap(pbmc(), dims = 1:15, cells = 500, balanced = TRUE)
    })
  })
  
  observeEvent(input$runUMAP, {
    req(pbmc())
    seurat_obj <- RunUMAP(pbmc(), dims = 1:20)
    pbmc <<- reactive(seurat_obj)
    output$umapPlot <- renderPlot({ DimPlot(pbmc(), reduction = "umap", label = TRUE) })
  })
  
  observeEvent(input$runTSNE, {
    req(pbmc())
    seurat_obj <- RunTSNE(pbmc(), dims = 1:20)
    pbmc <<- reactive(seurat_obj)
    output$tsnePlot <- renderPlot({ TSNEPlot(object = pbmc(), label = TRUE) })
  })
  
  observeEvent(input$findNeighbors, {
    req(pbmc())
    seurat_obj <- FindNeighbors(pbmc(), dims = 1:20)
    pbmc <<- reactive(seurat_obj)
    output$summary <- renderPrint({ pbmc() })
  })
  
  observeEvent(input$findClusters, {
    req(pbmc())
    seurat_obj <- FindClusters(pbmc(), resolution = 0.8)
    pbmc <<- reactive(seurat_obj)
    output$summary <- renderPrint({ pbmc() })
    output$qcPlots <- renderPlot({ DimPlot(pbmc(), group.by = "RNA_snn_res.0.8", label = TRUE) })
  })
  
  observeEvent(input$findMarkers, {
    req(pbmc())
    all.markers <- FindAllMarkers(pbmc(), only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
    output$markers <- renderPrint({ table(all.markers$cluster) })
    top3_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
    output$markersPlot <- renderPlot({ DotPlot(pbmc(), features = unique(top3_markers$gene)) + coord_flip() })
  })
  
  observeEvent(input$cellTypeAnnotation, {
    req(pbmc())
    monaco.ref <- celldex::MonacoImmuneData()
    sce <- as.SingleCellExperiment(DietSeurat(pbmc()))
    monaco.main <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.main)
    monaco.fine <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.fine)
    seurat_obj <- pbmc()
    seurat_obj@meta.data$monaco.main <- monaco.main$pruned.labels
    seurat_obj@meta.data$monaco.fine <- monaco.fine$pruned.labels
    pbmc <<- reactive(seurat_obj)
    srat <- SetIdent(pbmc(), value = "monaco.fine")
    output$annotationPlot <- renderPlot({ DimPlot(srat, label = TRUE, repel = TRUE, label.size = 3) + NoLegend() })
  })
  
  # QC Plots
  output$qcCorrelationPlot <- renderPlot({
    req(pbmc())
    metadata <- pbmc()@meta.data
    ggplot(metadata, aes(x = nUMI, y = nGene, color = mitoRatio)) + 
      geom_point() + 
      scale_colour_gradient(low = "gray90", high = "black") +
      stat_smooth(method = lm) +
      scale_x_log10() + 
      scale_y_log10() + 
      theme_classic() +
      geom_vline(xintercept = 500) +
      geom_hline(yintercept = 250)
  })
  
  output$qcViolinPlot <- renderPlot({
    req(pbmc())
    VlnPlot(object = pbmc(),
            features = c("nUMI", "nGene", "mitoRatio"),
            pt.size = 0.01, ncol = 3) &
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
            axis.title.x = element_blank())
  })
  
  # PCA Visualization
  output$pcaViz <- renderPlot({
    req(pbmc())
    VizDimLoadings(pbmc(), dims = 1:2, reduction = "pca")
  })
  
  # Non-linear Dimensional Reduction Plots
  output$umapPlot <- renderPlot({
    req(pbmc())
    DimPlot(pbmc(), reduction = "umap", label = TRUE)
  })
  
  output$tsnePlot <- renderPlot({
    req(pbmc())
    TSNEPlot(object = pbmc(), label = TRUE)
  })
  
  # Cluster Segregation by "nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio"
  output$featurePlot <- renderPlot({
    req(pbmc())
    metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")
    FeaturePlot(pbmc(), 
                reduction = "umap", 
                features = metrics,
                pt.size = 0.4, 
                order = TRUE,
                min.cutoff = 'q10',
                label = TRUE)
  })
  
  # Boxplot of nGene per cluster
  output$boxplot <- renderPlot({
    req(pbmc())
    ggplot(pbmc()@meta.data) +
      geom_boxplot(aes(x = RNA_snn_res.0.8, y = nGene, fill = RNA_snn_res.0.8)) +
      NoLegend()
  })
  
  # DotPlot for Marker Genes
  output$markerDotPlot <- renderPlot({
    req(pbmc())
    knownGenes <- c("CD34","CRHBP","GATA1", "CD14","IRF8","CD19",
                    "CD4","CD8B","GNLY")
    DotPlot(pbmc(), features = knownGenes) + coord_flip()
  })
  
  # Violin Plot for Marker Genes
  output$markerViolinPlot <- renderPlot({
    req(pbmc())
    knownGenes <- c("CD34","CRHBP","GATA1", "CD14","IRF8","CD19",
                    "CD4","CD8B","GNLY")
    VlnPlot(pbmc(), features = knownGenes, stack = TRUE, flip = TRUE)
  })
}

shinyApp(ui, server)

