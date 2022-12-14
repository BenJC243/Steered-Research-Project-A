library(shiny)
install.packages(factoextra)
library(shinythemes)
library(Seurat)
library(SeuratObject)
install.packages("uwot")
library(uwot)
library(factoextra)
#setwd('/home/solomon/Research_project/saverds')

set_1 <- readRDS('pipeline1.rds')

data_list = list(set_1=set_1)

ui <- fluidPage(theme = shinytheme('slate'),
                #titlePanel('Steered Research Project Group A'),
                navbarPage('Steered Research Project Group A',
                           tabPanel('About',
                                    titlePanel('Single cell RNA Sequencing Data Analysis'),
                                    br(),
                                    div(includeMarkdown('about.md'),
                                    )),
                           tabPanel('Damaris et al 2015 - Pipeline 1',
                                    mainPanel(
                                      #withTags({
                                      br(),
                                      tags$h4('As sbset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others'),
                                      br(),
                                      br(),
                                      plotOutput('expression_plot'),
                                      br(),
                                      tags$h4('(A) BIC classification using different models identifies ten clusters of cells. EII corresponds to a spherical model with equal volume and equal shape, VII corresponds to a spherical model with variable volume and equal shape, EEI 
                                      corresponds to a diagonal model with equal volume and equal shape, VEI corresponds to a diagonal model with variable volume and equal shape, EVI corresponds to a diagonal model with equal volume and variable shape, VVI corresponds to a diagonal 
                                      model with variable volume and variable shape, EEE corresponds to a ellipsoidal model with equal volume and equal shape, EEV corresponds to a ellipsoidal model with equal volume and equal shape, VEV corresponds to a ellipsoidal model with 
                                      variable volume and equal shape and VVV corresponds to a ellipsoidal model with variable volume and variable shape (B) Uncertainty of classification within members of each cluster.'),
                                      br(),
                                      textOutput('BIC_text'),
                                      br(),
                                      br(),
                                      plotOutput('BIC', width = '100%'),
                                      br(),
                                      tags$h4('We observed one less cluster than was reported in the original publication (Darmanis et al, 2015; Fig S2.B), however our method returned a lesser degree of uncertainty across all clusters. This suggests that our replication of the original 
                                              method may have merged two clusters that were of questionable definition in the original publication.'),
                                      br(),
                                      plotOutput('darmanis_cluster', width = '100%') 
                                    )),
                           
                           
                           tabPanel('Alternative Analysis - Pipeline 2 ',
                                    tags$h1('Seurat stuff'),
                                    sidebarPanel('Gene Selection',
                                                 selectInput('var', 'Choose a gene for more information',
                                                             choices = c('CD24' = ' CD24 is a signal transducer that is belived to play a role in cell differentiation. It also promotes AG-dependent proliferation of B cells and prevents the terminal differentiation into antibody forming cells ', 'TMSB15A'= 'is involved in cytoskeleton organisation via binding to actin monomors and inhibiting actin polymerisation', 'SOX11'= ' SOX11 is a transcriptional activator of TEAD2, binds to gene promoters with POU3F2/BRN2 or PIU3F1/OCT6 and plays a role in cell survival in developing tiiue. More specifically in tissues such as the neural tube, branchial arches and somites',
                                                                         'DCX'= 'DCX is a microtubule-associated protein invovled in neuronal dispersion and cortex lamination during cerebral cortex development. May also be involved in the promotion of neuronal migration', 'S100A11' ='Allows for the differentiation and cornification of keratinocytes', 'SRGN' = 'Involved in the formation of mast cell secretory granule. Necessary for storage of some proteases in connective tissue and mucosal mast cells, and for storage of granzyme B in T-lymphocytes. Also involved in the localisation of neutrophil elastase in azurophil granules of neutrophils. Mediates processing of MMP2. Involved in cytotoxic cell granule-mediated apoptosis.', 
                                                                         'BZM' = 'This is gene B2M', 'ARHGDIB'  = ' Regulates the GDP/GTP exchange reaction of the Rho proteins by inhibiting the dissociation of GDP from them, and the subsequent binding of GTP to them. Also plays a role in actin mediated cytoskeleton development.', 
                                                                         'CD74' = 'Encodes for the gamma chain of the HLA class II histocombatibility antigen. Stabilises the complex and transports the complex from the endoplasmic reticulum to the endosomal/lysosomal system to allow for antigen processing. Also acts as a cell surface receptor for MIF', 
                                                                         'GABRA1'  = 'Is involved in the heteropentameric receptor for GABA, a major neurotransmitter inhibitor in the brain. Is also involved in the formation of inhibitory GABAergic synapses.', 'GABRG2'  = 'Is involved in the heteropentameric receptor for GABA, a major neurotransmitter inhibitor in the brain. Is also involved in the formation of inhibitory GABAergic synapses.', 
                                                                         'SNAP25' = 'Is involved in regulation of neurotransmitter release. Also involved in the synaptic function of some neuronal systems', 'TMEM130' = 'Transmembrane protein expressed in the brain', 
                                                                         'SYNPR' = 'An intrinsic membran protein involved in small synaptic vesicles', 'CKMT1B' = 'Reversibly catalyses phosphate transfer between ATP and phosphogens. Therefore plays an important part in energy transduction in tissues such as the brain', 
                                                                         'MEG3' = 'Maternally expressed RNA gene, with expression found in the brain. May play a role in apoptosis', 'UNC80' = 'Subunit of the NALCN sodium channel complex. This is a voltage gated ion channel that controls Na+ levels, in turn controlling neuronal excitability.', 
                                                                         'AHI1' = 'Involved in vesicle trafficing. Also involved in cillary signalling during cerebellum development in embroyonic de velopment', 
                                                                         'EEF1A1' = 'Promotes GTP-dependent binding of aminoacyl-tRNA to ribsosomes. Also plays a role in regulaiton of IFNG transcription. ', 'ACSBG1' = 'Ligase that catalyses the converstion of fatty acids into acyl-CoAs. Present in all regions of the brain excluding the pituitary gland', 
                                                                         'SLC39A12' = 'Zinc influx transporter protein expressed in the brain', 'AGXT2L1' = 'Enzyme that catalyses the pyridoxal-phosphate-dependent brakdown of phosphoethanolamine', 'GJB6' = 'Encodes for a gap junction.',
                                                                         'PRODH' = 'Causes the conversion of proline to delta-1-pyrroline-5-carboxylate.', 'MAG' = 'This is gene MAG', 'KLK6' = 'This is gene KLK6', '9JB1' = 'This is gene 9JB1', 'CNDP1' = 'This is gene CNDP1', 'CARN1' = 'This is gene CARN1',
                                                                         'PRODH' = 'Causes the conversion of proline to delta-1-pyrroline-5-carboxylate.', 'MAG' = 'Mediates interaction between myleinating cells and neurons by binding to neuronal sialic acid-containing gangliosides. Plays a role in axon myelinaiton maintenance.', 
                                                                         'KLK6' = 'Is a serine protease which shows activity against multiple proteins including: amyloid precursor protein, myelin basic protein, gelatin, casein, fibronectin, laminin, vitronectin and collagen. It is also thought to play a potential role in axon outgrowth after spinal cord injury.', 
                                                                         'CNDP1' = 'Enzyme which catalyses peptide bond hydrolysis in Xaa-His dipeptides. Found in the adult central nervous system.', 'CARNS1' = 'Enzymes that catalyses the synthesis of carnosine and homocarnosine.' , 'VSTM2A' = 'Plays a role in the regulation of the early stage of white and brown preadipocyte cell differentiation.', 'NLGN4Y' = ' Neuronal cell surface protein invovled in cell-cell interactions', 
                                                                         'DLX6' = ' Has roles including DNA-binding transcription factor activity, and sequence specific DNA binding. DLX6 is expressed at day 9 of embryonic development, with expression on day 12.5 being in the brain and bones. From this, expression levels decline from day 15. Primarily ecpressed in the forebrain.'),
                                                             selected = '', multiple = FALSE),
                                                 radioButtons('dataset', label = ('Pipeline 2 - Seurat'),
                                                              choices = list('sample1' = 'set_1'),
                                                              selected = 'set_1'),
                                                 #helpText('Enter gene names in CAPS.  Gene names must be exact.'),
                                                 #textInput('gene2', label = 'Gene Name', value = ""),
                                                 selectInput('gene2', 'Select top 5 expressed genes in each cluster',
                                                             choices = c('CD24', 'TMSB15A', 'SOX11','DCX' , 'S100A11', 'SRGN', 'CD74',    'GABRA1' , 'GABRG2' , 'SNAP25', 'TMEM130' ,  'SYNPR' , 'CKMT1B' ,  'MEG3', 'UNC80' ,   'AHI1',    'EEF1A1'  , 'ACSBG1',      'SLC39A12' , 'AGXT2L1' , 'GJB6' ,       
                                                                         'PRODH', 'MAG' , 'CNDP1', 'PRODH' ,'MAG' ,'KLK6',  'CNDP1', 'CARNS1' , 'VSTM2A', 'NLGN4Y','DLX6')
                                                 ),
                                                 
                                                 column(8,
                                                        column(6,
                                                               mainPanel(
                                                                 #h1("First level title", align = "center"),
                                                                 br(),
                                                                 textOutput('dropdown'), #tableOutput('tabledata'),
                                                                 textOutput('firstgene'),
                                                                 plotOutput('dimPlot1', width = '1500%'),
                                                                 br(),
                                                                 tags$h4('We observed a total of 7 clusters seen here, commpared to 10 identified by Darmanis et al (2015).'),
                                                                 br(),
                                                                 br(),
                                                                 br(),
                                                                 plotOutput('genePlot2', width = '1500%'),
                                                                 #plotOutput('dimPlot3', width =  '2000%'),
                                                                 br(),
                                                                 tags$h4('This graph allows for the selection of individual genes within the cluster. The genes we idnetified can be found in the corresponding report.'),
                                                                 br(),
                                                               )
                                                        )
                                                 )
                                    )
                           ),
                           
                           tabPanel("References",
                                    div(includeMarkdown("references.md"))
                           )
                )
)








server <- function(input, output, session) {
  
  output$firstgene <- renderText(input$var)
  output$dropdown <- renderText(input$barCD24)
  #output$expression_text <- renderText(input = "subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly 
  #expressed in some cells, and lowly expressed in others'")
  
  output$txtout <- renderText({
    paste( input$txt1, input$txt2, sep = ' ' )
  })
  
  
  datasetInput <- reactive({
    df <- data_list[[input$dataset]]
  })
  
  output$dimPlot1 <- renderPlot({
    DimPlot(pbmc, reduction = 'umap', label =T)
  })
  
  output$darmanis_cluster <- renderPlot({
    plot(fit, what= 'classification')  
  })
  
  output$BIC <- renderPlot({
    plot(fit, what= 'BIC')
  })
  
  output$expression_plot <- renderPlot({
    plot4
  })
  
  output$genePlot2 <- renderPlot({
    FeaturePlot(pbmc, features = (input$gene2), reduction = 'umap')
  })
  
  # output$dimPlot3 <- renderPlot({
  #   Dimplot(pbmc, reduction = 'umap')
  # })
}

shinyApp(ui, server)  