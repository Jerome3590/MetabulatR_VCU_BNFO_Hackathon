#title: "Metabulator"
#author: "BNFOThon Team 3 - Cassidy Coates, Riley Fanus, Jerome Dixon, Mike Hall"
#date: "4/2/2022"


library(shiny)
library(paws)
library(aws.signature)
library(dplyr)
library(rlist)
library(data.table)
library(tidyr)
library(kableExtra)
library(DT)
library(here)
library(httr)
library(readr)
library(stringr)
library(jsonlite)
library(rvest)
library(reshape2)



run_query <- function(query_input) {
  
  
  
  query_input = str_replace_all(query_input, " ", "%20")
  r <- GET(paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", query_input, "/synonyms/TXT", sep = ""))
  # Save to file
  bin <- content(r, "text")
  # Data cleanup
  bin = str_replace_all(bin, ";", "\n")
  writeBin(bin, "Data/pubchem.txt")
  
  # Read pubchem synonyms as data frame
  PubchemSynonyms_df = read_delim("Data/pubchem.txt", delim = "\n", col_names = 'Synonyms')
  
  
  
  # Get the PubChem cid
  cid <- GET(paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", query_input, "/cids/TXT", sep = ""))
  bin = content(cid, "text")
  cid = str_replace_all(bin, "\n", "")
  
  # Get the CAS ID from Pubchem
  cas_get <- GET(paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/", query_input, "/json", sep = ""))
  # Save to file
  bin <- content(cas_get, "raw")
  writeBin(bin, "Data/pubchem_CAS.json")
  # Access the CAS value from the JSON
  # From an experimental sample, it appears to always be item 1 under 0->synonyms
  cas_json <- read_json("Data/pubchem_CAS.json")
  CAS <- (cas_json[["PC_Substances"]][[1]][["synonyms"]][[2]])
  
  # Get the molecular formula from pubchem
  pubchem_formula <- str_trim(as.character(GET(paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", query_input, "/property/MolecularFormula/txt", sep = ""))))
  
  # Wipe the data storage files
  if(file.exists("Data/kegg.txt")) {
    file.remove("Data/kegg.txt")
  }
  file.create("Data/kegg.txt")
  
  if(file.exists("Data/mbw.txt")) {
    file.remove("Data/mbw.txt")
  }
  file.create("Data/mbw.txt")
  
  if(file.exists("Data/metacyc.txt")) {
    file.remove("Data/metacyc.txt")
  }
  file.create("Data/metacyc.txt")
  
  # We don't need to iterate over all of the databases - many we can get the data from just one api call (mostly using Pubchem cid)
  # HMDB 
  # (requires the use of the CAS identifier obtained from Pubchem, see above in code)
  # Read the table from the HMDB database page
  hmdb_page <- read_html(paste("https://hmdb.ca/unearth/q?utf8=%E2%9C%93&query=", CAS, "&searcher=metabolites", sep = ""))
  hmdb_table <- html_table(html_element(hmdb_page, css="table"))
  # Extract the the HMDB ID
  hmdb_id = as.character(hmdb_table[5, 2])
  
  # Get the molecular formula from HMDB
  hmdb_formula <- str_trim(as.character(hmdb_table[33, 2]))
  
  # Metabolomics Workbench
  mbw <- GET(paste("https://www.metabolomicsworkbench.org/rest/compound/pubchem_cid/", cid, "/regno/text", sep= "")) 
  # Save to file
  bin <- content(mbw, "text")
  # Write
  cat(bin, file="Data/mbw.txt", append=TRUE)
  mbw_regno = as.character(read_delim("Data/mbw.txt", delim = "\t", col_names = TRUE)[2])
  
  # Get the molecular formula from Metabolomics Workbench
  mbw <- GET(paste("https://www.metabolomicsworkbench.org/rest/compound/pubchem_cid/", cid, "/formula/text", sep= "")) 
  # Save to file
  bin <- content(mbw, "text")
  # Write
  cat(bin, file="Data/mbw_formula.txt", append=FALSE)
  mbw_formula = as.character(read_delim("Data/mbw_formula.txt", delim = "\t", col_names = TRUE)[2])
  
  # MetaCyc (pathways)
  metacyc <- GET(paste("https://websvc.biocyc.org/META/foreignid?ids=Pubchem:", cid, sep= "")) 
  # Save to file
  bin <- content(metacyc, "text")
  # Write
  cat(bin, file="Data/metacyc.txt", append=TRUE)
  meta_name = as.character(read_delim("Data/metacyc.txt", delim = "\t", col_names = FALSE)[3])
  
  keggID = -1
  
  # Iterate over the synonyms
  keggID = -1
  cat("1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 12; 13; 14; 15; 16; 17; 18; 19; 20\n", file="Data/kegg.txt", append=FALSE)
  
  # Iterate over the synonyms
  for (row in 1:(nrow(PubchemSynonyms_df) - 1)) {
    #for (row in 1:2) {
    name = as.character(PubchemSynonyms_df[row, "Synonyms"])
    
    # Send API requests for all secondary databases
    # Metabolite chemical databases
    
    # KEGG
    kegg <- GET(paste("http://rest.kegg.jp/find/compound/", str_replace_all(name, " ", "%20"), sep= "")) 
    # Save to file
    bin <- content(kegg, "text")
    # Data cleanup
    bin <- str_replace_all(bin, "\t", "; ")
    # Write
    cat(bin, file="Data/kegg.txt", append=TRUE)
    
    # If we had an exact match, set the kegg ID to it
    kegg___df = fread(file="Data/kegg.txt", sep = ";", sep2 = "\n", fill = TRUE, data.table=FALSE)
    kegg___df[is.na(kegg___df)] <- 0
    for (row2 in 1:nrow(kegg___df)) {
      for (col in 1:ncol(kegg___df)) {
        if(tolower(name) == tolower(as.character(kegg___df[row2,col]))) {
          keggID = substring(tolower(as.character(kegg___df[row2, 1])), first = 5, last = 10)
          print(name)
          break;
        }
      }
    }
    
    if(keggID != -1) {
      break;
    }
    
    
    # TODO LIPID MAPS
    
    # Articles
    # TODO PubMed
    
    # Sleep, let the remote servers recover
    Sys.sleep(0.33)
  }
  
  # Push everything to the DynamoDB table
  svc <- dynamodb(
    config = list(
      endpoint = "https://dynamodb.us-east-1.amazonaws.com",
      region = "us-east-1"
    )
  )
  svc$batch_write_item(RequestItems = list(
    Metabulator = list(
      list(
        PutRequest = list(
          Item = list(
            Search_Term = list(
              S = query_input
            ),
            Query_Result = list(
              S = paste("Pubchem, ", cid, ", https://pubchem.ncbi.nlm.nih.gov/compound/", cid, ", ", pubchem_formula, sep = "")
            )
          )
        )
      ),
      list(
        PutRequest = list(
          Item = list(
            Search_Term = list(
              S = query_input
            ),
            Query_Result = list(
              S = paste("HMDB, ", hmdb_id,", https://hmdb.ca/metabolites/", hmdb_id, ", ", hmdb_formula, sep = "")
            )
          )
        )
      ),
      list(
        PutRequest = list(
          Item = list(
            Search_Term = list(
              S = query_input
            ),
            Query_Result = list(
              S = paste("KEGG, ", keggID,", https://www.kegg.jp/entry/", keggID, ", ", "NO_FORMULA", sep = "")
            )
          )
        )
      ),
      list(
        PutRequest = list(
          Item = list(
            Search_Term = list(
              S = query_input
            ),
            Query_Result = list(
              S = paste("Metabolomics Workbench, ", mbw_regno,", https://www.metabolomicsworkbench.org/data/StructureData.php?RegNo=", mbw_regno, ", ", mbw_formula, sep = "")
            )
          )
        )
      ),
      list(
        PutRequest = list(
          Item = list(
            Search_Term = list(
              S = query_input
            ),
            Query_Result = list(
              S = paste("MetaCyc, ", meta_name,", https://metacyc.org/compound?orgid=META&id=", meta_name, ", ", "NO_FORMULA", sep = "")
            )
          )
        )
      )
    )
  )
  )
}



svc <- paws::dynamodb()


# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("MetabolatR"),
  
  textInput('input_text', 'Enter Metabolite for Search:', value="ethanol"),
  
  verbatimTextOutput('output'),
  mainPanel(
  tableOutput("items_kable")
  )
 )


# Define server logic required to draw a histogram
server <- function(input, output) {

  output$items_kable <- function() {
    query_input <- as.character(input$input_text)
    
    run_query(query_input)
    
    dynamoDB <- 
      svc$scan(
        TableName = "Metabulator"
      )
    
    items <- bind_rows(list.clean(dynamoDB$Items, function(x) length(x) == 0L, recursive = TRUE))
    
    items %>%
      filter(Search_Term == tolower(query_input)) %>% 
      knitr::kable("html") %>%
      kable_styling("striped", full_width = F) 
    }
}

# Run the application 
shinyApp(ui = ui, server = server)
