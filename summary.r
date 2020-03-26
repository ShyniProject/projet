# Read data
summary = tabItem(tabName = "summary",
  
  h1("Summary"),
  
  ## FILE INPUT
  fileInput("data", label = NULL,
    buttonLabel = "Browse...",
    placeholder = "No file selected"),
  
  ## pvalue
  sliderInput("pv", label = h3("p-value threshold : "), min = 0, 
    max = 0.1, value = 0.05),
  textOutput("text"),  
  
  h2("Volcano plot")
)