


# Load the necessary library
library(officer)

# Read the existing Word document
doc <- read_docx("Rshiny_Report_Template.docx")

# Modify the document
# Assuming the placeholders are in the text and need to be replaced
doc <- doc %>%
  body_replace_all_text("Project Name: [Enter Project Name]", "Project Name: Your Project Name", fixed = TRUE) %>%
  body_replace_all_text("Date: [Enter Date]", "Date: " %>% format(Sys.Date()), fixed = TRUE) %>%
  body_replace_all_text("Number of tables: xxx", "Number of tables: [Your Number Here]", fixed = TRUE) %>%
  body_replace_all_text("Number of Figures: xxx", "Number of Figures: [Your Number Here]", fixed = TRUE)

# Save the modified document
print(doc, target = "Modified_Rshiny_Report.docx")
