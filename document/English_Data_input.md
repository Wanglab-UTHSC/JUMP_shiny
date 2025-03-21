---
title: "Steps for Data Import"
output: 
  html_document:
     css: ../www/style/style.css
---

# Exploratory Analysis

**Exploratory Analysis** provides a user-friendly interface for uploading and visualizing datasets. It summarizes key dataset characteristics, helping users understand data distribution and identify underlying patterns. This section provides effective quality control of your data.

---

## Steps for Exploratory Analysis

1. **Navigate to Exploratory Analysis Tab**

    Click on the `Exploratory Analysis` tab in the left sidebar of this page.

    ![Data Import Tab](../www/images/import_tab.png){width="90%"}

2. **Upload Data**

    Users can seamlessly upload a proteomic dataset, which should be in `tab-delimited text` or `csv` file format. Please ensure your file follows the correct format shown below.
    <mark>For a large dataset (over 30 Mb), please use JUMP Shiny locally.</mark>

    In addition, please click `Download example input expression table` for an example illustration.  
    
    The input data should be organized into three columns: Protein Accession Number (i.e., UniProt), Gene Name, and Protein Description, followed by as many samples as needed. JUMP Shiny supports `raw` abundance values as well as `log2` conversion values.   
    <mark>Note: if your data is already log2 transformed, please choose the `log2` option.</mark>  

    An example of input data is shown below:
         
    ![Jumpshiny](../www/images/example_jumpshiny_input.png){width=90%}

    If successfully uploaded, the input data will be displayed in the `Protein expression table` panel (see below).

    ![Import Data](../www/images/example_input.png){width=90%}



---

**JUMP Shiny Format:**

If your data contains Accession number, Gene Name, Description, and samples, use `jumpshiny` to upload. **The first column is required.**



---

**JUMPq Result Format:**

If your data sample starts from the 24th column, use `jumpq` to upload. Please remove the header rows of jumpq data.

![Jumpq](../www/images/example_jumpq_input.png){width=90%}

---

**JUMPq Batch Result Format:**

If your data contains batch info, use `jump_batch` to upload.

![Jumpq Batch](../www/images/example_jumpqbatch_input.png){width=90%}


3. **Group Assignment**

    After loading the dataset, input your grouping in the `[Meta information]` panel.

    ![Group Assignment](../www/images/group_info.png){width="30%"}  
  
    **Group Information File**  
    The Group Information File is required for the data analysis. It should adhere to the following structure:  
    
    - <ins>Sample Name Column</ins>: The first column should contain your sample names. These names must exactly match the corresponding column names in your input expression table. Only the columns specified in this file will be included in the analysis, so ensure that they are correct, complete, and matched.  
    
    - <ins>Grouping Name Column</ins>: After the sample name column, you can include one or more grouping columns. Each grouping column can represent different categories or factors relevant to your analysis (e.g., "control" vs. "treatment," "male" vs. "female," etc.). You can add as many grouping columns as necessary to capture all relevant grouping factors.  
  
    *_Note_*: Headers are required in this file to clearly identify each column. The header of the first column is required to be named as `Sample` or `sample`.
    
    Below is an example of `Group Information` file:  
 
    ![Group Info](../www/images/group_info2.png){width="30%"}  

4. **Confirm and Analyze**

    Click the `[Assign group information]` button and wait for the `[Summary]` section to display additional information about your dataset. You can download and save the plots in .svg format for further analysis or publication. All plots can be zoomed in and out for a closer examination of the data. 

    **Intensity Distribution Plot**  

    By clicking the `[Intensity Distribution]` tab, you can view box plots for all uploaded samples, with each group highlighted in different colors. You can filter out proteins with low intensity and customize the title, X-axis, and Y-axis labels as needed.

    ![Distribution](../www/images/example_distribution.png){width="90%"} 
    


    **PCA Plot**  

    The PCA Plot visualizes the distribution of selected groups based on Principal Component Analysis (PCA). This plot helps in identifying patterns and trends in your dataset by reducing the dimensionality and highlighting the differences and similarities between groups. We include 2D and 3D PCA plots. Each point in the plot represents a sample, and the position of the points indicates their relative similarity or difference based on the principal components. The axes represent the first two/three principal components, which capture the most variance in the data. This visualization can be useful for identifying outliers, clusters, and potential relationships between groups.  

    You can define the number of top variable proteins included in the PCA. The results may vary depending on the number of proteins selected.   Additionally, you can toggle the buttons to display or hide labels on the plot.  

    ![PCA Plot](../www/images/example_PCA.png){width="90%"}


    **Sample Correlation**

    The sample correlation heatmap visualizes the correlation within and between groups. This method organizes samples and features into a hierarchical tree, known as a dendrogram, based on their similarity or dissimilarity. The heatmap uses color gradients to represent the intensity of the correlation, with closely related samples or groups appearing closer together on the dendrogram. This visualization can help identify clusters of similar samples, reveal patterns in the data, and highlight differences between groups. It's a valuable tool for understanding the relationships and structure within your dataset.  

    Similarly, you can select the number of proteins to include in the cluster analysis. The percentage indicates the ratio of the selected top variable proteins to the total number of proteins in the dataset. You can also choose from various agglomeration and distance methods to customize the clustering process. These options allow you to refine the analysis and tailor the clustering approach based on the characteristics of your data and your specific research needs. 

    ![Heatmap](../www/images/example_heatmap.png){width="90%"}   



    **Group Selection**
  
    Navigate to the `[Group selection]` panel to explore different ways of grouping your data for visualization. You can select variables or categories to group your data, such as experimental conditions, genes, or sexes. This flexibility allows you to customize the visualization and highlight specific aspects of your data for more detailed analysis. Use the options in the panel to easily switch between different groupings and gain insights from various perspectives.  

    ![Group Selection](../www/images/group_select.png){width="40%"}  
