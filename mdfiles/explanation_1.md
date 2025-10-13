<font size="+1">

#### **Analyzing Datasets in HRonBC**
Multiple (up to 4) MCF7 and/or T47D RNA-seq datasets can be chosen for simultaneous differential gene expression analyses to compare samples of interest. It can be useful to check the *Data Explanation* tab to see details of each dataset included.
<br><br>Note that the app allows for the user to utilize the preferred tool (Limma-Voom/DESeq2) regarding its convenience for the desired downstream analyses and choose among its parameters to display as columns. A meta analysis will be conducted for the p-values across all selected datasets with the desired method and the resulting value for each gene will be added as an additional column. All meta analysis methods are provided by the [poolr](https://cran.r-project.org/web/packages/poolr/index.html) package. 
<br><br>**Example**
<br>HRonBC RNA-seq datasets of MR overexpressing MCF7 cells and T47D cells are chosen for identifying differentially expressed genes upon treatment with drospirenone (DRSP) only, and with drospirenone in combination with estradiol (DRSP+E2).
<br><br><center><img src="screenshot1-1.png" style="height: 263.5px; width:941.5px;")/></center>
<br>

-------------------------------

#### **Visualizing The Chosen Experiments**
When scrolled down the page, further insight into the attributes of samples within each of the chosen datasets are provided in separate tab panels. 
<br><br>The read counts of each sample can be checked by the *Library Size Bar Plot* and the deviation among biological replicates can be seen by the *Multidimensional Scaling Plot*. It is also possible to check out and download the raw count data from the tab beside.
<br><br>**Example**
</font>
<br><br><center><img src="screenshot1-2.png" style="height: 314.5px; width:942px;">
<br>(MCF7 experiment)
<br><img src="screenshot1-3.png" style="height: 298.5px; width:937.5px;">
<br>(T47D experiment) </center>