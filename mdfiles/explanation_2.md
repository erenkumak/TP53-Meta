<font size="+1">

#### **Filtering and Visualizing Differentially Expressed Genes**
In the *Results* section, genes can be filtered by any of the displayed columns. It is also possible to constrict the gene set to a certain pathway of interest among KEGG pathways for a specific analysis.
<br><br>All selected rows will be added up to the resulting gene list displayed down the page (previous selections will be stored upon a new filtration). The *select all rows* button helps store all the genes that suit the given filtrations. The number of genes that has been stored in the list can be checked by the selected row count that shows up below the table.
<br><br>**Example**
<br>Genes that are upregulated upon both DRSP and DRSP+E2 in the selected MCF7 samples are filtered out to have logFC values above 0.5, and p-values within the specified range of significance. 
<br><br><center><img src="screenshot2-1.png" style="height: 594px; width:1151.4px;"/></center>
<br><br> Several more sets of genes are stored for analysis by different filtrations that are all done with 0.5 logFC cutoff and within the 95% confidence interval: downregulated genes upon both DRSP and DRSP+E2 in MCF7 cells; up/downregulated genes upon both DRSP and DRSP+E2 in T47D cells; genes that are downregulated upon DRSP while being upregulated upon DRSP+E2 -and vice versa- in MCF7 and T47D cells separately; genes that are up/downregulated upon DRSP in both cell lines; genes that are up/downregulated upon DRSP+E2 in both cell lines. 
<br><br>The comparison of T47D and MR overexpressing MCF7 data can help distinguish between gene regulations associated with ER/PR/MR status and cell specific gene regulations upon the treatments. 
<br><br>The complete gene list is accumulated down the page as row selections are done after each filtration.
<br><br><center><img src="screenshot2-2.png" style="height: 2633.8px; width:1236.3px;"/></center>
<br><br> It is seen that a total of 637 genes has been selected. The list can be exported as a file for later uses as well.
<br><br><center><img src="screenshot2-3.png" style="height: 199.2x; width:1137px;"/></center>
<br>

-------------------------------

The volcano plots down the page provide a quick visualization of differentially expressed genes for all of the selected treatment groups, which might as well help with gene selection.
<br><br>**Example**
<br>Several genes such as *ZBTB16*, *PARM1*, *TUBA3E*, *AQP3*, and *SGK1* stand out when the volcano plots are compared. Also, it can be interpreted that the number of genes that exhibit significant differential expression is higher in the T47D experiment than in the MCF7 experiment. 
<br><br><center><img src="screenshot2-4.png" style="height: 719px; width:939px;"/></center>
<br>

-------------------------------

On the right, a scatter plot can be drawn to visualize how the experiments of interest correlate with one another in terms of the differential expression profiles of either the whole gene set or the filtered gene list.
<br><br>**Example**
<br><br><center><img src="screenshot2-5.png" style="height: 338px; width:1039px;"/></center>
</font>

<br>

-------------------------------

<font size="-1">
Please note that any example is for representation purposes only and is not conclusive.
</font>


