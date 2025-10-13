<font size="+1">
Genome-wide chromatin interactions can be investigated by the Hi-C square heatmap for which two distinct data can be selected as upper and bottom matrices separated by the diagonal where chromosomal positions are mapped along. The physical interactions between genomic regions are scored by colour intensities. Higher intensities indicate that the regions are spatially in close proximity such as in chromatin loops. 
Self interactions can be seen on the diagonal that generally comprises of high colour values unless there is contrary data available. The color values on the sides of the diagonal provide insight into topologically associating domains (TADs) which have high interaction frequencies. Examining TAD boundaries among different samples can help interpret their roles in gene regulation by altering promoter-enhancer interactions.
<br><br>Analyzing Hi-C data along with the RNA-seq results might reveal the possible roles of the steroid hormone receptors ER, PR and MR as inducers of changes in the spatial genome organization apart from their roles as sequence-specific transcription factors.
<br><br>Note that the matrix data should be chosen equally with a maximum of 3 for each.
<br><br>**Example**
<br>Several genes were selected from the RNA-seq analysis results regarding their logFC valaues and their biological relevance to visualize their genomic regions over Hi-C squares formed by *in situ* and capture Hi-C data of control and E2 treated T47D samples. All control samples were selected for the upper matrices.
<br><br><center><img src="screenshot9-1.png" style="height: 1307.5px; width:949.5px;"/></center>
<br><br>It might provide a better understanding of TADs to increase the resolution and the distance to the genes of interest and visualize a larger region.
<br><br><center><img src="screenshot9-2.png" style="height: 1664px; width:1093.3px;"/></center>
<br>The second Hi-C squares show that *PGR* and *SGK1* are located in TAD boundaries whereas *ZBTB16* is licated within a contact enriched region.

-------------------------------

New selections can be done for the experimental groups to compare the expression patterns of selected genes in other samples. 
<center><img src="screenshot9-3.png" style="height: 822px; width:951px;"/></center>
<br><br>Colour codes for logFC values can be preferentially changed. The default palette is the same as in the *ChIP* section. &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <img src="screenshot9-4.png" style="height: 222.5px; width:118.5px;"/>
</font>
<br>

-------------------------------

<font size="-1">
Please note that any example is for representation purposes only and is not conclusive.
</font>

