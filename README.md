Streamlining Reverse Transcription-Quantitative Polymerase Chain Reaction Results

This is my first coding project that streamlines RT-qPCR results. Understanding the values produced by RT-qPCR software can be difficult. In this project, I generated a 
python script that would take these values and produce the meaningful results users wanted from their data.

This script was created  for an undergradaute lab course that would teach the basics of RT-qPCR. In this lab, E. coli was treated with hydrogen peroxide, or treated with 
water. Two reference genes (housekeeping genes) and the gene of interest were utilized for this lab. After being treated, reverse transcription-quantitative polymerase chain
reaction was performed and their relative gene expressions were calculated with this script. Below is the equation utilized for calculation:

Relative Gene Expression = 

                    (E of gene of interest) ^ deltaCT GOI
                    -------------------------------------
                    Geometric Mean[(E Ref) ^ deltaCT Ref]

E = Primer efficiency
GOI = Gene of interest
Ref = Reference genes

The information on how to utilize these statistics was produced by Steven Bradburn, PhD.
Link to his website: https://toptipbio.com/qpcr-multiple-reference-genes/
