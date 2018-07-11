Create a txt file similar to "input data.txt" 
Save it in the C:/R/R3.1/library/estimate/extdata/input data.txt


# Choose directory to save results.

Run code: 

library(estimate)

CancerExpr <- system.file("extdata", "input data.txt", package="estimate")

filterCommonGenes(input.f=CancerExpr, output.f="genes.gct", id="GeneSymbol")
	
estimateScore("genes.gct", "estimate_score.gct", platform="affymetrix")

plotPurity(scores="estimate_score.gct", samples=" GSM447614.CEL", platform="affymetrix")

Results: 
Open estimate_score.gct in a excel sheet to view scores for individual samples.
Total tumour purity is generated as a plot in the estimated_purity_plots folder.
