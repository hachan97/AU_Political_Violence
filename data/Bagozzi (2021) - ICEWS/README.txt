REPLICATION INFORMATION

These replication files allow one to replicate the figures and results reported in Bagozzi, Berliner, & Welch's
	"The Diversity of Repression: Measuring State Repressive Repertoires with Events Data," 
	as well as the figures and tables reported in its Supplemental Appendix.

This folder includes:

1. A "JPR_replication_code_final.R" script, which (1) formats all relevant variables from the raw data, (2) produces all relevant plots for the main article and
	the Supplemental Appendix, and (3) produces all relevant tables reported in the Supplemental Appendix. Note:
	this R script loads and utilizes the following five R packages: e1071, entropy, foreign, sandwich, and stargazer.

2. A "JPR_countrydata.csv" datataset, which corresponds to the country-year covariates that are used in Tables A.II-A.V in the Supplemental Appendix. Note: the 
	"JPR_replication_code_final.R" script mentioned above reads-in this csv file and merges these country-year covariates to our primary event data measures.

3. An "eventdata_icews_threedigit.dta" Stata dataset, which corresponds to the repression-oriented ICEWS event data inputs that are used for all primary analyses.
	Note: "JPR_replication_code_final.R" script reads-in this dta file and formats these raw repression-oriented event data inputs into our final variables.

4. An "eventdata_icews_threedigitNOFILTER" Stata dataset, which corresponds to the repression-oriented ICEWS event data inputs that have been created without using
	one-a-day filtering. Note: "JPR_replication_code_final.R" script reads-in this dta file and formats these raw repression-oriented event data inputs 
	into our final variable for the first model in Table A.III.
 
5. An "eventdata_icews_threedigit_publisher_1995" Stata dataset, which corresponds to the repression-oriented ICEWS event data inputs that have been created 
	for only ICEWS events that were coded from news sources that had been coded by ICEWS since the its very first year of coverage (1995). 
	Note: "JPR_replication_code_final.R" script reads-in this dta file and formats these raw repression-oriented event data inputs 
	into our final variable for the final model in Table A.III.

6. An "eventdata_icews_goldstein" Stata dataset, which corresponds to the Goldstein-scored repression-oriented ICEWS event data inputs that are used for the various
	Goldstein-score based assessments in the Supplemental Appendix.

7. A "Diversity_of_Repression_Appendix.pdf" file, which corresponds to the final Supplemental Appendix for our article.
 

Please contact Daniel Berliner (danberliner@gmail.com) or Benjamin E. Bagozzi (bagozzib@udel.edu) for any questions about these data or files.

