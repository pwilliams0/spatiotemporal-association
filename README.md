# spatiotemporal-association
## Author
Peter Williams
## About
This project analyses fine-scale spatiotemporal associations among terrestrial vertebrates in Malaysian Borneo. This approach consists of calculating the time differences between the detections of one species (the inducer) and detections of another species (the responder), simulating what the time differences would be if the two species ignored each other, and comparing the observed and expected distributions of time differences to determine whether the responder species is associating with the inducer species

There are 3 sections in the R code: **Load Data for analyses**, **Single test of spatiotemporal association**, and **Pairwise comparisions of spatiotemporal association for list of species**.

* **Load Data for analyses** loads required libraries and the two attached CSV files. These files are camera trap capture data (Captures.csv) and a table of the periods in which the cameras were active (CT_periods.csv).

* **Single test of spatiotemporal association** analyses the spatiotemporal response of one species (responder) to the presence of one other species (inducer). Any species can be chosen as responder or inducer. The roles are for interpreting results -- they are not inherently biologically meaningful. The script compiles the observed distribution of time differences, simulates the expected distribution of time differences, and compares the two distrubtions. The distributions used are trucanted Weibull distributions. The two parameters of a Weibull distribution are k and lambda. The output from this analysis includes: the number of observed responder detections after inducer within truncation limit, p-value for comparing observed vs. expected, k observed, SE of k observed, lambda observed, SE of lambda observed, k expected, SE of k expected, lambda expected, SE of lambda expected, and the % change in k comparing observed vs. expected. P-value and % in k are especially valuable for interpreting results.

* **Pairwise comparisions of spatiotemporal association for list of species** performs the same analysis as the above section, but repeats the analysis for all possible pairwise combinations of a list of species. This includes all combinations of inducer/responder pairing for all species. The ouput from this analysis is the same as the single test, but results are saved as a table of values rather than single values. In these tables, row names are the responder, and column names are the inducer.

## Variables
The following variables may be changed for different versions of the analyses:
### Single test of spatiotemporal association
* *resp* (line 42) - Responder species, set to "sambar"
* *indu* (line 43) - Inducer species, set to "bearded_pig"
* *Periods* (line 46) - Subset of data used in analysis (e.g. mast years, logged forest, etc.), set to all data
* *N_sim* (line 76) - Number of simulated observations of responder, set to 50,000
* *trunc* (line 152) - Truncation point for analyses, set to 14 days
### Pairwise comparisions of spatiotemporal association for list of species
* *N_sim* (line 195) - Number of simulated observations of responder, set to 50,000
* *trunc* (line 198) - Truncation point for analyses, set to 14 days
* *Periods_subset* (line 239) - Subset of data used in analysis (e.g. mast years, logged forest, etc.), set to all data
* *species_list* (line 245) - List of species to include in pairwise comparisons, set to c("bearded_pig", "mousedeer", "yellow_muntjac", "sambar", "argus", "pig_tailed_macaque", "fireback", "malay_civet", "banded_civet")
