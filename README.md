# PhD - Advancing the measurement of cognitive ability: Developing a Cattell-Horn-Carroll computer adaptive screening test

This repository includes the data and code required for the cleaning and analyses in the PhD thesis.

## Abstract

The current thesis explored the application of a Computer Adaptive Test (CAT) in the measurement of intelligence, the efficiency and validity of a CAT when measuring intelligence across a range of ages and compared this newly developed CAT with the Wechsler Intelligence Scale for Children – Fifth Edition. 

Any measurement tool in psychology must be developed based on a good theoretical framework and sound measurement principles. The background of this thesis (Chapter 2) demonstrated that the Cattell-Horn-Carroll (CHC) theory is the most contemporary and suitable theory of intelligence to form the theoretical basis of a cognitive ability CAT. It is also important that researchers implement appropriate statistical models when developing and implementing a CAT, and thus the background of this thesis introduced many basic concepts related to different IRT models and characteristics of CAT that some readers may be unfamiliar with. A review of the literature demonstrated that despite CHC, IRT and CAT all being well known within the literature, there has been only limited integration of CHC with CAT in practice. Nearly all studies reviewed failed to describe the CAT used with enough detail for their studies to be replicated by researchers, or the CATs to be implemented by practicing psychologists. Additionally, none of the CATs investigated demonstrated utility with an Australian sample.

The four studies in this thesis investigated the use of four sets of items developed from the perspective of CHC theory, designed to measure Lexical Knowledge, Induction, Visualisation and Working Memory. Both the data and statistical analyses for all four studies included in this thesis are accessible at github.com/jakekraska/phd.
Early studies included an Item Tryout Study (ITOS; Chapter 3) and Item Calibration Study (ICS; Chapter 4). The initial ITOS demonstrated Lexical Knowledge items tended to be too easy, Induction items were not always predictable in their ordering, that there were possible unidimensionality issues with the Visualisation items, and there were design issues with the Working Memory item stimuli. The ICS took advantage of further item development and addressed problems with the item set as identified by the ITOS and the analysis resulted in retaining 47 Lexical Knowledge items, 23 Induction items, 30 Visualisation items, 25 Working Memory items. Item parameters were exported for subsequent CAT simulation.

Each item set was included in a CAT simulation (Chapter 5) that made use of 1,000 simulated participants across 52 simulations. It was found that that item administration could be reduced by approximately 50% and that the item sets best measured at very low ability levels through to average abilities. If such a test were implemented in practice, testing time would be approximately 18 minutes depending on the ability of the examinee. Classification of those with deficits in cognitive abilities could be achieved, however differentiation between those with average ability and above average ability may not be possible with the current psychometric outcomes. 

Analysis of the convergent validity of the items is investigated (Chapter 6). is mounting evidence of the psychometric validity of the Lexical Knowledge and Induction item sets. Further analysis of the Working Memory item sets in comparison to other cognitively complex Working Memory tasks and further development of Visualisation items that provide a better fit to the Rasch model are both recommended. 
A discussion (Chapter 7) of the implications and limitations of this research is also presented. Future research opportunities are identified surrounding multidimensional item response theory, improvements and standardisation in item stimuli via funded multimedia design, further item development, and implementation of items that measure other CHC abilities.

## Authors and Contributions

Jake Kraska wrote the thesis that this code was developed for. Supervision for this PhD was provided by Dr Shane Costello, Dr John Roodenburg and Associate Professor Wendy McKenzie.

## Correspondence

All correspondence regarding this analysis and thesis can be directed to [Jake Kraska](mailto:jake.kraska@monash.edu)

## Installation and Running

* Clone this repo to your local machine using https://github.com/jakekraska/phd
* Open each chapters R code in your preferred IDE, starting with the first chapter
* Set the working directory to the source file within the IDE or using `setwd(dirname(rstudioapi::getActiveDocumentContext($path)`
* Install the necessary packages using the `install.packages("PACKAGENAME")` command
* Each chapter may have used older versions of packages and thus you may need to specify the URL. For example:
    * `packageurl <- "http://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_0.9.1.tar.gz"`
    * `install.packages(packageurl, repos = NULL, type = "source")`
    * If you don't know the URL, you can look for it in the [CRAN binary archive](https://cran.r-project.org/src/contrib/Archive)
* The code is moderately commented and minimal SEM, IRT and CAT knowledge should allow understanding of the output. Some R knowledge is required to understand the code.