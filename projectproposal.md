# project proposal


## generel topics
- weekly meetings on Wednesdays, more meetings if necessary


## analysis steps
### part 1 - describing data, reprobuccibility, normalization
- describing structure of data frame
- normalization of data: the amount of each protein should be 100 for every replicate (e.g. amount of protein 1 in all 25 fractions of Ctrl_Rep1 should be 100)
- plotting distribution of random proteins to show data
- reproducibility: why is it important to show reproducibility
- prove reproducibility with qq-plots and correlation between replicates
- last step: mean values of the three replicates 
### part 2 - find maxima for each protein in Ctrl and RNase
- find global maxima -> highest value across all fractions
- find local maxima -> first criteria: each value of the two neighbouring pairs has to be lower as the local maxima (avoid fluctuation bias)
- find local maxima -> second criteria: value of local maxima has to be greater than 4% of total protein amount (100) (avoid background noise)
- maxima for each protein stored in an R object (e.g. data frame)
- criteria for maxima could change during project work 
### part 3 - differences of Ctrl and RNase maxima and selection criteria for R-DeeP
- evalute differences between Ctrl and RNase condition maxima with statistical test (e.g. t-test) using the three replicates
- selection criteria 1: p-value of statistical test < 0.05
- selection criteria 2: shift > 1 fraction
- selection criteria 3: 
### part 4 - apply slection criteria and clustering/PCA
- apply selection criteria to identify R-DeeP 
- group proteins with the same shifting behaviour via kmeans (e.g. left shift, right shift, no shift)
### part 5 - comparison with data bases and further analysis
- comparing our results with R-DeeP or RBPs of data bases
- linear regression (?)
- if we find new candidates for R-DeeP, search for their cellular functions and connection to diseases 


## project proposal presentation
- start presentation with example of RBPs or R-DeeP linked to disease
- describe RBPs / R-DeeP in generel and how our data was obtained
- describe our data
- present the 5 analysis steps including milestones and deliverables 
- approx. timeline for our project


## questions
- What does "describing the dataset" mean?
- structure of the report
- how to evalute differences of maxima between Ctrl and RNase condition
- which maxima shifted (global or local) and does a shift of a local maximum mean R-DeeP?
- why are there more than one peak for each condition? -> possible answer: different 3D-orientation
- why are there left shifts and right shifts? (komplex should always be in higher fraction than a single protein)
- regression analysis? modell for identifying false positive R-Deep?
- Can we upload our data on GitHub? (we should treat it confidential)
