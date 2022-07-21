% Morphospace-Disparity Analysis  
%===========================================================================================
% MDA: Morphospace-Disparity Analysis Version 1.2.
%   Examines occupational patterns of morphospace using disparity metrics.
%   Estimates of metrics are based on bootstrap resampling or rarefied bootstrap resampling.
%   Three temporal analysis which extract disparity pattern are available, and one performs
%   a bootstrap tail test.
%
% MDA: main function
%   User interface for the choice of data and type of analysis
%
%   Usage: mda in the Matlab workspace without any argument
%   Necessary arguments are asked using user-interface during program execution according
%   to needs.
%
%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
%   SGA: Single group analysis.
%      Morphospace are subdivided by bin contains in occurrence file
%      10 metrics of disparity are available:
%      Range-based:
%            Sum of ranges:             Rarefied Sum of univariates ranges
%            Product of ranges:         Rarefied N Root of product of univariates ranges
%            Range:                     Rarefied Maximal Euclidean pairwise distance
%            Area CHull:                Rarefied Area of the convex hull 
%                                       (calculate only on the two firsts dimensions)
%            PCO Volume:                Product of the two largest eigenvalues of
%                                       the cross-distance matrix, divided by (N)^2
%                               
%      Variance-based:
%            Sum of variances:          Sum of univariates variances
%            Product of variances:      N Root of product of univariates variances
%            Mean Pairwise Distance:    Mean Euclidean pairwise distance
%            Median Pairwise Distance:  Median Euclidean pairwise distance
%            Mean Distance Centroid:    Mean Euclidean distance species and subspace 
%                                       centroid
%      
%      2 metrics of location are available:
%            Minima:                    Rarefied Minimum on each component
%            Maxima:                    Rarefied maximum on each component
%
%-------------------------------------------------------------------------------------------
%   MGA: Multiple groups analysis.
%      Similar analysis and metrics than UGA but previously morphospace   
%      is divided following data contains in group file.
%
%-------------------------------------------------------------------------------------------
%   PDA: Partial disparity analysis.
%      Overall disparity in t-time is partitionned into groups (contain in group file):
%      Use particular metric:
%            Partial Disparity:          Sum of square Euclidean distance between species 
%                                        and the overall centroid in the subspace 
%
%-------------------------------------------------------------------------------------------
%   Bootstrap resampling: 
%      For these previous analyses, estimates of metrics correspond to the mean of the 
%      bootstrap samples. Number of resampling is a user-choice
%           Standard deviation:          On bootstrap samples
%           Confindence interval:        Can be performed based on percentile method with an 
%                                        user-choice of the level of confidence
%
%-------------------------------------------------------------------------------------------
%   Rarefaction: 
%      Range-based metrics are very sensitive to sample size, for supress this bias a 
%      rarefaction can be used. Rarefaction is perform using bootstrapping.
%
%-------------------------------------------------------------------------------------------
%   BTailTest: Bootstrap Tail Test.
%      All variance-based metrics plus partial disparity are available
%      Compare observed disparity of groups contain in group file to the bootstrap 
%      distribution of the metric choose for the first group
%      Probabilities done correspond to the two-tail, upper one-tail and lower one-tail case
%      These probabilities equal to Number of bootstrap replicates > and/or < to the 
%      observed value divided by the number of boostrap resampling
%
%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
% MDA_prefs       - contain default name of input files and format
% MDA_version     - contain address author, version number, copyright, form of citation
%                   display in the workspace when mda run.
% Modif_version   - Major modifications from previous version
%  
%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
% MDA_1D: one-dimensional data version 
%        This function is restricted to 10 variables 
%        Just SGA and BTailTest function are available
%        For SGA:
%                Mean, Min, Max, Variance, Range of each variable are calculated
%
%        For BTailTest:
%                Analysis is performed on variance estimated
%
%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
%IntrMDA: Mat file for the introductory image 
%
%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
%Help folder contains:
%        1. A pdf file 'A user guide to MDA' explaining format of input data, analyses, options
%        ouput file, figures and examples of the text user-interface
%        2. Two folders of examples: Ex_Input_files and Ex_Output_files
%        This two folders correspond to the data and results obtain from simulation and use in 
%        the introductory paper of the program in Computers & Geosciences
%        3. The randclade folder contains function for stochastic simulation of the 
%        evolution of clade with extraction of subclades and alteration of fossil record. 
%        This function was used for created previous example data.
%
%===========================================================================================

