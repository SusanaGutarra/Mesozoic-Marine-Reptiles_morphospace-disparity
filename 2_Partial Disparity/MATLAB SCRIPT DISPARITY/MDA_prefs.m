% MDA_prefs
%------------------------------------------------------------------------
% Initialize global variable for loading data
%    If files not exist in MDA directory, a file open
%    dialog box is open for select directory and file
%    Default names:
%      Data_matrix.txt:     coordinates of species in morphospace
%      line = species; columns = axes; tab-delimited without heading
%      ex:     1.50   3.67 ....
%             -2.43   1.76 ....
%              ....
%         Analysis: All
%
%      eigenvalues.txt:     array with eigenvalues 
%                           [lamb1 lamb2 lamb3 ...]
%         Can be necessary only with data coming from PCA 
%         on correlation matrix where variance on axes equal to 1, 
%         rescaling axes on eigenvalues permit to recover initial 
%         part of variance
%         Analysis: All if necessary
%
%      occurrence.txt:      binary matrix of presence/absence 
%                           of species in temporal level
%         line = species; columns = level; 
%         tab-delimited without heading
%         ex:      0 1 1 0 ....
%                  1 1 1 0 ....
%                  ....
%         Analysis: SGA, MGA, PDA
%         => level can be other nature that temporal 
%             (e.g., geographical, taxinomic)
%
%      group.txt:           binary matrix of presence/absence of species 
%                           in group
%         line = species; columns = group; 
%         tab-delimited without heading
%         ex:      0 1 1 0 ....
%                  1 1 1 0 ....
%                  .... 
%         Analysis: MGA, PDA, BTailTest
%------------------------------------------------------------------------

global filename                      % Default name of input file (scores on principal components)
global Eigenfile                     % Default name of input file (eigenvalues for standardization of PCs)
global Occurrencefile                % Default name of input file (temporal or geographical...occurrences)
global Groupfile                     % Default name of input file (group member)

filename='data_matrix.txt'; 
Eigenfile='eigenvalues.txt';
Occurrencefile='occurrence.txt';
Groupfile='group.txt';