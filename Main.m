% set seed for random matrix
seed=9575;
rng(seed);

% generate 10x25 random binary matrix
M = randi([0 1], 10,25);

% Specify the motif positions to calculate. Here we request all 148
% positions
V_motifs=1:148;

% Calculate the number of times each node occurs in each requested position
% within motifs. The output is a table with 148 rows, one for each
% unique position within motifs. There are also 36 columns: the first column 
%(IDmotif) gives the ID of the motif position, the next 10 columns 
% (A1 - A10) give the frequency with which each row node occurs
% in each motif position, the final 25 columns (B1 - B25) 
% give the frequency with which each column node occurs in each motif 
%position. NaN is returned when a node cannot occur in a given position 
% because that position can only be occupied by nodes in the other level. 
% For example, position 1 (IDmotif 1) can only be occupied by column nodes 
% and therefore a count of NaN is returned for all row nodes.
Motifp = Position_motifs(M,V_motifs)

% Specify the motifs to calculate. Here we request all 44 motifs
T_motifs=1:44;

% Calculate the number of times each requested motif occurs in the network.
% The output is a table with one row for each motif. The first column 
% (MotifID) gives the ID of the motif, the second column gives the 
% frequency with which that motif occurs in the network.
MotifT = motifs(M,T_motifs)
