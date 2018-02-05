# bmotif: counting motifs in bipartite networks

## Overview

`bmotif` is software to count occurrences of motifs in bipartite networks, as well as the number of times each node appears in each unique position within motifs. `bmotif` considers all 44 unique bipartite motifs up to six nodes and all 148 unique positions within these motifs. It is available in [R](https://github.com/SimmonsBI/bmotif), [MATLAB](https://github.com/SimmonsBI/bmotif-matlab) and [Python](https://github.com/SimmonsBI/bmotif-python). `bmotif` was originally developed to analyse bipartite species interaction networks in ecology but its methods are general and can be applied to any bipartite graph. The code is released under the MIT license.

## Scripts
### User-facing functions
- **Main.m**: Shows how to use the two main functions on a simple example network
- **motifs.m**: Counts the number of times motifs occur in a network
- **Position_motifs.m**: Counts the number of times nodes occur in unique positions within motifs

### Internal functions
- **Check_mot.m**: Internal function for checking that the input arguments for motifs.m are valid
- **Check_pos.m**: Internal function for checking that the input arguments for Position_motifs.m are valid
- **tensor_make.m**: Internal function for making tensor
- **tensorR.m**: Internal function for calculating the tensor product of arrays
