# bmotif: counting motifs in bipartite networks

## Overview

`bmotif` is software to count occurrences of motifs in bipartite networks, as well as the number of times each node appears in each unique position within motifs. `bmotif` considers all 44 unique bipartite motifs up to six nodes and all 148 unique positions within these motifs. It is available in [R](https://github.com/SimmonsBI/bmotif-release), [MATLAB](https://github.com/SimmonsBI/bmotif-matlab) and [Python](https://github.com/SimmonsBI/bmotif-python). `bmotif` was originally developed to analyse bipartite species interaction networks in ecology but its methods are general and can be applied to any bipartite graph.

## Motif and motif position dictionary
The motifs corresponding to each motif ID and the positions corresponding to each motif position ID can be found in **Simmons, B. I., Sweering, M. J. M., Dicks, L. V., Sutherland W. J., Di Clemente, R. bmotif: a package for counting motifs in bipartite networks. bioRxiv. doi: 10.1101/302356**. These were defined following Baker et al (2015) Appendix 1 Figure A27.

## Scripts
### Guide
- **Main.m**: Shows how to use the two main functions on a simple example network

### User-facing functions
- **motifs.m**: Counts the number of times motifs occur in a network
- **Position_motifs.m**: Counts the number of times nodes occur in unique positions within motifs

### Internal functions
- **Check_mot.m**: Internal function for checking that the input arguments for motifs.m are valid; called by motifs.m
- **Check_pos.m**: Internal function for checking that the input arguments for Position_motifs.m are valid; called by Position_motifs.m
- **tensor_make.m**: Internal function for making tensor; called by motifs.m and Position_motifs.m
- **tensorR.m**: Internal function for calculating the tensor product of arrays; called by motifs.m and Position_motifs.m

## License
The code is released under the MIT license (see LICENSE file).

## Citation
If you use the package in your work, please cite:
Simmons, B. I., Sweering, M. J. M., Dicks, L. V., Sutherland W. J., Di Clemente, R. bmotif: a package for counting motifs in bipartite networks. bioRxiv. doi: 10.1101/302356

## References
Baker, N.J., Kaartinen, R., Roslin, T. and Stouffer, D.B., 2015. Species’ roles in food webs show fidelity across a highly variable oak forest. Ecography, 38(2), pp.130-139.
