# ROAD: RNA Origami Automated Design

The ROAD software package automates the design of RNA nanostructures using the RNA origami method and consists of the following main programs:

1.	**RNAbuild** - generates 3D PDB coordinate files from a blueprint .
1.	**Revolvr** - repeatedly mutates a randomized RNA starting sequence until it folds into a target 2D structure under the ViennaRNA folding model.
1.	**Batch_revolvr** - runs Revolvr multiple times to create different candidate designs for the same blueprint.
1.	**RNAvis** - analyzes and animates one potential folding path from a blueprint.
1.	**flip_trace** - outputs a blueprint in all four possible orientations.
1.	**trace_analysis** - parses a blueprint and outputs analyses and helpful variants of the blueprint. 
1.	**trace_pattern** - parses a blueprint and outputs the dot-paren notation for the structure in the blueprint, as well as the sequence constraints. It creates an input for Revolvr.
1.	**trace** - threads a sequence onto a blueprint. Only used by Revolvr.

## Citing ROAD

When using this software in publications please cite:

[Geary, Grossi, McRae, Rothemund & Andersen; *RNA origami design tools enable cotranscriptional folding of kilobase-sized nanoscaffolds*; Nature Chemistry 2021](https://bion.au.dk/publications/)  

## Installation requirements:
To run RNAbuild you will need to first install a Perl environment.
Next, to run Revolvr, you will additionally need to install the ViennaRNA package.
Finally, you will need a plaintext editor that can copy/paste blocks of fixed-width text (using so-called “block” or “column” mode) such as TextEdit on MacOS, or Notepad++ and Vim on Windows.

This software was written and tested with Perl v5.18.2 for MacOS, and ViennaRNA versions 2.1.9 and 2.4.10. Installation for Perl and ViennaRNA should be a few minutes. Software was tested on a Mid 2012 Macbook Pro running El Capitan 10.11.6 and a Mid 2015 Macbook Pro running Sierra 10.12.6.

## Video Tutorials
[Part 1: installation](https://www.youtube.com/watch?v=lvTj_qRppME)
