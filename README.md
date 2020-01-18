OverallChirality
================

This code will compute the chirality described in (A. Pietropaolo,
L. Muccioli, R. Berardi and C. Zannoni, Proteins 70, 667, 2008),
extending its application from backbone atoms of proteins to any
sequence of beads. It will read an xyz file and two parameters need to
be defined: a cutoff distance and the number of atoms to be taken into
account.


BUILD (Linux)
-------------
g++ -lm walking_chirality.cpp -o WalkingChirality


USAGE
-----
Usage:
        ./WalkingChirality <xyz> <cutoff> <Na>

The code will compute an average chirality value and write two
files. The two columns in "chirality_index.txt" contain the id of the
atom and its own chiral value. Plotting this file, the chirality trend
along the system can be followed. The "chirality_index.dump" file can
be viewed with [Ovito](www.ovito.org) software and beads can be
colored with their own chirality index value.


EXAMPLE
-------
        ./WalkingChirality conical_spiral.xyz 15 15

A snapshot of the conical spiral colored with its chirality index and
the output plot are attached.