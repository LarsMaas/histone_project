The text file named

<celltype>.<antibody>.densities

lists the density of <antibody> enriched fragments in <celltype>
at 25 bp resolution across each chromosome (mouse, assembly version mm8),
starting at position 0. The format is

<chr> <pos> <density>

A non-negative value indicates the number of uniquely aligned 
fragments oriented towards, and within 300 bp of, a particular 
position. Reads within 200 bp add 1 to the density, and reads 
within 300 bp add 0.25. (This is an approximation of the expected
number of fragments overlapping the postion, since the exact length
of each is unknown).

A value of -1 indicates a masked (non-unique) position.
