## Output 

The variables commonly stored in LBM simulations are the density and velocity. 
We refer to these collectively as the *fluid* variables. 

Output formats currently supported include:

* VTK legacy format (both ASCII and binary)
  - `STRUCTURED_POINTS`
  - `RECTILINEAR_GRID`
* gnuplot-compatible text files

In case of VTK, the variables are interpreted as `CELL_DATA`, meaning that the PDFs and macroscopic variables are located at the cell centers.