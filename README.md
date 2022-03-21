# periodic-lbm

A collection of research codes for LBM in periodic domains

## Data layout 

When organizing an LBM code, there are two main choices for the data layout
- structure of arrays (SoA), and
- array of structures (AoS)
The chosen layout has important performance implications. For optimum performance the kernels should be specialized for the chosen layout. 

In the SoA layout, the PDF's in a given direction at different spatial points are stored contiguosly. In the AoS layout the PDF's pointing in different directions at a single spatial point are stored contiguously.



## See also

- [Medusa](https://e6.ijs.si/medusa/wiki/index.php/Medusa)
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [Boost](https://www.boost.org/)
- [nanoflann](https://github.com/jlblancoc/nanoflann)
- [psblas3](https://github.com/sfilippone/psblas3)