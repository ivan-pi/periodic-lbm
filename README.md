# periodic-lbm

A personal collection of research codes for lattice Boltzmann methods in periodic 2D domains.

For actual simulation work take a look at [other codes](/other_codes).

## Data layout 

When organizing an LBM code, there are two main choices for the data layout
- structure of arrays (SoA), and
- array of structures (AoS)
The chosen layout has important performance implications. For optimum performance the collision and streaming kernels should be specialized for the chosen layout. 

In the SoA layout, the PDF's pointing in a given direction at different spatial points are stored contiguously. 
In the AoS layout, the PDF's pointing in different directions at a given spatial point are stored contiguously.

## See also

- [Medusa](https://e6.ijs.si/medusa/wiki/index.php/Medusa)
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [Boost](https://www.boost.org/)
- [nanoflann](https://github.com/jlblancoc/nanoflann)
- [psblas3](https://github.com/sfilippone/psblas3)
- [librsb](http://librsb.sourceforge.net/)

## Other codes

For large-scale simulations consider using one of the following freely availables codes
* [waLBerla](https://walberla.net/)
* [Palabos](https://palabos.unige.ch/)
* [TCLB](https://github.com/CFD-GO/TCLB)
* [Musubi](https://geb.inf.tu-dresden.de/doxy/musubi/index.html)
* [OpenLB](https://www.openlb.net/)
* [VirtualFluids](https://git.rz.tu-bs.de/irmb/virtualfluids)

For industrial purpose simulations the following solvers can also be used
* Dassault Systèmes Simulia [PowerFLOW](https://www.3ds.com/products-services/simulia/products/powerflow/) and [XFlow](https://www.3ds.com/products-services/simulia/products/xflow/) (formerly from Exa Corporation)
* Altair [ultraFluidX](https://www.altair.com/altair-cfd-capabilities/#lbm) (formerly from FluiDyna GmbH)
* [ProLB](http://www.prolb-cfd.com/)

A bunch of other (mostly research oriented) codes are described in [the list](https://github.com/sthavishtha/list-lattice-Boltzmann-codes) by @sthavishtha.
