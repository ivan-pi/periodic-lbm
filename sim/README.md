

In principle one can mix and match. All of the following LBM methods use a regular grid.
(Some methods are also extendable to other kinds of mesh topologies.)

- [x] Standard LBM
- [ ] Thread-safe LBM
- [ ] Upwind (Lucienne Vienne)
- [x] Lax-Wendroff with FDM (Bardow's method)
- [ ] Lax-Wendroff with FEM (also known as Taylor-Galerkin)
- [x] DUGKS
- [ ] Interpolation-based
- [ ] SL-FEM
- [ ] FDM-T2S2
- [ ] Horstmann
- [ ] Shrestha
- [ ] Characteristic-based scheme (Sofonea)

There is freedom to tune the weights. The question is can this be used somehow.
Perhaps tuned dynamically, to improve stability, or reduce errors.

Spatial error, using same grid and safe CFL number.

Temporal error, using same CFL number, but varying tau. 