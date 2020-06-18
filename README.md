# ![Demo](/resources/vx_logo_v1.png?raw=true) VoronoiX (IN DEV) 

__VoronoiX__ (or VX) - __C++__ library, with __JS__ interface on 
top of it, that allows to work with different data structures and 
algorithms related to Voronoi diagram, Delaunay triangulation and 
other geometric structures. 

### Features
- [x] __QHull__ algorithm (k-dimensional convex hull)
- [x] Intersection of half-spaces (using plane-point duality: 
[blog](https://11011110.github.io/blog/2011/03/16/halfspace-intersections-and.html),
[wiki](https://en.wikipedia.org/wiki/Duality_(projective_geometry))
- [ ] Support of various spaces (even with curvature)
    - [x] 2D Plane
    - [x] 3D Space
    - [x] kD Hyperspace
    - [x] Sphere
    - [ ] Cube
- [ ] Various geometric structures and algorithms
    - [x] Voronoi diagram
    - [x] Delaunay triangulation
    - [x] Euclidean minimal spanning tree
    - [x] Minimal empty circle
    - [x] Nearest point
    - [x] Safest path
    - [x] Convex hull
- [ ] Robustness (HOPEFULLY :D )
- [x] C++ library
- [ ] JS binding
### Images and GIFs
3D Voronoi diagram inside a cube
![3D Voronoi](/resources/3d_voronoi_1.png?raw=true)

Convex hull of points randomly distributed on sphere 
![3D Sphere](/resources/3d_convex.gif?raw=true)

Convex hull of points projected from a plane onto a paraboloid. 
The edges of this convex hull corresponds to the edges of Delaunay triangulation 
on 2D plane.
![Paraboloid](/resources/3d_paraboloid.gif?raw=true)
### Build
### Use
### Contribution

 