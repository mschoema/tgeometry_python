# Temporal Geometry Algorithms in Python

![Animation of distance computation](/img/poly_poly_dist.gif)

This repository contains the python implementation of algorithms for temporal geometries (point and polygons).  
To reproduce all the experiments, run the file [src/main.py](/src/main.py)

```
src$ python main.py
```

The file structure is as follows:
* src/
  * main.py - file used to run the code
  * viz.py - file used to create the above animation of the distance between two moving polygon
  * classes/ - contains a simple model for temporal points and fixed-shape temporal polygons
  * distance/ - contains algorithms to compute the temporal distance between temporal geometries
  * test/ - contains experiments used to test different aspects of the impemented algorithms
* img/ - contains the images output by the tests
