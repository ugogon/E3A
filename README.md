# Efficient Embeddings in Exact Arithmetic
Official repository of the SIGGRAPH 2023 journal paper [[PDF](https://cybertron.cg.tu-berlin.de/projects/EEEA/media/paper.pdf)] ["Efficient Embeddings in Exact Arithmetic"](https://cybertron.cg.tu-berlin.de/projects/EEEA/).
## Code will be available soon

# Cite
```
@article{10.1145/3592445,
  author = {Finnendahl, Ugo and Bogiokas, Dimitrios and Robles Cervantes, Pablo and Alexa, Marc},
  title = {Efficient Embeddings in Exact Arithmetic},
  year = {2023},
  issue_date = {August 2023},
  publisher = {Association for Computing Machinery},
  address = {New York, NY, USA},
  volume = {42},
  number = {4},
  issn = {0730-0301},
  url = {https://doi.org/10.1145/3592445},
  doi = {10.1145/3592445},
  abstract = {We provide a set of tools for generating planar embeddings of triangulated topological spheres. The algorithms make use of Schnyder labelings and realizers. A new representation of the realizer based on dual trees leads to a simple linear time algorithm mapping from weights per triangle to barycentric coordinates and, more importantly, also in the reverse direction. The algorithms can be implemented so that all coefficients involved are 1 or -1. This enables integer computation, making all computations exact. Being a Schnyder realizer, mapping from positive triangle weights guarantees that the barycentric coordinates form an embedding. The reverse direction enables an algorithm for fixing flipped triangles in planar realizations, by mapping from coordinates to weights and adjusting the weights (without forcing them to be positive). In a range of experiments, we demonstrate that all algorithms are orders of magnitude faster than existing robust approaches.},
  journal = {ACM Trans. Graph.},
  month = {jul},
  articleno = {71},
  numpages = {17},
  keywords = {integer coordinates, schnyder labeling, parametrization}
}
```
