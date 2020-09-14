# CCNB

#### Description
A highly efficient and informative method to identify ion transport networks in fast ion conductors
By He Bing, Mi Penghui, Ye Anjiang, Chi Shuting 
we propose a highly efficient and informative method to identify interstices and connecting segments 
constructing an ion transport network by combining topological pathway network and BVSE landscape, 
which enables to obtain both the geometry and energy profiles of nonequivalent ion transport pathways between adjacent lattice sites. 
These pathways can be further used as the input for nudged elastic band calculations with automatically generated chains of images.

#### Software Architecture
1.The geometric crystal structure analysis method (CAVD model)
2.BVSE method (bvse model)
3.The combined method comprising of the geometric crystal structure analysis and BVSE methods (bvse_cavd model)
4.Automated NEB calculation model

#### Installation


  1.python 3.6
  2.ase 3.14.0
  3.pymatgen 2019.1.24
  4.scipy 1.2.1
  5.numpy 1.18.1
  6.matplotlib 3.1.3
  7.cavd 0.2.8 The package can be accessed in the repository (https://gitee.com/shuhebing/cavd)
