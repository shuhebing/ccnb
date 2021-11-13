# CCBN

#### Description
A highly efficient and informative method to identify ion transport networks in fast ion conductors
By He Bing, Mi Penghui, Ye Anjiang, Chi Shuting 
we propose a highly efficient and informative method to identify interstices and connecting segments 
constructing an ion transport network by combining topological pathway network and BVSE landscape, 
which enables to obtain both the geometry and energy profiles of nonequivalent ion transport pathways between adjacent lattice sites. 
These pathways can be further used as the input for nudged elastic band calculations with automatically generated chains of images.

#### Software Architecture



(1)	BVAnalysis

This module implements three models: BVS, BVSE and BVEL. 
Through this module, the BVSE potential field required for the combined method comprising of the geometric crystal structure analysis and BVSE methods can be calculated.

(2)	cavd_channel

This module is a cavd interface module，CCBN provides this interface to shield the implementation of the underlying cavd,
which is convenient for relevant researchers to get the interstitial network required for the combined method comprising of the geometric crystal structure analysis and BVSE methods.

(3)	MergeCluster

This module implements the algorithm of merging interstice clusters.

(4)	MigrationPath 

This module implements the algorithm for finding the MEP based on the BVSE landscape.

(5)	MigrationNetwork

This module implements the algorithm for calculateing all nonequivalent transport pathways between adjacent lattice sites.

(6)	neb_packages

This module implements a FP-NEB automation calculation algorithm by utilizing these pathways information found by our combined method. 

#### Requirements


  1.python 3.6
  
  2.ase 3.14.0
  
  3.pymatgen 2019.1.24
  
  4.scipy 1.2.1
  
  5.numpy 1.18.1
  
  6.matplotlib 3.1.3
  
  7.cavd 0.2.8 The package can be accessed in the repository (https://gitee.com/shuhebing/cavd)
  
  After installing these dependency packages, you can use it by downloading the source code.For specific usage, you can refer to the examples in the example folder.
