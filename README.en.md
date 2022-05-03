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

 The dependency package information is kept in the requirements.txt file After installing these dependency packages, you can use it by downloading the source code.For specific usage, you can refer to the examples in the example folder.

#### Installation Tutorial

   Dependency package information is kept in the requirements.txt file.

   Installation command:

```bash
pip install ccnb
````
#### Instructions

After the installation is complete, execute the   **ccnb -h**  command to query the usage.

If you use the results calculated by ccnb in your paper. please acknowledge the use of the code and cite the following papers:

> 1. He, B.; Mi, P.; Ye, A.; Chi, S.; Jiao, Y.; Zhang, L.; Pu, B.; Zou, Z.; Zhang, W.; Avdeev, M.; Adams, S.; Zhao, J.; Shi, S. A Highly Efficient and Informative Method to Identify Ion Transport Networks in Fast Ion Conductors.*Acta Materialia* **2021** ,  *203* , 116490. [https://doi.org/10.1016/j.actamat.2020.116490](https://doi.org/10.1016/j.actamat.2020.116490).
> 2. He, B.; Ye, A.; Chi, S.; Mi, P.; Ran, Y.; Zhang, L.; Zou, X.; Pu, B.; Zhao, Q.; Zou, Z.; Wang, D.; Zhang, W.; Zhao, J.; Avdeev, M.; Shi, S. CAVD, towards Better Characterization of Void Space for Ionic Transport Analysis.*Sci Data* **2020** , *7* (1), 153. [https://doi.org/10.1038/s41597-020-0491-x](https://doi.org/10.1038/s41597-020-0491-x).
> 3. He, B.; Chi, S.; Ye, A.; Mi, P.; Zhang, L.; Pu, B.; Zou, Z.; Ran, Y.; Zhao, Q.; Wang, D.; Zhang, W.; Zhao, J.; Adams, S.; Avdeev, M.; Shi, S. High-Throughput Screening Platform for Solid Electrolytes Combining Hierarchical Ion-Transport Prediction Algorithms.*Sci Data* **2020** , *7* (1), 151. [https://doi.org/10.1038/s41597-020-0474-y](https://doi.org/10.1038/s41597-020-0474-y).
