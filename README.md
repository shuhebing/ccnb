# CCBN

#### 介绍

使用Python3编程语言开发了一个Python库CCBN来实现了几何分析和BVSE的融合方法的所有功能

#### 软件架构

(1)	BVAnalysis

在该模块中实现了键价和方法，提供了BVS、BVSE、BVEL三种模型的计算功能。通过该模块可以计算出几何分析和BVSE的融合方法计算所需要的BVSE势场文件。

(2)	cavd_channel

该模块是一个CAVD接口模块，调用了CAVD计算间隙网络功能。CCBN提供了该接口用以屏蔽底层CAVD的实现，方便相关研究人员产生几何分析和BVSE融合方法计算所需要的间隙网络输入文件。

(3)	MergeCluster

该模块实现了3.2.1节所叙述的合并间隙簇算法。

(4)	MigrationPath

该模块实现了3.1节所叙述的基于BVSE势场的寻找最小能量路径算法。

(5)	MigrationNetwork

该模块实现了3.2.2节所叙述的计算离子输运网络和3.3节叙述的相邻晶格位之间的非等价路径计算算法。

(6)	neb_packages

该模块实现了3.4节叙述的自动化DFT-NEB计算算法，该模块可根据相邻晶格位之间的非等价路径计算算法的计算结果自动产生第一性原理DFT-NEB方法计算所需的POSCAR文件，以及自动生成其他的例如POTCAR、INCAR、KPOINTS等文件。

#### 安装教程

  依赖包信息保放在requirements.txt文件中.

  安装命令：

```bash
pip intall ccnb
```

#### 使用方法

安装完成以后，执行ccnb -h 命令查询使用方法。

如果您在论文中使用 ccnb 计算的结果。 请确认代码的使用并引用以下论文：

> 1. He, B.; Mi, P.; Ye, A.; Chi, S.; Jiao, Y.; Zhang, L.; Pu, B.; Zou, Z.; Zhang, W.; Avdeev, M.; Adams, S.; Zhao, J.; Shi, S. A Highly Efficient and Informative Method to Identify Ion Transport Networks in Fast Ion Conductors.*Acta Materialia* **2021** ,  *203* , 116490. [https://doi.org/10.1016/j.actamat.2020.116490](https://doi.org/10.1016/j.actamat.2020.116490).
> 2. He, B.; Ye, A.; Chi, S.; Mi, P.; Ran, Y.; Zhang, L.; Zou, X.; Pu, B.; Zhao, Q.; Zou, Z.; Wang, D.; Zhang, W.; Zhao, J.; Avdeev, M.; Shi, S. CAVD, towards Better Characterization of Void Space for Ionic Transport Analysis.*Sci Data* **2020** , *7* (1), 153. [https://doi.org/10.1038/s41597-020-0491-x](https://doi.org/10.1038/s41597-020-0491-x).
> 3. He, B.; Chi, S.; Ye, A.; Mi, P.; Zhang, L.; Pu, B.; Zou, Z.; Ran, Y.; Zhao, Q.; Wang, D.; Zhang, W.; Zhao, J.; Adams, S.; Avdeev, M.; Shi, S. High-Throughput Screening Platform for Solid Electrolytes Combining Hierarchical Ion-Transport Prediction Algorithms.*Sci Data* **2020** , *7* (1), 151. [https://doi.org/10.1038/s41597-020-0474-y](https://doi.org/10.1038/s41597-020-0474-y).
