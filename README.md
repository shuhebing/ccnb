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
  所需依赖包
  1.python 3.6
  2.ase 3.14.0
  3.pymatgen 2019.1.24
  4.scipy 1.2.1
  5.numpy 1.18.1
  6.matplotlib 3.1.3
  7.cavd 0.2.8 The package can be accessed in the repository (https://gitee.com/shuhebing/cavd)

