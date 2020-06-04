# CCNB

#### 介绍
A highly efficient and informative method to identify ion transport networks in fast ion conductors
By He Bing, Mi Penghui, Ye Anjiang, Chi Shuting 
我们提出了一种更高效的方法通过结合拓扑路径网络和BVSE能量势场构建来离子输运网络，从而得到相邻晶格点间非等效离子输运路径的几何和能量属性。
这些路径信息可进一步用作输入，用于自动生成images的NEB计算。

#### 软件架构

(1) CAVD计算几何拓扑路径网络
(2) bvse模块
(1) 融合几何分析方法和键价和方法等经验方法来计算和分析晶体结构的离子输运通道网络，并得到迁移离子晶格位之间的输运路径的几何和能量属性。
(2) 为了实现并加速自动NEB计算，将沿着融合几何分析和键价和方法计算出的路径自动生成过渡态，并自动生成NEB计算需要的配置文件。


#### 安装教程
  所需依赖包
  1.python 3.6
  2.ase 3.14.0
  3.pymatgen 2019.1.24
  4.scipy 1.2.1
  5.numpy 1.18.1
  6.matplotlib 3.1.3
  7.cavd 0.2.8 The package can be accessed in the repository (https://gitee.com/shuhebing/cavd)

#### 使用说明

1.  xxxx
2.  xxxx
3.  xxxx

#### 参与贡献

1.  Fork 本仓库
2.  新建 Feat_xxx 分支
3.  提交代码
4.  新建 Pull Request


#### 码云特技

1.  使用 Readme\_XXX.md 来支持不同的语言，例如 Readme\_en.md, Readme\_zh.md
2.  码云官方博客 [blog.gitee.com](https://blog.gitee.com)
3.  你可以 [https://gitee.com/explore](https://gitee.com/explore) 这个地址来了解码云上的优秀开源项目
4.  [GVP](https://gitee.com/gvp) 全称是码云最有价值开源项目，是码云综合评定出的优秀开源项目
5.  码云官方提供的使用手册 [https://gitee.com/help](https://gitee.com/help)
6.  码云封面人物是一档用来展示码云会员风采的栏目 [https://gitee.com/gitee-stars/](https://gitee.com/gitee-stars/)
