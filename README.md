# Abasys

Abasys a basic NURBS-based IsoGeometric Analysis package written in C++.


Abasys 是一个用C++编写实现的一个轻量级的等几何分析库。现阶段已经实现了对二维弹性问题的等几何算法，

[example](./example/SPSheetCircHole.cpp)中给出了一个iga经典算例的分析流程，

其中所有等几何分析相关的函数都放在**namespace abab**中，使用前请先阿巴阿巴。

采用Abasys进行等几何分析的流程可以先看[example](./example/SPSheetCircHole.cpp)中一个经典算例。

其中有关nurbs的数值算法主要依赖了[tinynurbs](https://github.com/pradeep-pyro/tinynurbs)库，并对其进行了一些修改和补充。

目前这个库还在发展和完善中，欢迎提交测试bug和修改意见，贡献出你的智慧（源码）。


