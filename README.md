# Abasys

Abasys a basic NURBS-based IsoGeometric Analysis package written in C++.


Abasys 是一个用C++编写实现的轻量级的等几何分析库。现阶段可以对二维线弹性问题进行求解，具体可以看

[example](./example/SPSheetCircHole.cpp)中给出的一个iga经典算例的分析程序。

Abasys主要由三部分组成，首先最上层便是进行等几何分析的iga子模块，其所有程序放在命名空间abab（阿巴阿巴）下。

其中有关nurbs的数值算法主要参考了[tinynurbs](https://github.com/pradeep-pyro/tinynurbs)库，并对其进行了调整和补充。

在nurbs下的detail文件夹内，是可查阅到的nurbs经典算法的函数实现，我们将其都放在命名空间funs下。

这些函数最为稳定且高效，是上层封装的基础，最重要的是这一部分对学习和理解nurbs提供了很好的帮助。

目前这个库尚有许多不足，欢迎大家提交bug和修改意见，也可以提供算例程序，一起发展和完善Abasys。


