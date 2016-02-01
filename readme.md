##目的
用来熟悉一些数学算法，所有的代码都没有进行优化，以便于更好的理解算法本身。
##实现一些代数计算
##编译方法
windows下需要cygwin
```
autoreconf --install
./configure
make
```
###LU分解
考虑到数值精确性，每次将列中绝对值最大的作为主元。
###crout LU方法
crout方法比标准的LU分解要快一倍左右。
###Cholesky分解
[算法](http://mathfaculty.fullerton.edu/mathews/n2003/CholeskyMod.html)
[算法](http://www.netlib.org/utk/papers/factor/node9.html)
###解逆矩阵
###解线性方程Ax=b
##微分方程的数值解
###欧拉法数值解
###中点法数值解