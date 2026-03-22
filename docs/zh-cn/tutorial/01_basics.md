# 几何代数基础教程

本教程介绍几何代数的基本概念和 nblade 的使用方法。

## 目录

1. [几何代数简介](#1-几何代数简介)
2. [多重向量](#2-多重向量)
3. [几何积](#3-几何积)
4. [外积](#4-外积)
5. [内积](#5-内积)
6. [对合运算](#6-对合运算)
7. [对偶](#7-对偶)

---

## 1. 几何代数简介

### 什么是几何代数？

几何代数（Geometric Algebra）是由 Hermann Grassmann 和 William Kingdon Clifford 在 19 世纪发展起来的数学体系。它将：

- **标量**（0 维）
- **向量**（1 维）
- **二重向量**（2 维，表示平面）
- **三重向量**（3 维，表示体积）
- ... 等等

统一到 **多重向量（Multivector）** 的概念中。

### 为什么使用几何代数？

| 传统方法 | 几何代数方法 |
|---------|-------------|
| 使用矩阵表示旋转 | 使用转子，更简洁 |
| 叉积只在 3D 有效 | 外积在任意维度有效 |
| 复数、四元数是独立概念 | 都是几何代数的子代数 |
| 需要多种运算规则 | 几何积统一一切 |

## 2. 多重向量

多重向量是几何代数中的基本对象，是不同阶次元素的线性组合。

### 阶次（Grade）

- **0 阶**：标量，如 `5`
- **1 阶**：向量，如 `e1`, `2e1 + 3e2`
- **2 阶**：二重向量，如 `e1∧e2`（表示平面）
- **3 阶**：三重向量，如 `e1∧e2∧e3`（表示体积）

### 在 nblade 中创建多重向量

```python
import nblade

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

# 标量
s = alg.scalar(5.0)

# 向量
v = alg.vector([1.0, 2.0, 3.0])

# 二重向量
B = e1 ^ e2

# 混合多重向量
mixed = s + v + B
```

### 提取特定阶次

```python
# 获取标量部分
scalar_part = mixed.grade(0)

# 获取向量部分
vector_part = mixed.grade(1)

# 获取二重向量部分
bivector_part = mixed.grade(2)
```

## 3. 几何积

几何积是几何代数最核心的运算。

### 定义

对于两个向量 `a` 和 `b`：

```
ab = a·b + a∧b
```

- `a·b`：内积，是对称部分（标量）
- `a∧b`：外积，是反对称部分（二重向量）

### 在 nblade 中使用

```python
a = alg.vector([1.0, 2.0, 0.0])
b = alg.vector([3.0, 1.0, 0.0])

# 几何积
product = a * b
print(product)  # 标量 + 二重向量
```

### 性质

1. **结合律**：`(ab)c = a(bc)`
2. **分配律**：`a(b + c) = ab + ac`
3. **非交换**：`ab ≠ ba`（一般情况）

### 基向量的几何积

```python
# 正交向量
print(e1 * e2)  # e12（二重向量）

# 平行向量
print(e1 * e1)  # 1（标量）
```

## 4. 外积

外积（楔积）表示两个向量张成的子空间。

### 几何意义

- 向量 `a` 和 `b` 的外积 `a∧b` 是一个 **二重向量**
- 表示 `a` 和 `b` 张成的平面
- 其范数等于平行四边形的面积

### 反对称性

```python
print(e1 ^ e2)   # e12
print(e2 ^ e1)   # -e12
print(e1 ^ e1)   # 0
```

### 高阶外积

```python
# 三重向量（体积）
volume = e1 ^ e2 ^ e3
print(volume)  # e123

# 在 3D 中，这是伪标量
```

### 计算面积

```python
a = alg.vector([2.0, 0.0, 0.0])
b = alg.vector([1.0, 1.0, 0.0])

area = (a ^ b).norm()
print(f"平行四边形面积: {area}")  # 2.0
```

## 5. 内积

内积表示投影关系。

### 向量的内积

```python
a = alg.vector([1.0, 2.0, 3.0])
b = alg.vector([4.0, 5.0, 6.0])

# 内积（结果是标量）
inner = a | b
print(inner)  # 32
```

### 左收缩和右收缩

```python
# 左内积（左收缩）
left = a | b

# 右内积（右收缩）
right = a.right_inner(b)
```

### 几何意义

内积可以用来：

1. **计算投影**：`(a|b)/|b|² × b` 是 `a` 在 `b` 上的投影
2. **判断正交**：如果 `a|b = 0`，则 `a` 和 `b` 正交
3. **计算夹角**：`cos(θ) = (a|b)/(|a||b|)`

## 6. 对合运算

对合是将多重向量映射到自身的运算。

### 阶次对合（Grade Involution）

对于 r 阶次元素，乘以 `(-1)^r`：

```python
v = alg.vector([1, 2, 3])
print(v.grade_involution())  # -v（向量是 1 阶，(-1)^1 = -1）
```

### 反转（Reversion）

对于 r 阶次元素，乘以 `(-1)^(r(r-1)/2)`：

```python
print(v.reversion())  # v（向量不变）
```

### Clifford 共轭

Clifford 共轭 = 阶次对合 + 反转：

```python
print(v.clifford_conjugate())
```

### 应用

对合运算主要用于计算：

- **范数**：`|A|² = ⟨A†A⟩₀`
- **逆**：`A⁻¹ = A†/|A|²`（对于可逆元素）

## 7. 对偶

对偶是将 k-向量映射到 (n-k)-向量的运算。

### 定义

在 n 维空间中，向量 `v` 的对偶是：

```
v* = v·I⁻¹
```

其中 `I` 是伪标量（最高阶基向量）。

### 在 nblade 中使用

```python
v = alg.vector([1, 2, 3])
dual_v = v.dual()
print(dual_v)  # 二重向量

# 逆对偶
original = dual_v.inverse_dual()
```

### 几何意义

在 3D 中：

- 向量的对偶是二重向量（平面）
- 二重向量的对偶是向量
- 这建立了 **叉积** 和 **外积** 的联系：`a × b = (a ∧ b)*`

### 伪标量

```python
# 获取伪标量
I = alg.config.volume_element()
print(I)  # e123

# 伪标量的平方
print(I * I)  # -1（在 3D 欧几里得空间）
```

---

## 下一步

- [旋转教程](../tutorials/05_rotations.py) - 学习如何使用转子进行旋转
- [API 参考](../api/algebra.md) - 查看完整 API 文档
- [示例代码](https://github.com/nblade/nblade/tree/main/examples) - 更多实用示例