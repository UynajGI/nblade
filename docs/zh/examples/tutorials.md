# 教程示例概览

本页介绍 nblade 的教程示例，帮助你快速掌握几何代数的基本概念和操作。

## 示例列表

### 01_quickstart.py - 快速入门

**难度**: ⭐ 初级

**内容**:
- 创建几何代数
- 基向量与向量创建
- 三种基本积
- 旋转操作简介

**运行**:
```bash
python examples/tutorials/01_quickstart.py
```

**预期输出**:
```
============================================================
nblade 快速入门 - 5分钟掌握核心功能
============================================================

【1. 创建几何代数】
代数类型: 几何代数 G(3, 0, 0) (3D)
维度: 3
符号 (p, q, r): (3, 0, 0)

【2. 基向量与向量创建】
基向量 e1 = 1.000000e1
...
```

**关键概念**:
- `nblade.Algebra.euclidean(n)` - 创建 n 维欧几里得代数
- `alg.vector([...])` - 创建向量
- `*`, `^`, `|` - 几何积、外积、内积

---

### 02_basic_operations.py - 基本运算

**难度**: ⭐⭐ 基础

**内容**:
- 几何积详解
- 外积详解
- 内积详解
- 三种积的关系
- 实际应用示例

**运行**:
```bash
python examples/tutorials/02_basic_operations.py
```

**关键概念**:
- 核心公式: `ab = a·b + a∧b`
- 外积表示张成的子空间
- 内积表示投影关系

---

### 03_vectors.py - 向量操作

**难度**: ⭐⭐ 基础

**内容**:
- 向量的多种创建方法
- 向量属性（范数、系数等）
- 代数运算（加、减、乘）
- 几何运算（投影、反射）
- 对合运算
- 高维向量

**运行**:
```bash
python examples/tutorials/03_vectors.py
```

**关键概念**:
- `alg.basis_vectors()` - 获取所有基向量
- `v.norm()` - 计算范数
- `v.project_to(blade)` - 投影
- `v.reflect_in(n)` - 反射

---

### 04_dual_operations.py - 对偶运算

**难度**: ⭐⭐ 基础

**内容**:
- 对偶运算的基本概念
- Hodge 对偶的计算方法
- 对偶与叉积的关系
- 逆对偶运算
- 不同维度下的对偶特性

**运行**:
```bash
python examples/tutorials/04_dual_operations.py
```

**关键概念**:
- 对偶将 r 阶元素映射到 (n-r) 阶元素
- 3D 中: `a × b = (a ∧ b)*`
- `*` 表示对偶运算
- 向量的对偶表示其正交补空间

---

### 05_rotations.py - 旋转教程

**难度**: ⭐⭐⭐ 进阶

**内容**:
- 2D 旋转与复数
- 3D 旋转与转子
- 组合多个旋转
- 旋转插值（SLERP）
- 转子与四元数的关系
- 实际应用示例

**运行**:
```bash
python examples/tutorials/05_rotations.py
```

**关键概念**:
- 转子: `R = exp(-B·θ/2)`
- 旋转公式: `v' = R·v·R†`
- 组合旋转: `R_combined = R2·R1`

---

### 06_reciprocal_frame.py - 互逆标架

**难度**: ⭐⭐⭐ 进阶

**内容**:
- 什么是互逆标架
- 如何计算互逆标架
- 验证互逆条件
- 度量张量的计算
- 互逆标架的应用

**运行**:
```bash
python examples/tutorials/06_reciprocal_frame.py
```

**关键概念**:
- 互逆标架满足 `aⁱ⌋aⱼ = δⁱⱼ`
- 正交基的互逆标架就是其本身
- 非正交基的互逆标架需要计算
- 互逆标架用于在任意标架下分解向量
- `reciprocal_frame()` 函数计算互逆标架
- `metric_tensor()` 函数计算度量张量

---

### 07_grade_operations.py - 阶次操作

**难度**: ⭐⭐ 基础

**内容**:
- 阶次的概念（0阶到n阶）
- 阶次提取
- 偶部和奇部分解
- 阶次对合（Grade Involution）
- 反序（Reversion）
- Clifford 共轭

**运行**:
```bash
python examples/tutorials/07_grade_operations.py
```

**关键概念**:
- 阶次表示基向量的数量
- `mv.grade(r)` - 提取 r 阶部分
- `mv.even_part()` / `mv.odd_part()` - 偶部/奇部
- `~mv` - 阶次对合
- `mv.reversion()` - 反序
- `mv.clifford_conjugate()` - Clifford 共轭

---

## 学习路径

```
01_quickstart.py
       ↓
02_basic_operations.py
       ↓
03_vectors.py
       ↓
04_dual_operations.py
       ↓
07_grade_operations.py
       ↓
05_rotations.py
       ↓
06_reciprocal_frame.py
       ↓
   应用示例
```

## 常见问题

### Q: 为什么外积结果是二重向量？

A: 外积 `a ^ b` 表示两个向量张成的子空间（平面）。在几何代数中，平面由二重向量表示，而不是像传统向量分析那样由法向量表示。

### Q: 几何积与内/外积的关系是什么？

A: 几何积是基本运算，分解为内积和外积: `ab = a·b + a∧b`。内积是对称部分（标量），外积是反对称部分（二重向量）。

### Q: 转子与四元数的关系是什么？

A: 四元数是几何代数的偶子代数。四元数单位 i, j, k 对应于二重向量 e₂∧e₃, e₃∧e₁, e₁∧e₂。

### Q: 什么是对偶运算？

A: 对偶运算将 r 阶元素映射到 (n-r) 阶元素，其中 n 是空间维度。在 3D 中，向量的对偶是二重向量，表示垂直于原向量的平面。

### Q: 什么是互逆标架？

A: 互逆标架是一组向量，与原标架满足 `aⁱ⌋aⱼ = δⁱⱼ`（Kronecker delta）。它用于在任意标架下分解向量，特别是在非正交坐标系中。

---

## 下一步

完成教程示例后，你可以继续学习:

- [物理示例](physics.md) - 刚体动力学
- [图形学示例](cg.md) - 几何变换

或查看 [API 参考](../api/algebra.md) 获取更多详细信息。
