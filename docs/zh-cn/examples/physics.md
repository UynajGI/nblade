# 物理示例概览

本页介绍 nblade 在物理学中的应用示例。

## 示例列表

### 01_rigid_body.py - 刚体动力学

**难度**: ⭐⭐⭐ 高级

**内容**:
- 角动量作为二重向量
- 惯性张量
- 刚体旋转动力学
- 几何代数形式的欧拉方程
- 转动动能
- 简单刚体模拟

**运行**:
```bash
python examples/physics/01_rigid_body.py
```

**预期输出**:
```
============================================================
nblade Rigid Body Physics Example
Application of Geometric Algebra in Classical Mechanics
============================================================

============================================================
Section 1: Angular Momentum as Bivector
============================================================

In traditional vector analysis, angular momentum L = r × p
In geometric algebra, angular momentum naturally represents as bivector L = r ∧ p

Position vector r = 3.000000e1 + 4.000000e2
Momentum vector p = 2.000000e1 + -1.000000e2 + 5.000000e3

Angular momentum bivector L = r∧p = -11.000000e1∧e2 + 15.000000e1∧e3 + 20.000000e2∧e3
...
```

---

### 02_kepler_problem.py - 开普勒问题

**难度**: ⭐⭐⭐ 高级

**内容**:
- 角动量作为二重向量 L = r ∧ p
- 中心力场运动
- 轨道方程的几何代数推导
- 偏心率矢量（拉普拉斯-龙格-楞次矢量）
- 数值模拟与守恒量验证

**运行**:
```bash
python examples/physics/02_kepler_problem.py
```

**预期输出**:
```
============================================================
nblade 开普勒问题示例
几何代数在轨道力学中的应用
============================================================

============================================================
第一节：角动量作为二重向量
============================================================

在中心力场中，角动量守恒是核心性质
传统方法：L = r × p (叉积，仅限3D)
几何代数：L = r ∧ p (外积，任意维度)

位置向量 r = 5.000000e1
动量向量 p = m·v = 1.000000e2

角动量 L = r ∧ p = 5.000000e1∧e2
...
```

---

## 核心概念

### 角动量的几何表示

在传统向量分析中，角动量定义为叉积：

$$\mathbf{L} = \mathbf{r} \times \mathbf{p}$$

在几何代数中，角动量自然地表示为二重向量：

$$L = r \wedge p$$

这种表示的优势：
- 适用于任意维度（叉积仅在3D有效）
- 几何意义清晰（张成的平面）
- 与几何代数的其他运算保持一致

### 与传统叉积的关系

在三维空间中，角动量的对偶对应于传统叉积：

```python
L_bivector = r ^ p        # 几何代数表示
L_vector = L_bivector.dual()  # 对应于传统叉积结果
```

### 刚体旋转

刚体姿态用转子描述：

```python
# 角速度平面
omega_plane = omega.dual()

# 更新转子
delta_rotor = alg.rotor(omega_plane, omega_magnitude * dt)
rotor = delta_rotor * rotor
```

### 开普勒问题中的守恒量

在中心力场 F = -k/r² · r̂ 中，有两个重要的守恒量：

**角动量**:
$$L = r \wedge p$$

角动量二重向量表示位置和动量张成的有向面积，其方向垂直于轨道平面。

**偏心率矢量**（拉普拉斯-龙格-楞次矢量）:
$$\mathbf{e} = \frac{\mathbf{L} \times \mathbf{v}}{k} - \hat{\mathbf{r}}$$

偏心率矢量指向轨道的近日点方向，其大小等于轨道偏心率。

---

## 理论背景

### 欧拉方程

无外力矩刚体的欧拉方程：

$$\frac{d\mathbf{L}}{dt} + \boldsymbol{\omega} \times \mathbf{L} = 0$$

在几何代数中：

$$\frac{dL}{dt} + \omega \wedge L = 0$$

### 转动动能

$$T = \frac{1}{2} \omega \cdot L = \frac{1}{2} \sum_i I_i \omega_i^2$$

在 nblade 中：

```python
T = 0.5 * (omega | L).scalar_part()
```

### 轨道方程

开普勒问题的极坐标轨道方程：

$$r(\theta) = \frac{l}{1 + e \cos \theta}$$

其中 l 是半正焦弦，e 是偏心率。

半正焦弦与角动量的关系：

$$l = \frac{L^2}{mk}$$

---

## 常见问题

### Q: 为什么用二重向量表示角动量？

A: 二重向量明确表示旋转平面，比传统向量（表示旋转轴）更加几何直观。此外，二重向量在任何维度都有效。

### Q: 如何从角速度计算转子？

A: 角速度的对偶就是旋转平面，然后用该平面创建转子：

```python
omega = alg.vector([0, 0, 1])  # 绕 z 轴
omega_plane = omega.dual()      # e1∧e2
rotor = alg.rotor(omega_plane, angle)
```

### Q: 几何代数如何处理开普勒问题？

A: 几何代数通过以下方式简化开普勒问题：
- 角动量 L = r ∧ p 作为二重向量自然表示
- 偏心率矢量作为守恒量统一处理
- 所有运算在任意维度保持一致
- 避免了复杂的坐标系变换

---

## 参考文献

- [Geometric Algebra for Physicists](https://www.cambridge.org/core/books/geometric-algebra-for-physicists/) - Doran & Lasenby
- [New Foundations for Classical Mechanics](https://www.springer.com/gp/book/9780792355141) - David Hestenes

---

## 下一步

- [图形学示例](cg.md) - 几何变换
- [教程示例](tutorials.md) - 基础概念
