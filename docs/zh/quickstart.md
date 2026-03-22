# nblade 快速入门

## 安装

### 使用 pip 安装

```bash
pip install nblade
```

### 从源码安装

```bash
# 安装 maturin
pip install maturin

# 克隆仓库
git clone https://github.com/UynajGI/nblade.git
cd nblade

# 构建并安装
maturin develop --release
```

## 基本概念

### 什么是几何代数？

几何代数（Geometric Algebra，又称 Clifford Algebra）是传统向量代数的推广。它提供了：

1. **统一的运算系统**：几何积统一了内积和外积
2. **几何对象表示**：标量、向量、平面、体积等都用同一种对象表示
3. **坐标系无关**：运算结果不依赖于特定坐标系

### 核心公式

几何积是最基本的运算：

```
ab = a·b + a∧b
```

其中：
- `ab` 是几何积
- `a·b` 是内积（对称部分）
- `a∧b` 是外积（反对称部分）

## 快速开始

### 创建几何代数

```python
import nblade

# 创建 3D 欧几里得几何代数 G(3,0,0)
alg = nblade.Algebra.euclidean(3)

# 查看代数信息
print(alg)  # Geometric Algebra G(3, 0, 0) (3D)
```

### 创建向量

```python
# 从列表创建向量
v = alg.vector([1.0, 2.0, 3.0])

# 获取基向量
e1, e2, e3 = alg.basis_vectors()

# 创建标量
s = alg.scalar(5.0)
```

### 基本运算

```python
# 创建两个向量
a = alg.vector([1.0, 0.0, 0.0])
b = alg.vector([1.0, 1.0, 0.0])

# 几何积
geom = a * b

# 外积（楔积）
outer = a ^ b

# 内积
inner = a | b
```

### 旋转操作

```python
import math

# 创建旋转平面
plane = e1 ^ e2  # xy 平面

# 创建转子（旋转 45 度）
rotor = alg.rotor(plane, math.pi / 4)

# 旋转向量
rotated = e1.rotate_by(rotor)
```

## 下一步

- [基础教程](tutorial/01_basics.md) - 深入了解几何代数基础
- [API 参考](api/algebra.md) - 查看完整 API 文档
- [示例代码](https://github.com/UynajGI/nblade/tree/main/examples) - 更多实用示例