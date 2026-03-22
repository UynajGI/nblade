# nblade - N-dimensional Blade

**高性能几何代数 Python 库** | **High-Performance Geometric Algebra Library**

[![PyPI](https://img.shields.io/pypi/v/nblade.svg)](https://pypi.org/project/nblade/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## 概述 | Overview

nblade 是一个基于 Rust 实现的高性能几何代数库，支持任意维度（最多 64 维）和任意度量签名 G(p, q, r)。

nblade is a high-performance geometric algebra library implemented in Rust, supporting arbitrary dimensions (up to 64D) and arbitrary metric signatures G(p, q, r).

## 特性 | Features

- **任意维度 | Arbitrary Dimensions**: 支持最多 64 维向量空间 | Support for up to 64-dimensional vector spaces
- **任意签名 | Arbitrary Signatures**: 支持 G(p, q, r) 任意度量签名 | Support for G(p, q, r) metric signatures
- **高性能 | High Performance**: Rust 实现，支持并行计算 | Rust backend with parallel computing
- **双表示 | Dual Representation**: 自动选择密集或稀疏表示 | Automatic dense/sparse selection
- **NumPy 集成 | NumPy Integration**: 直接从 NumPy 数组创建向量 | Create vectors directly from NumPy arrays

## 快速开始 | Quick Start

### 安装 | Installation

```bash
pip install nblade
```

### 示例 | Example

```python
import nblade
import math

# 创建 3D 欧几里得几何代数 | Create 3D Euclidean geometric algebra
alg = nblade.Algebra.euclidean(3)

# 获取基向量 | Get basis vectors
e1, e2, e3 = alg.basis_vectors()

# 创建向量 | Create a vector
v = alg.vector([1.0, 2.0, 3.0])

# 几何积 | Geometric product
product = e1 * e2  # Results in bivector e12

# 外积 | Outer product
wedge = e1 ^ e2  # Results in bivector e12

# 内积 | Inner product
inner = e1 | e2  # Results in 0 (orthogonal vectors)

# 旋转 | Rotation
plane = e1 ^ e2
rotor = alg.rotor(plane, math.pi / 2)  # 90-degree rotation
rotated = e1.rotate_by(rotor)
```

## 文档 | Documentation

- [快速入门 | Quick Start](zh/quickstart.md)
- [基础教程 | Basics Tutorial](zh/tutorial/01_basics.md)
- [API 参考 | API Reference](api/algebra.md)

## 示例代码 | Examples

| 示例 | 说明 | Description |
|------|------|-------------|
| [quickstart.py](https://github.com/nblade/nblade/blob/main/examples/tutorials/01_quickstart.py) | 5 分钟入门 | 5-minute introduction |
| [basic_operations.py](https://github.com/nblade/nblade/blob/main/examples/tutorials/02_basic_operations.py) | 三种基本乘积 | Three fundamental products |
| [rotations.py](https://github.com/nblade/nblade/blob/main/examples/tutorials/05_rotations.py) | 旋转操作 | Rotation operations |
| [rigid_body.py](https://github.com/nblade/nblade/blob/main/examples/physics/01_rigid_body.py) | 刚体物理 | Rigid body physics |
| [transformations.py](https://github.com/nblade/nblade/blob/main/examples/cg/01_transformations.py) | 计算机图形学变换 | CG transformations |

## 支持的运算 | Supported Operations

| 运算 | 符号 | 说明 |
|------|------|------|
| 几何积 | `a * b` | 几何代数的基本运算 |
| 外积 | `a ^ b` | 表示张成的子空间 |
| 内积 | `a \| b` | 表示投影关系 |
| 对偶 | `a.dual()` | Hodge 对偶 |
| 旋转 | `a.rotate_by(R)` | 使用转子旋转 |
| 投影 | `a.project_to(b)` | 投影到子空间 |
| 反射 | `a.reflect_in(n)` | 关于平面反射 |

## 数学背景 | Mathematical Background

几何代数定义了统一的多重向量运算：

```
ab = a·b + a∧b
```

- **几何积 (Geometric Product)**: `ab` - 包含完整信息的乘法
- **内积 (Inner Product)**: `a·b` - 对称部分，表示投影
- **外积 (Outer Product)**: `a∧b` - 反对称部分，表示张成的子空间

## 许可证 | License

MIT License