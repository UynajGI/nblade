# Algebra 类 API 参考

`Algebra` 类是 nblade 的核心入口，用于创建和配置几何代数空间。

## 类定义

```python
class Algebra:
    def __init__(self, dimension: int, p: int = 0, q: int = 0, r: int = 0) -> None
```

### 参数

| 参数 | 类型 | 说明 |
|------|------|------|
| `dimension` | `int` | 向量空间维度 (1-64) |
| `p` | `int` | 正平方基向量数量 (e_i² = +1) |
| `q` | `int` | 负平方基向量数量 (e_i² = -1) |
| `r` | `int` | 零平方基向量数量 (e_i² = 0) |

### 示例

```python
import nblade

# 创建自定义签名代数 G(2,1,0)
alg = nblade.Algebra(dimension=3, p=2, q=1, r=0)
```

---

## 工厂方法

### `euclidean(dimension)`

创建欧几里得几何代数 G(n, 0, 0)。

```python
@classmethod
def euclidean(cls, dimension: int) -> Algebra
```

**参数**: `dimension` - 空间维度

**返回**: 配置为欧几里得签名的 Algebra 实例

**示例**:
```python
# 3D 欧几里得几何代数
alg = nblade.Algebra.euclidean(3)

# 所有基向量满足 e_i² = +1
```

---

### `spacetime(dimension)`

创建时空代数 G(1, n-1, 0)。

```python
@classmethod
def spacetime(cls, dimension: int) -> Algebra
```

**参数**: `dimension` - 时空维度 (通常为 4)

**返回**: 配置为时空签名的 Algebra 实例

**示例**:
```python
# 相对论时空代数 G(1,3,0)
sta = nblade.Algebra.spacetime(4)

# e0² = +1 (时间方向)
# e1² = e2² = e3² = -1 (空间方向)
```

---

### `cga()`

创建共形几何代数 G(4, 1, 0)。

```python
@classmethod
def cga(cls) -> Algebra
```

**返回**: 配置为共形几何代数的 Algebra 实例

**示例**:
```python
# 共形几何代数，用于表示点、线、圆、球等几何对象
cga = nblade.Algebra.cga()
```

---

## 属性

### `dimension`

```python
@property
def dimension(self) -> int
```

返回向量空间的维度。

**示例**:
```python
alg = nblade.Algebra.euclidean(3)
print(alg.dimension)  # 3
```

---

### `basis_count`

```python
@property
def basis_count(self) -> int
```

返回基的数量 (2^dimension)。

**示例**:
```python
alg = nblade.Algebra.euclidean(3)
print(alg.basis_count)  # 8 (标量, 3个向量, 3个二重向量, 1个三重向量)
```

---

### `signature`

```python
@property
def signature(self) -> Tuple[int, int, int]
```

返回度量签名 (p, q, r)。

**示例**:
```python
alg = nblade.Algebra.spacetime(4)
print(alg.signature)  # (1, 3, 0)
```

---

### `config`

```python
@property
def config(self) -> AlgebraConfig
```

返回底层 AlgebraConfig 对象，用于底层 API 调用。

---

## 向量创建方法

### `vector(data)`

从列表或元组创建向量 (1-向量)。

```python
def vector(self, data: Union[List[float], Tuple[float, ...]]) -> MultiVector
```

**参数**: `data` - 长度为 dimension 的列表或元组

**返回**: 表示向量的 MultiVector

**异常**: `ValueError` - 如果数据长度与维度不匹配

**示例**:
```python
alg = nblade.Algebra.euclidean(3)

# 从列表创建
v = alg.vector([1.0, 2.0, 3.0])  # 1*e1 + 2*e2 + 3*e3

# 从元组创建
w = alg.vector((4.0, 5.0, 6.0))
```

---

### `basis_vector(i)`

创建第 i 个基向量 e_i。

```python
def basis_vector(self, i: int) -> MultiVector
```

**参数**: `i` - 基向量索引 (0 到 dimension-1)

**返回**: 基向量 MultiVector

**示例**:
```python
alg = nblade.Algebra.euclidean(3)
e1 = alg.basis_vector(0)  # 第一个基向量
e2 = alg.basis_vector(1)  # 第二个基向量
e3 = alg.basis_vector(2)  # 第三个基向量
```

---

### `basis_vectors()`

获取所有基向量。

```python
def basis_vectors(self) -> List[MultiVector]
```

**返回**: 基向量列表

**示例**:
```python
alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()
```

---

### `scalar(value)`

创建标量多向量。

```python
def scalar(self, value: float) -> MultiVector
```

**参数**: `value` - 标量值

**返回**: 标量 MultiVector

**示例**:
```python
s = alg.scalar(5.0)  # 纯标量
```

---

### `one()`

创建单位多向量 (标量 1)。

```python
def one(self) -> MultiVector
```

**返回**: 单位标量 MultiVector

---

### `zeros()`

创建零多向量。

```python
def zeros(self) -> MultiVector
```

**返回**: 零值 MultiVector

---

### `from_coefficients(coefficients)`

从系数数组创建多向量。

```python
def from_coefficients(self, coefficients: List[float]) -> MultiVector
```

**参数**: `coefficients` - 所有基的系数数组，长度为 2^dimension

**返回**: MultiVector

**示例**:
```python
# 3D 代数的系数顺序: [1, e1, e2, e3, e12, e13, e23, e123]
coeffs = [1.0, 2.0, 3.0, 4.0, 0.5, 0.0, 0.0, 0.1]
mv = alg.from_coefficients(coeffs)
```

---

## 转子创建

### `rotor(plane, angle)`

创建转子 (Rotator)。

```python
def rotor(self, plane: MultiVector, angle: float) -> MultiVector
```

**参数**:
- `plane` - 旋转平面 (二重向量)
- `angle` - 旋转角度 (弧度)

**返回**: 转子 MultiVector

**示例**:
```python
import math

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

# 在 xy 平面旋转 45 度
plane = e1 ^ e2
rotor = alg.rotor(plane, math.pi / 4)

# 使用转子旋转向量
rotated = e1.rotate_by(rotor)
```

---

## 字符串表示

### `__repr__()`

返回详细的代数信息。

```python
alg = nblade.Algebra.euclidean(3)
print(repr(alg))  # Algebra(dimension=3, signature=(3, 0, 0))
```

### `__str__()`

返回人类可读的代数描述。

```python
alg = nblade.Algebra.euclidean(3)
print(str(alg))  # Geometric Algebra G(3, 0, 0) (3D)
```

---

## 完整示例

```python
import nblade
import math

# 创建 3D 欧几里得代数
alg = nblade.Algebra.euclidean(3)

# 获取基向量
e1, e2, e3 = alg.basis_vectors()

# 创建向量
v = alg.vector([1.0, 2.0, 3.0])

# 基本运算
geometric = e1 * e2  # 几何积
outer = e1 ^ e2      # 外积
inner = e1 | e2      # 内积

# 创建转子并旋转
rotor = alg.rotor(e1 ^ e2, math.pi / 4)
rotated = v.rotate_by(rotor)

print(f"原始向量: {v}")
print(f"旋转后: {rotated}")
```