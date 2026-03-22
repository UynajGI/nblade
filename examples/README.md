# nblade 示例目录

本目录包含几何代数库的各种示例，按功能分类组织：

## 目录结构

```
examples/
├── tutorials/       # 教程示例（推荐入门）
│   ├── 01_quickstart.py        # 快速入门
│   ├── 02_basic_operations.py  # 基本运算
│   ├── 03_vectors.py           # 向量操作
│   ├── 05_rotations.py         # 旋转操作
│   └── 06_reciprocal_frame.py  # 互逆标架
├── basic/           # 基础运算示例
│   ├── 01_basic_operations.py  # 基本运算演示
│   ├── 02_vector_creation.py   # 向量创建演示
│   └── vector_from_array.py    # 从数组创建向量
├── advanced/        # 高级运算示例
│   └── 02_advanced_operations.py
├── applications/    # 应用示例
│   └── 03_applications.py      # 几何应用演示
├── physics/         # 物理应用示例
│   └── 01_rigid_body.py        # 刚体动力学
└── cg/              # 计算机图形学示例
    └── 01_transformations.py   # 几何变换
```

## 示例说明

### tutorials/ (推荐入门)
- **01_quickstart.py**: 快速入门，5 分钟上手 nblade
- **02_basic_operations.py**: 几何积、外积、内积等基本运算
- **03_vectors.py**: 向量创建和操作
- **05_rotations.py**: 使用转子进行旋转
- **06_reciprocal_frame.py**: 互逆标架、度量张量、向量分解

### basic/
- **01_basic_operations.py**: 演示基本几何代数运算
- **02_vector_creation.py**: 演示从 Python 列表/数组创建向量

### advanced/
- **02_advanced_operations.py**: 演示高级运算（对合、逆、对偶、交换子等）

### applications/
- **03_applications.py**: 演示实际应用（旋转、反射、投影等几何变换）

### physics/
- **01_rigid_body.py**: 刚体动力学模拟示例

### cg/
- **01_transformations.py**: 计算机图形学中的几何变换

## 运行示例

```bash
# 运行教程示例
python examples/tutorials/01_quickstart.py

# 运行基础示例
python examples/basic/01_basic_operations.py

# 运行高级示例
python examples/advanced/02_advanced_operations.py

# 运行应用示例
python examples/applications/03_applications.py
```

## 依赖

所有示例都需要安装 nblade 库：

```bash
pip install maturin
maturin develop  # 在项目根目录运行
```