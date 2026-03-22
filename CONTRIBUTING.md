# Contributing to nblade

感谢你对 nblade 项目感兴趣！本文档将帮助你了解如何为项目做出贡献。

Thank you for your interest in nblade! This guide will help you contribute to the project.

---

## 开发环境设置 | Development Setup

### 前置要求 | Prerequisites

- **Rust**: 1.70+ ([安装指南](https://www.rust-lang.org/tools/install))
- **Python**: 3.8+ (用于 Python 绑定测试)
- **maturin**: Python 构建工具

### 安装步骤 | Installation Steps

```bash
# 1. Fork 并克隆仓库
git clone https://github.com/YOUR_USERNAME/nblade.git
cd nblade

# 2. 安装 Rust（如果尚未安装）
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# 3. 构建 Rust 库
cargo build --release

# 4. 安装 Python 开发工具
pip install maturin pytest ruff

# 5. 构建 Python 绑定
maturin develop --release

# 6. 运行测试
cargo test --all-features
pytest python/tests/
```

---

## 代码风格 | Code Style

### Rust

- 运行 `cargo fmt` 格式化代码
- 运行 `cargo clippy --all-features` 检查代码质量
- 所有公共 API 必须有文档注释

```bash
cargo fmt
cargo clippy --all-features -- -D warnings
```

### Python

- 遵循 PEP 8
- 使用 `ruff format` 格式化
- 使用 type hints

```bash
ruff format python/
ruff check python/
```

---

## 提交规范 | Commit Convention

使用 [Conventional Commits](https://www.conventionalcommits.org/) 格式：

```
<type>(<scope>): <description>

[optional body]

[optional footer]
```

### 类型 | Types

| Type | Description |
|------|-------------|
| `feat` | 新功能 |
| `fix` | Bug 修复 |
| `docs` | 文档更新 |
| `style` | 代码格式（不影响功能） |
| `refactor` | 重构 |
| `perf` | 性能优化 |
| `test` | 测试相关 |
| `chore` | 构建/工具相关 |

### 示例 | Examples

```
feat(multivector): add grade projection operation
fix(products): correct sign in geometric product
docs(api): add Chinese documentation for Algebra class
test(operations): add tests for dual operation
```

---

## Pull Request 流程 | Pull Request Process

1. **创建分支**: 从 `main` 创建功能分支
   ```bash
   git checkout -b feat/your-feature-name
   ```

2. **编写代码**: 遵循代码风格，添加测试

3. **运行测试**:
   ```bash
   cargo test --all-features
   pytest python/tests/
   cargo clippy --all-features -- -D warnings
   ```

4. **提交更改**: 使用规范的提交信息

5. **推送分支**:
   ```bash
   git push origin feat/your-feature-name
   ```

6. **创建 PR**: 在 GitHub 上创建 Pull Request

### PR 检查清单 | PR Checklist

- [ ] 代码通过所有测试
- [ ] 新功能有对应的测试
- [ ] 公共 API 有文档注释
- [ ] 提交信息符合规范
- [ ] 更新了相关文档

---

## 项目结构 | Project Structure

```
nblade/
├── src/               # Rust 源码
│   ├── lib.rs         # 入口文件
│   ├── multivector/   # 多重向量类型
│   ├── products/      # 几何积、外积、内积
│   ├── operations/    # 对偶、逆等操作
│   └── python/        # PyO3 绑定
├── python/            # Python 包
│   └── tests/         # Python 测试
├── examples/          # 示例代码
│   ├── tutorials/     # 教程示例
│   ├── physics/       # 物理应用
│   └── cg/            # 计算机图形
├── docs/              # 文档
│   ├── en/            # 英文文档
│   └── zh/            # 中文文档
└── tests/             # Rust 集成测试
```

---

## 报告问题 | Reporting Issues

在 [GitHub Issues](https://github.com/UynajGI/nblade/issues) 报告问题时，请包含：

- 操作系统和版本
- Rust 版本 (`rustc --version`)
- Python 版本 (`python --version`)
- 复现步骤
- 期望行为 vs 实际行为
- 相关日志或错误信息

---

## 行为准则 | Code of Conduct

- 尊重所有贡献者
- 保持专业和友好的交流
- 接受建设性批评
- 关注对社区最有利的事情

---

## 许可证 | License

通过贡献代码，你同意你的代码将以 MIT 许可证发布。

By contributing, you agree that your contributions will be licensed under the MIT License.