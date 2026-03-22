//! GA-Rust: 任意维度通用几何代数库
//!
//! 高性能几何代数计算库，支持：
//! - 任意维度 n（最多 64 维）
//! - 任意签名 (p, q, r)
//! - 密集/稀疏双表示
//! - 并行计算优化
//!
//! # 示例
//!
//! ```rust
//! use nblade::{AlgebraConfig, Signature, MultiVector};
//! use std::sync::Arc;
//!
//! // 创建 3D 欧几里得几何代数
//! let config = Arc::new(AlgebraConfig::euclidean(3));
//!
//! // 创建基向量
//! let e1 = MultiVector::basis_vector(config.clone(), 0);
//! let e2 = MultiVector::basis_vector(config.clone(), 1);
//!
//! // 几何积
//! let product = e1.geometric_product(&e2);
//!
//! // 外积
//! let wedge = e1.outer_product(&e2);
//!
//! // 左内积
//! let inner = e1.left_inner(&e2);
//! ```

// Clippy allows for common patterns in this codebase
#![allow(clippy::needless_range_loop)]
#![allow(clippy::manual_is_multiple_of)]
#![allow(clippy::unnecessary_cast)]

pub mod algebra_config;
pub mod basis;
pub mod geometry;
pub mod multiplication_table;
pub mod multivector;
pub mod operations;
pub mod products;
pub mod signature;

pub use algebra_config::{AlgebraConfig, AlgebraConfigRef};
pub use multiplication_table::MultiplicationTables;
pub use multivector::{DenseMultiVector, MultiVector, SparseMultiVector};
pub use signature::Signature;

// 重新导出常用模块
pub mod prelude {
    pub use crate::basis::grade_of_index;
    pub use crate::basis::index::index_to_string;
    pub use crate::AlgebraConfig;
    pub use crate::MultiVector;
    pub use crate::Signature;
}

// PyO3 Python 绑定
#[cfg(feature = "python")]
pub mod python;
