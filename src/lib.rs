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

pub mod signature;
pub mod algebra_config;
pub mod basis;
pub mod multivector;
pub mod multiplication_table;
pub mod products;
pub mod operations;
pub mod geometry;

pub use signature::Signature;
pub use algebra_config::{AlgebraConfig, AlgebraConfigRef};
pub use multivector::{MultiVector, DenseMultiVector, SparseMultiVector};
pub use multiplication_table::MultiplicationTables;

// 重新导出常用模块
pub mod prelude {
    pub use crate::Signature;
    pub use crate::AlgebraConfig;
    pub use crate::MultiVector;
    pub use crate::basis::grade_of_index;
    pub use crate::basis::index::index_to_string;
}

// PyO3 Python 绑定
#[cfg(feature = "python")]
pub mod python;
