//! 范数和标量积模块
//!
//! 实现：
//! - 标量积：A * B = ⟨A†B⟩ (文档公式 71)
//! - 范数：|A|² = A * A (文档公式 85)
//!
//! # 公式总结
//!
//! | 公式 | 描述 | 文档位置 |
//! |------|------|---------|
//! | A * B = ⟨A†B⟩ | 标量积定义 | 公式 (71) |
//! | |A|² = A * A | 范数平方 | 公式 (85) |
//! | |A|² = |A*|² | 阶次对合保持范数 | 公式 (86) |
//! | |A|² = |A†|² | 反转保持范数 | 公式 (87) |
//! | |A|² = |A‡|² | Clifford 共轭保持范数 | 公式 (88) |

use crate::multiplication_table::GeometricProductTable;
use crate::multivector::{DenseMultiVector, MultiVector, SparseMultiVector};
use crate::operations::involution::reversion_dense;
use crate::products::geometric::geometric_product_dense;

/// 密集多向量的标量积
///
/// # 公式
/// A * B = ⟨A†B⟩ (文档公式 71)
///
/// # 参数
/// * `a` - 左操作数
/// * `b` - 右操作数
///
/// # 返回
/// 标量积结果
///
/// # 性质
/// - 对称性：A * B = B * A (公式 77)
/// - 对合不变：A * B = A* * B* = A† * B† = A‡ * B‡ (公式 78-80)
pub fn scalar_product_dense(a: &DenseMultiVector, b: &DenseMultiVector) -> f64 {
    // 计算 A†
    let a_rev = reversion_dense(a);

    // 计算 A†B (几何积)
    let table = GeometricProductTable::new(&a.config);
    let product = geometric_product_dense(&a_rev, b, &table);

    // 取标量部分 ⟨A†B⟩
    product.coefficients[0]
}

/// 稀疏多向量的标量积
///
/// 对于稀疏表示，使用等效公式：
/// A * B = Σ_I A†_I B_I (对所有基 I 求和)
pub fn scalar_product_sparse(a: &SparseMultiVector, b: &SparseMultiVector) -> f64 {
    let mut sum = 0.0;

    for (&idx, &coef_a) in &a.coefficients {
        if let Some(&coef_b) = b.coefficients.get(&idx) {
            // 计算 A† 的系数：需要应用反转符号
            // 公式：A†_r = (-1)^{r(r-1)/2} A_r
            let grade = a.config.grade_of_index(idx as usize);
            let sign_change = if grade == 0 {
                0
            } else {
                (grade * (grade - 1) / 2) % 2
            };
            let sign = if sign_change == 0 { 1.0 } else { -1.0 };
            sum += sign * coef_a * coef_b;
        }
    }

    sum
}

/// 密集多向量的范数平方
///
/// # 公式
/// |A|² = A * A (文档公式 85)
///
/// # 参数
/// * `mv` - 多向量
///
/// # 返回
/// 范数平方值（可能为负数，取决于签名）
pub fn norm_squared_dense(mv: &DenseMultiVector) -> f64 {
    scalar_product_dense(mv, mv)
}

/// 稀疏多向量的范数平方
pub fn norm_squared_sparse(mv: &SparseMultiVector) -> f64 {
    scalar_product_sparse(mv, mv)
}

/// 计算范数（可能为负数的平方根）
///
/// # 参数
/// * `mv` - 多向量
///
/// # 返回
/// 如果 norm² >= 0，返回 √(norm²)
/// 如果 norm² < 0，返回 -√(-norm²)
pub fn norm(mv: &MultiVector) -> f64 {
    let norm_sq = mv.norm_squared();
    if norm_sq >= 0.0 {
        norm_sq.sqrt()
    } else {
        -(-norm_sq).sqrt()
    }
}

impl MultiVector {
    /// 标量积 A * B
    ///
    /// # 公式
    /// A * B = ⟨A†B⟩ (文档公式 71)
    ///
    /// # 示例
    /// ```
    /// use nblade::{AlgebraConfig, MultiVector};
    /// use std::sync::Arc;
    ///
    /// let config = Arc::new(AlgebraConfig::euclidean(3));
    /// let e1 = MultiVector::basis_vector(config.clone(), 0);
    /// let e2 = MultiVector::basis_vector(config.clone(), 1);
    ///
    /// // 正交向量的标量积为 0
    /// assert!(e1.scalar_product(&e2).abs() < 1e-10);
    ///
    /// // 向量自身的标量积等于其平方
    /// assert!((e1.scalar_product(&e1) - 1.0).abs() < 1e-10);
    /// ```
    pub fn scalar_product(&self, other: &Self) -> f64 {
        match (self, other) {
            (MultiVector::Dense(a), MultiVector::Dense(b)) => scalar_product_dense(a, b),
            (MultiVector::Sparse(a), MultiVector::Sparse(b)) => scalar_product_sparse(a, b),
            (MultiVector::Dense(a), MultiVector::Sparse(b)) => {
                scalar_product_dense(a, &b.to_dense())
            }
            (MultiVector::Sparse(a), MultiVector::Dense(b)) => {
                scalar_product_dense(&a.to_dense(), b)
            }
        }
    }

    /// 范数平方 |A|²
    ///
    /// # 公式
    /// |A|² = A * A (文档公式 85)
    ///
    /// # 示例
    /// ```
    /// use nblade::{AlgebraConfig, MultiVector};
    /// use std::sync::Arc;
    ///
    /// let config = Arc::new(AlgebraConfig::euclidean(3));
    /// let e1 = MultiVector::basis_vector(config.clone(), 0);
    ///
    /// // |e1|² = 1
    /// assert!((e1.norm_squared() - 1.0).abs() < 1e-10);
    ///
    /// // |2e1|² = 4
    /// let scaled = e1.scale(2.0);
    /// assert!((scaled.norm_squared() - 4.0).abs() < 1e-10);
    /// ```
    pub fn norm_squared(&self) -> f64 {
        match self {
            MultiVector::Dense(d) => norm_squared_dense(d),
            MultiVector::Sparse(s) => norm_squared_sparse(s),
        }
    }

    /// 范数 |A|
    ///
    /// # 示例
    /// ```
    /// use nblade::{AlgebraConfig, MultiVector};
    /// use std::sync::Arc;
    ///
    /// let config = Arc::new(AlgebraConfig::euclidean(3));
    /// let e1 = MultiVector::basis_vector(config.clone(), 0);
    ///
    /// // |e1| = 1
    /// assert!((e1.norm() - 1.0).abs() < 1e-10);
    /// ```
    pub fn norm(&self) -> f64 {
        norm(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AlgebraConfig, MultiVector};
    use std::sync::Arc;

    fn test_config() -> crate::AlgebraConfigRef {
        Arc::new(AlgebraConfig::euclidean(3))
    }

    #[test]
    fn test_scalar_product_vectors() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        // e1 * e1 = 1
        assert!((e1.scalar_product(&e1) - 1.0).abs() < 1e-10);

        // e1 * e2 = 0 (正交)
        assert!(e1.scalar_product(&e2).abs() < 1e-10);
    }

    #[test]
    fn test_scalar_product_symmetry() {
        // 公式 (77): A * B = B * A
        let config = test_config();
        let mut a = MultiVector::zeros(config.clone());
        a = a.add(&MultiVector::from_scalar(config.clone(), 1.0));
        a = a.add(&MultiVector::basis_vector(config.clone(), 0));

        let mut b = MultiVector::zeros(config.clone());
        b = b.add(&MultiVector::basis_vector(config.clone(), 1));
        b = b.add(&MultiVector::basis_vector(config.clone(), 2));

        let ab = a.scalar_product(&b);
        let ba = b.scalar_product(&a);
        assert!((ab - ba).abs() < 1e-10);
    }

    #[test]
    fn test_norm_squared() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        // |e1|² = 1
        assert!((e1.norm_squared() - 1.0).abs() < 1e-10);

        // |2e1|² = 4
        let scaled = e1.scale(2.0);
        assert!((scaled.norm_squared() - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_norm() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        // |e1| = 1
        assert!((e1.norm() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_norm_squared_bivectors() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);

        // |e12|² = 1 (欧几里得空间)
        assert!((e12.norm_squared() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_norm_squared_trivector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);
        let e123 = e1.outer_product(&e2).outer_product(&e3);

        // |e123|² = 1 (欧几里得空间)
        assert!((e123.norm_squared() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_norm_squared_mixed() {
        let config = test_config();
        let scalar = MultiVector::from_scalar(config.clone(), 3.0);
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let mv = scalar.add(&e1.scale(4.0));

        // |3 + 4e1|² = 3² + 4² = 25
        assert!((mv.norm_squared() - 25.0).abs() < 1e-10);
    }

    #[test]
    fn test_scalar_product_grade_orthogonality() {
        // 不同阶次的多重向量标量积为 0
        let config = test_config();
        let scalar = MultiVector::from_scalar(config.clone(), 1.0);
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);

        // 标量 * 向量 = 0
        assert!(scalar.scalar_product(&e1).abs() < 1e-10);

        // 标量 * 二重向量 = 0
        assert!(scalar.scalar_product(&e12).abs() < 1e-10);

        // 向量 * 二重向量 = 0
        assert!(e1.scalar_product(&e12).abs() < 1e-10);
    }

    #[test]
    fn test_norm_invariance_under_involution() {
        // 公式 (86-88): |A|² = |A*|² = |A†|² = |A‡|²
        let config = test_config();
        let mut mv = MultiVector::zeros(config.clone());
        mv = mv.add(&MultiVector::from_scalar(config.clone(), 1.0));
        mv = mv.add(&MultiVector::basis_vector(config.clone(), 0));
        mv = mv.add(&MultiVector::basis_vector(config.clone(), 1));
        mv = mv.add(
            &MultiVector::basis_vector(config.clone(), 0)
                .outer_product(&MultiVector::basis_vector(config.clone(), 1)),
        );

        let norm_sq = mv.norm_squared();
        let norm_sq_star = mv.grade_involution().norm_squared();
        let norm_sq_dagger = mv.reversion().norm_squared();
        let norm_sq_conj = mv.clifford_conjugate().norm_squared();

        assert!((norm_sq - norm_sq_star).abs() < 1e-10, "阶次对合应保持范数");
        assert!((norm_sq - norm_sq_dagger).abs() < 1e-10, "反转应保持范数");
        assert!(
            (norm_sq - norm_sq_conj).abs() < 1e-10,
            "Clifford 共轭应保持范数"
        );
    }

    #[test]
    fn test_scalar_product_spacetime() {
        // 测试时空代数 G(1,3,0)
        let config = Arc::new(AlgebraConfig::new(4, crate::Signature::new(1, 3, 0)));
        let e0 = MultiVector::basis_vector(config.clone(), 0);
        let e1 = MultiVector::basis_vector(config.clone(), 1);

        // e0² = +1
        assert!((e0.scalar_product(&e0) - 1.0).abs() < 1e-10);

        // e1² = -1
        assert!((e1.scalar_product(&e1) + 1.0).abs() < 1e-10);
    }
}
