//! 对偶运算模块 / Dual Operations Module
//!
//! 实现：
//! - 对偶：A⊥ = A·I (右几何积) / Dual: A⊥ = A·I (right geometric product)
//! - 逆对偶：A⁻⊥ = I⁻¹·A (左几何积) / Inverse Dual: A⁻⊥ = I⁻¹·A (left geometric product)
//!
//! 其中 I 是体积元（最高阶次基向量） / where I is the volume element (highest grade basis)

use crate::multivector::{DenseMultiVector, MultiVector, SparseMultiVector};
use crate::products::inner::{left_inner_dense, left_inner_sparse};

/// 获取体积元索引
///
/// 体积元 I = e₁∧e₂∧...∧eₙ 的索引是所有位都为 1
#[inline]
fn volume_element_index(dimension: u32) -> u64 {
    (1u64 << dimension) - 1
}

/// 创建体积元多向量 I = e₁∧e₂∧...∧eₙ
fn create_volume_element(config: &crate::AlgebraConfigRef) -> MultiVector {
    let n = config.dimension();
    let idx = volume_element_index(n);

    match config.basis_count() {
        size if size <= 1024 => {
            // 小维度使用密集表示
            let mut coeffs = vec![0.0f64; config.basis_count()];
            coeffs[idx as usize] = 1.0;
            MultiVector::Dense(DenseMultiVector::from_coefficients(config.clone(), coeffs))
        }
        _ => {
            // 大维度使用稀疏表示
            let mut coeffs = std::collections::HashMap::new();
            coeffs.insert(idx, 1.0);
            MultiVector::Sparse(SparseMultiVector::from_map(config.clone(), coeffs))
        }
    }
}

/// 计算体积元的逆 I⁻¹
///
/// # 公式
/// I⁻¹ = I†/|I|²
///
/// 对于欧几里得空间，|I|² = ±1，所以 I⁻¹ = ±I
#[allow(dead_code)]
fn volume_element_inverse(config: &crate::AlgebraConfigRef) -> MultiVector {
    let n = config.dimension();
    let idx = volume_element_index(n) as usize;

    // I² = (-1)^{n(n-1)/2 + q} 其中 q 是负平方基向量数
    // 对于欧几里得空间 G(n,0,0)：I² = (-1)^{n(n-1)/2}
    let n_neg = config.signature.negative;
    let sign_exp = (n * (n - 1) / 2) + n_neg;
    let sign = if sign_exp % 2 == 0 { 1.0 } else { -1.0 };

    // I⁻¹ = I/I² = I * sign (因为 |I|² = ±1)
    match config.basis_count() {
        size if size <= 1024 => {
            let mut coeffs = vec![0.0f64; config.basis_count()];
            coeffs[idx] = sign;
            MultiVector::Dense(DenseMultiVector::from_coefficients(config.clone(), coeffs))
        }
        _ => {
            let mut coeffs = std::collections::HashMap::new();
            coeffs.insert(idx as u64, sign);
            MultiVector::Sparse(SparseMultiVector::from_map(config.clone(), coeffs))
        }
    }
}

/// 密集多向量的对偶
///
/// A⊥ = A†·I = (-1)^{r(r-1)/2} A·I (右几何积配反转符号)
#[allow(non_snake_case)]
pub fn dual_dense(mv: &DenseMultiVector) -> DenseMultiVector {
    let config = mv.config.clone();
    let I = create_volume_element(&config);

    let I_dense = match I {
        MultiVector::Dense(d) => d,
        MultiVector::Sparse(s) => s.to_dense(),
    };

    let mv_rev = crate::operations::involution::reversion_dense(mv);

    let mv_mv = MultiVector::Dense(mv_rev);
    let I_mv = MultiVector::Dense(I_dense);
    let result = mv_mv.geometric_product(&I_mv);

    match result {
        MultiVector::Dense(d) => d,
        MultiVector::Sparse(s) => s.to_dense(),
    }
}

/// 稀疏多向量的对偶
#[allow(non_snake_case)]
pub fn dual_sparse(mv: &SparseMultiVector) -> SparseMultiVector {
    let config = mv.config.clone();
    let I = create_volume_element(&config);

    let I_sparse = match I {
        MultiVector::Sparse(s) => s,
        MultiVector::Dense(d) => SparseMultiVector::from_dense(&d),
    };

    let mv_rev = crate::operations::involution::reversion_sparse(mv);

    let mv_mv = MultiVector::Sparse(mv_rev);
    let I_mv = MultiVector::Sparse(I_sparse);
    let result = mv_mv.geometric_product(&I_mv);

    match result {
        MultiVector::Sparse(s) => s,
        MultiVector::Dense(d) => SparseMultiVector::from_dense(&d),
    }
}

/// 密集多向量的逆对偶
///
/// # 公式
/// A⁻⊥ = A⌋I
#[allow(non_snake_case)]
pub fn inverse_dual_dense(mv: &DenseMultiVector) -> DenseMultiVector {
    let config = mv.config.clone();
    let n = config.dimension();
    let idx = volume_element_index(n) as usize;

    let mut I = DenseMultiVector::zeros(config.clone());
    I.coefficients[idx] = 1.0;

    left_inner_dense(mv, &I)
}

/// 稀疏多向量的逆对偶
#[allow(non_snake_case)]
pub fn inverse_dual_sparse(mv: &SparseMultiVector) -> SparseMultiVector {
    let config = mv.config.clone();
    let n = config.dimension();
    let idx = volume_element_index(n);

    let mut I = SparseMultiVector::zeros(config.clone());
    I.coefficients.insert(idx, 1.0);

    left_inner_sparse(mv, &I)
}

impl MultiVector {
    /// 对偶 A⊥ / Dual A⊥
    ///
    /// A⊥ = A·I (右几何积) / A⊥ = A·I (right geometric product)
    ///
    /// # 示例 / Example
    /// ```
    /// use nblade::{AlgebraConfig, MultiVector};
    /// use std::sync::Arc;
    ///
    /// let config = Arc::new(AlgebraConfig::euclidean(3));
    /// let e1 = MultiVector::basis_vector(config.clone(), 0);
    ///
    /// // e1⊥ = e2∧e3 (在 3D 中) / e1⊥ = e2∧e3 (in 3D)
    /// let dual = e1.dual();
    /// ```
    pub fn dual(&self) -> Self {
        match self {
            MultiVector::Dense(d) => MultiVector::Dense(dual_dense(d)),
            MultiVector::Sparse(s) => MultiVector::Sparse(dual_sparse(s)),
        }
    }

    /// 逆对偶 A⁻⊥ / Inverse Dual A⁻⊥
    ///
    /// # 公式 / Formula
    /// A⁻⊥ = A⌋I (左内积与体积元素) / A⁻⊥ = A⌋I (left inner product with volume element)
    ///
    /// # 示例 / Example
    /// ```
    /// use nblade::{AlgebraConfig, MultiVector};
    /// use std::sync::Arc;
    ///
    /// let config = Arc::new(AlgebraConfig::euclidean(3));
    /// let e1 = MultiVector::basis_vector(config.clone(), 0);
    ///
    /// // e1⁻⊥ = e23 (在 3D 中) / e1⁻⊥ = e23 (in 3D)
    /// let inv_dual = e1.inverse_dual();
    /// ```
    pub fn inverse_dual(&self) -> Self {
        match self {
            MultiVector::Dense(d) => MultiVector::Dense(inverse_dual_dense(d)),
            MultiVector::Sparse(s) => MultiVector::Sparse(inverse_dual_sparse(s)),
        }
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
    fn test_dual_scalar() {
        let config = test_config();
        let scalar = MultiVector::from_scalar(config.clone(), 1.0);

        // 1⊥ = e123 (在 3D 中)
        let dual = scalar.dual();
        assert!((dual.get_coefficient(7) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_dual_vector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        // e1⊥ = e23 (在 3D 中)
        let dual = e1.dual();
        assert!((dual.get_coefficient(6) - 1.0).abs() < 1e-10); // e23 = 0b110 = 6
    }

    #[test]
    fn test_dual_bivector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);

        // (e1∧e2)⊥ = e3 (在 3D 中)
        let dual = e12.dual();
        assert!((dual.get_coefficient(4) - 1.0).abs() < 1e-10); // e3 = 0b100 = 4
    }

    #[test]
    fn test_dual_trivector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);
        let e123 = e1.outer_product(&e2).outer_product(&e3);

        // (e1∧e2∧e3)⊥ = 1 (在 3D 中)
        let dual = e123.dual();
        assert!((dual.get_coefficient(0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_dual_double() {
        // 双重对偶：A⊥⊥ = ±A
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        let dual = e1.dual();
        let dual_dual = dual.dual();

        // 在 3D 欧几里得空间中，e1⊥⊥ = e1
        assert!((dual_dual.get_coefficient(1) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_inverse_dual() {
        let config = test_config();
        let scalar = MultiVector::from_scalar(config.clone(), 1.0);

        // 1⁻⊥ = I (volume element) because scalar ⌋ anything = that thing
        let inv_dual = scalar.inverse_dual();
        // I = e123 in 3D, coefficient at index 7
        assert!((inv_dual.get_coefficient(7) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_dual_4d() {
        // 测试 4D 空间的对偶
        let config = Arc::new(AlgebraConfig::euclidean(4));
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        // e1⊥ = e2∧e3∧e4 (在 4D 中)
        let dual = e1.dual();
        // e234 = 0b1110 = 14
        assert!((dual.get_coefficient(14) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_dual_bivector_4d() {
        // 测试 4D 空间二重向量的对偶
        let config = Arc::new(AlgebraConfig::euclidean(4));
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);

        // (e1∧e2)⊥ = e3∧e4 (在 4D 中)
        let dual = e12.dual();
        // e34 = 0b1100 = 12
        assert!((dual.get_coefficient(12) - 1.0).abs() < 1e-10);
    }

    // ==================== inverse_dual 综合测试 ====================

    #[test]
    fn test_inverse_dual_scalar() {
        // 标量的逆对偶: 1⁻⊥ = 1⌋I = I (体积元)
        let config = test_config();
        let scalar = MultiVector::from_scalar(config.clone(), 2.0);

        let inv_dual = scalar.inverse_dual();
        // 2⌋I = 2I (在 3D 中，I = e123，索引为 7)
        assert!((inv_dual.get_coefficient(7) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_inverse_dual_vector() {
        // 向量的逆对偶: e1⁻⊥ = e1⌋I
        // 在 3D 中，I = e123，e1⌋e123 = e23 (索引 6 = 0b110)
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        let inv_dual = e1.inverse_dual();
        // e1⌋e123 = e23
        assert!((inv_dual.get_coefficient(6) - 1.0).abs() < 1e-10); // e23 = 0b110 = 6
    }

    #[test]
    fn test_inverse_dual_bivector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);

        let inv_dual = e12.inverse_dual();
        assert!((inv_dual.get_coefficient(4) + 1.0).abs() < 1e-10); // e12⌋e123 = -e3
    }

    #[test]
    fn test_inverse_dual_trivector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);
        let e123 = e1.outer_product(&e2).outer_product(&e3);

        let inv_dual = e123.inverse_dual();
        assert!((inv_dual.get_coefficient(0) + 1.0).abs() < 1e-10); // e123⌋e123 = -1
    }

    #[test]
    fn test_inverse_dual_4d_vector() {
        // 4D 空间向量的逆对偶
        // e1⁻⊥ = e1⌋e1234 = e234 (索引 14 = 0b1110)
        let config = Arc::new(AlgebraConfig::euclidean(4));
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        let inv_dual = e1.inverse_dual();
        assert!((inv_dual.get_coefficient(14) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_inverse_dual_bivector_4d() {
        let config = Arc::new(AlgebraConfig::euclidean(4));
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);

        let inv_dual = e12.inverse_dual();
        assert!((inv_dual.get_coefficient(12) + 1.0).abs() < 1e-10); // e12⌋e1234 = -e34
    }

    #[test]
    fn test_dual_inverse_dual_relationship() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        // e1⊥ = e23, e23⁻⊥ = e23⌋e123 = -e1
        let dual_e1 = e1.dual();
        let inv_dual_of_dual = dual_e1.inverse_dual();

        // e23⌋e123 = -e1
        assert!((inv_dual_of_dual.get_coefficient(1) + 1.0).abs() < 1e-10);

        // Reverse: inverse_dual then dual
        let inv_dual_e1 = e1.inverse_dual(); // e1⌋e123 = e23
        let _dual_of_inv_dual = inv_dual_e1.dual();
    }

    #[test]
    fn test_inverse_dual_zero() {
        // 零向量的逆对偶应该为零
        let config = test_config();
        let zero = MultiVector::zeros(config);

        let inv_dual = zero.inverse_dual();
        assert!(inv_dual.is_zero());
    }

    #[test]
    fn test_inverse_dual_multivector() {
        // 测试混合多向量的逆对偶
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let scalar = MultiVector::from_scalar(config.clone(), 3.0);
        let combined = e1.add(&scalar);

        // (e1 + 3)⌋I = e1⌋I + 3⌋I = e23 + 3*e123
        let inv_dual = combined.inverse_dual();
        assert!((inv_dual.get_coefficient(6) - 1.0).abs() < 1e-10); // e23 部分
        assert!((inv_dual.get_coefficient(7) - 3.0).abs() < 1e-10); // 3*e123 部分
    }

    #[test]
    fn test_inverse_dual_sparse() {
        // 测试稀疏表示的逆对偶
        let config = Arc::new(AlgebraConfig::euclidean(3));
        // 创建一个稀疏多向量（低密度）
        let mut coeffs = vec![0.0; 8];
        coeffs[1] = 1.0; // e1

        let sparse_mv = MultiVector::from_coefficients(config, coeffs);
        let inv_dual = sparse_mv.inverse_dual();

        // e1⌋e123 = e23
        assert!((inv_dual.get_coefficient(6) - 1.0).abs() < 1e-10);
    }
}
