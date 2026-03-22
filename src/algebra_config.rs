//! 代数配置模块 / Algebra Configuration Module
//!
//! 定义几何代数的完整配置，包括维度、签名和预计算的度量因子
//!
//! Defines the complete configuration for geometric algebra, including dimension,
//! signature, and precomputed metric factors.

use crate::multivector::MultiVector;
use crate::signature::Signature;
use std::sync::Arc;

/// 几何代数配置 / Geometric Algebra Configuration
///
/// 包含维度、签名和预计算的度量因子
/// 使用 Arc 实现线程安全的共享
///
/// Contains dimension, signature, and precomputed metric factors.
/// Uses Arc for thread-safe sharing.
pub type AlgebraConfigRef = Arc<AlgebraConfig>;

/// 几何代数的完整配置
#[derive(Debug, Clone)]
pub struct AlgebraConfig {
    /// 度量签名
    pub signature: Signature,
    /// 向量空间维度 n
    dimension: u32,
    /// 基向量数量 2^n
    basis_count: usize,
    /// 预计算的度量因子（用于加速几何积）
    metric_factors: Vec<f64>,
}

impl AlgebraConfig {
    /// 创建新的代数配置
    ///
    /// # 参数
    /// * `dimension` - 向量空间维度 n
    /// * `signature` - 度量签名
    ///
    /// # 示例
    /// ```
    /// use nblade::{AlgebraConfig, Signature};
    ///
    /// // 3D 欧几里得几何代数 G(3,0,0)
    /// let config = AlgebraConfig::new(3, Signature::euclidean(3));
    ///
    /// // 4D 时空代数 G(1,3,0)
    /// let spacetime = AlgebraConfig::new(4, Signature::spacetime(4));
    /// ```
    pub fn new(dimension: u32, signature: Signature) -> Self {
        assert!(
            signature.dimension() == dimension,
            "Signature dimension {} must match algebra dimension {}",
            signature.dimension(),
            dimension
        );

        let basis_count = 1 << dimension;
        let metric_factors = Self::precompute_metric_factors(&signature);

        Self {
            signature,
            dimension,
            basis_count,
            metric_factors,
        }
    }

    /// 创建欧几里得几何代数配置
    pub fn euclidean(dimension: u32) -> Self {
        Self::new(dimension, Signature::euclidean(dimension))
    }

    /// 获取向量空间维度 n
    #[inline]
    pub fn dimension(&self) -> u32 {
        self.dimension
    }

    /// 获取基向量数量 2^n
    #[inline]
    pub fn basis_count(&self) -> usize {
        self.basis_count
    }

    /// 获取预计算的度量因子
    #[inline]
    pub fn metric_factors(&self) -> &[f64] {
        &self.metric_factors
    }

    /// 获取第 index 个基的度量因子
    #[inline]
    pub fn get_metric_factor(&self, index: usize) -> f64 {
        self.metric_factors[index]
    }

    /// 预计算所有基向量的度量因子
    ///
    /// 对于基索引 i，度量因子是所有设置位对应的基向量平方值的乘积
    fn precompute_metric_factors(signature: &Signature) -> Vec<f64> {
        let n = signature.dimension() as usize;
        let mut factors = Vec::with_capacity(1 << n);

        for i in 0..(1 << n) {
            let mut factor = 1.0;
            for bit in 0..n {
                if i & (1 << bit) != 0 {
                    factor *= signature.basis_square_sign(bit as u32);
                    if factor == 0.0 {
                        break;
                    }
                }
            }
            factors.push(factor);
        }

        factors
    }

    /// 计算基索引的阶次（popcount）
    #[inline]
    pub fn grade_of_index(&self, index: usize) -> u32 {
        (index as u64).count_ones()
    }

    /// 检查索引是否有效
    #[inline]
    pub fn is_valid_index(&self, index: usize) -> bool {
        index < self.basis_count
    }

    /// 创建体积元素（伪标量）I = e₁∧e₂∧...∧eₙ / Create Volume Element (Pseudoscalar) I
    ///
    /// 体积元素是最高阶的单位 n-刃，基索引为所有位都为 1（即 2^n - 1）
    ///
    /// The volume element is the unit n-blade of highest grade with basis index
    /// where all bits are set (i.e., 2^n - 1).
    ///
    /// # 示例
    /// ```
    /// use nblade::{AlgebraConfig, MultiVector};
    /// use std::sync::Arc;
    ///
    /// let config = Arc::new(AlgebraConfig::euclidean(3));
    /// let I = config.volume_element();
    ///
    /// // I = e₁∧e₂∧e₃，基索引为 0b111 = 7
    /// assert!((I.get_coefficient(7) - 1.0).abs() < 1e-10);
    /// ```
    pub fn volume_element(&self) -> MultiVector {
        // 体积元素的基索引 = 所有位都为 1 = 2^n - 1
        let index = self.basis_count - 1;
        let mut coeffs = vec![0.0; self.basis_count];
        coeffs[index] = 1.0;
        MultiVector::from_coefficients(Arc::new(self.clone()), coeffs)
    }

    /// 计算体积元素的平方 I² / Compute Volume Element Squared I²
    ///
    /// 对于签名 G(p,q,r)，有：
    /// - I² = (-1)^(n(n-1)/2 + q)，其中 n = p + q + r
    ///
    /// For signature G(p,q,r):
    /// - I² = (-1)^(n(n-1)/2 + q) where n = p + q + r
    ///
    /// # 返回 / Returns
    /// - `Some(value)` 如果签名非退化（r = 0）/ if signature is non-degenerate (r = 0)
    /// - `None` 如果签名退化（r > 0，即存在零向量）/ if signature is degenerate (r > 0, null vectors exist)
    ///
    /// # 示例
    /// ```
    /// use nblade::AlgebraConfig;
    /// use std::sync::Arc;
    ///
    /// // 3D 欧几里得：I² = (-1)^(3*2/2 + 0) = (-1)^3 = -1
    /// let config = Arc::new(AlgebraConfig::euclidean(3));
    /// assert_eq!(config.volume_element_squared(), Some(-1.0));
    ///
    /// // 2D 欧几里得：I² = (-1)^(2*1/2 + 0) = (-1)^1 = -1
    /// let config2d = Arc::new(AlgebraConfig::euclidean(2));
    /// assert_eq!(config2d.volume_element_squared(), Some(-1.0));
    ///
    /// // 4D 欧几里得：I² = (-1)^(4*3/2 + 0) = (-1)^6 = 1
    /// let config4d = Arc::new(AlgebraConfig::euclidean(4));
    /// assert_eq!(config4d.volume_element_squared(), Some(1.0));
    /// ```
    pub fn volume_element_squared(&self) -> Option<f64> {
        // 如果存在零向量（r > 0），则 I² = 0，逆不存在
        if self.signature.null > 0 {
            return None;
        }

        let n = self.dimension as i64;
        let q = self.signature.negative as i64;

        // I² = (-1)^(n(n-1)/2 + q)
        let exponent = n * (n - 1) / 2 + q;

        if exponent % 2 == 0 {
            Some(1.0)
        } else {
            Some(-1.0)
        }
    }

    /// 创建体积元素的逆 I⁻¹ / Create Volume Element Inverse I⁻¹
    ///
    /// 由于 I² = ±1（对于非退化签名），所以 I⁻¹ = I / I² = I * I²
    ///
    /// Since I² = ±1 (for non-degenerate signatures), I⁻¹ = I / I² = I * I²
    ///
    /// # 返回 / Returns
    /// - `Some(I⁻¹)` 如果签名非退化 / if signature is non-degenerate
    /// - `None` 如果签名退化（存在零向量）/ if signature is degenerate (null vectors exist)
    ///
    /// # 示例
    /// ```
    /// use nblade::{AlgebraConfig, MultiVector};
    /// use std::sync::Arc;
    ///
    /// let config = Arc::new(AlgebraConfig::euclidean(3));
    /// let I = config.volume_element();
    /// let I_inv = config.volume_element_inverse().unwrap();
    ///
    /// // I * I⁻¹ = 1
    /// let product = I.geometric_product(&I_inv);
    /// assert!((product.scalar_part() - 1.0).abs() < 1e-10);
    /// ```
    #[allow(non_snake_case)]
    pub fn volume_element_inverse(&self) -> Option<MultiVector> {
        match self.volume_element_squared() {
            Some(i_squared) => {
                let I = self.volume_element();
                Some(I.scale(i_squared))
            }
            None => None,
        }
    }
}

impl PartialEq for AlgebraConfig {
    fn eq(&self, other: &Self) -> bool {
        self.dimension == other.dimension && self.signature == other.signature
    }
}

impl Eq for AlgebraConfig {}

impl std::hash::Hash for AlgebraConfig {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.dimension.hash(state);
        self.signature.hash(state);
    }
}

// 实现 Send 和 Sync 以支持多线程共享
unsafe impl Send for AlgebraConfig {}
unsafe impl Sync for AlgebraConfig {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_algebra_config_creation() {
        let config = AlgebraConfig::euclidean(3);
        assert_eq!(config.dimension(), 3);
        assert_eq!(config.basis_count(), 8);
    }

    #[test]
    fn test_metric_factors_euclidean() {
        let config = AlgebraConfig::euclidean(2);
        // G(2,0,0): 基为 [], [e1], [e2], [e1e2]
        // 度量因子：1, 1, 1, 1
        assert_eq!(config.metric_factors.len(), 4);
        assert_eq!(config.metric_factors[0], 1.0); // 标量
        assert_eq!(config.metric_factors[1], 1.0); // e1
        assert_eq!(config.metric_factors[2], 1.0); // e2
        assert_eq!(config.metric_factors[3], 1.0); // e1e2
    }

    #[test]
    fn test_metric_factors_spacetime() {
        let config = AlgebraConfig::new(2, Signature::new(1, 1, 0));
        // G(1,1,0): e0²=+1, e1²=-1
        // 基为：[], [e0], [e1], [e0e1]
        // 度量因子：1, 1, -1, -1
        assert_eq!(config.metric_factors.len(), 4);
        assert_eq!(config.metric_factors[0], 1.0); // 标量
        assert_eq!(config.metric_factors[1], 1.0); // e0
        assert_eq!(config.metric_factors[2], -1.0); // e1
        assert_eq!(config.metric_factors[3], -1.0); // e0e1
    }

    #[test]
    fn test_grade_of_index() {
        let config = AlgebraConfig::euclidean(3);
        assert_eq!(config.grade_of_index(0), 0); // 0b000 -> 标量
        assert_eq!(config.grade_of_index(1), 1); // 0b001 -> 向量
        assert_eq!(config.grade_of_index(3), 2); // 0b011 -> 二重向量
        assert_eq!(config.grade_of_index(7), 3); // 0b111 -> 三重向量
    }

    #[test]
    fn test_arc_sharing() {
        let config1 = Arc::new(AlgebraConfig::euclidean(3));
        let config2 = Arc::clone(&config1);

        assert_eq!(config1.dimension(), config2.dimension());
        assert_eq!(Arc::strong_count(&config1), 2);
    }

    // ========== 体积元素测试 ==========

    #[test]
    fn test_volume_element_2d() {
        let config = Arc::new(AlgebraConfig::euclidean(2));
        let I = config.volume_element();

        // I = e1∧e2，基索引 0b11 = 3
        assert_eq!(I.get_coefficient(3), 1.0);
        assert!((I.get_coefficient(0) - 0.0).abs() < 1e-15);
        assert!((I.get_coefficient(1) - 0.0).abs() < 1e-15);
        assert!((I.get_coefficient(2) - 0.0).abs() < 1e-15);
    }

    #[test]
    fn test_volume_element_3d() {
        let config = Arc::new(AlgebraConfig::euclidean(3));
        let I = config.volume_element();

        // I = e1∧e2∧e3，基索引 0b111 = 7
        assert!((I.get_coefficient(7) - 1.0).abs() < 1e-15);
        assert!((I.get_coefficient(0) - 0.0).abs() < 1e-15);

        // 检查阶次为 3
        assert_eq!(config.grade_of_index(7), 3);
    }

    #[test]
    fn test_volume_element_4d() {
        let config = Arc::new(AlgebraConfig::euclidean(4));
        let I = config.volume_element();

        // I = e1∧e2∧e3∧e4，基索引 0b1111 = 15
        assert!((I.get_coefficient(15) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_volume_element_squared_euclidean() {
        // 2D: I² = (-1)^(2*1/2 + 0) = (-1)^1 = -1
        let config2d = AlgebraConfig::euclidean(2);
        assert_eq!(config2d.volume_element_squared(), Some(-1.0));

        // 3D: I² = (-1)^(3*2/2 + 0) = (-1)^3 = -1
        let config3d = AlgebraConfig::euclidean(3);
        assert_eq!(config3d.volume_element_squared(), Some(-1.0));

        // 4D: I² = (-1)^(4*3/2 + 0) = (-1)^6 = 1
        let config4d = AlgebraConfig::euclidean(4);
        assert_eq!(config4d.volume_element_squared(), Some(1.0));

        // 5D: I² = (-1)^(5*4/2 + 0) = (-1)^10 = 1
        let config5d = AlgebraConfig::euclidean(5);
        assert_eq!(config5d.volume_element_squared(), Some(1.0));

        // 6D: I² = (-1)^(6*5/2 + 0) = (-1)^15 = -1
        let config6d = AlgebraConfig::euclidean(6);
        assert_eq!(config6d.volume_element_squared(), Some(-1.0));
    }

    #[test]
    fn test_volume_element_squared_spacetime() {
        // Spacetime G(1,3,0): n=4, q=3
        // I² = (-1)^(4*3/2 + 3) = (-1)^(6+3) = (-1)^9 = -1
        let config = AlgebraConfig::new(4, Signature::spacetime(4));
        assert_eq!(config.volume_element_squared(), Some(-1.0));
    }

    #[test]
    fn test_volume_element_squared_degenerate() {
        // 退化签名 G(2,0,1): 存在零向量
        let config = AlgebraConfig::new(3, Signature::new(2, 0, 1));
        assert_eq!(config.volume_element_squared(), None);
    }

    #[test]
    fn test_volume_element_inverse_3d() {
        let config = Arc::new(AlgebraConfig::euclidean(3));
        let I = config.volume_element();
        let I_inv = config.volume_element_inverse().unwrap();

        // I² = -1，所以 I⁻¹ = -I
        assert!((I_inv.get_coefficient(7) - (-1.0)).abs() < 1e-15);

        // I * I⁻¹ = 1
        let product = I.geometric_product(&I_inv);
        assert!((product.scalar_part() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_volume_element_inverse_4d() {
        let config = Arc::new(AlgebraConfig::euclidean(4));
        let I = config.volume_element();
        let I_inv = config.volume_element_inverse().unwrap();

        // I² = 1，所以 I⁻¹ = I
        assert!((I_inv.get_coefficient(15) - 1.0).abs() < 1e-15);

        // I * I⁻¹ = 1
        let product = I.geometric_product(&I_inv);
        assert!((product.scalar_part() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_volume_element_inverse_degenerate() {
        let config = Arc::new(AlgebraConfig::new(3, Signature::new(2, 0, 1)));
        assert!(config.volume_element_inverse().is_none());
    }

    #[test]
    fn test_volume_element_product_verification() {
        // 验证 I² 的几何积计算与公式结果一致
        let config = Arc::new(AlgebraConfig::euclidean(3));
        let I = config.volume_element();

        // 计算 I * I
        let I_squared_geom = I.geometric_product(&I);
        let expected = config.volume_element_squared().unwrap();

        // I * I 应该等于 I²（一个标量）
        assert!((I_squared_geom.scalar_part() - expected).abs() < 1e-10);
    }
}
