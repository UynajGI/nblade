//! 稀疏多向量模块
//!
//! 使用 HashMap 表示多向量，适用于高维代数或稀疏数据

use crate::algebra_config::{AlgebraConfig, AlgebraConfigRef};
use crate::basis::index::BasisIndex;
use std::collections::HashMap;

/// 稀疏多向量表示的阈值
/// 当密度低于此值时使用稀疏表示
pub const SPARSE_THRESHOLD: f64 = 0.1;

/// 稀疏多向量
///
/// 只存储非零系数，适用于高维代数
#[derive(Debug, Clone)]
pub struct SparseMultiVector {
    /// 代数配置（共享引用）
    pub config: AlgebraConfigRef,
    /// 非零系数映射：索引 -> 系数
    pub coefficients: HashMap<BasisIndex, f64>,
}

impl SparseMultiVector {
    /// 创建新的零多向量
    pub fn zeros(config: AlgebraConfigRef) -> Self {
        Self {
            config,
            coefficients: HashMap::new(),
        }
    }

    /// 创建单位多向量（标量部分为 1）
    pub fn one(config: AlgebraConfigRef) -> Self {
        let mut mv = Self::zeros(config);
        mv.coefficients.insert(0, 1.0);
        mv
    }

    /// 从 HashMap 创建多向量
    pub fn from_map(config: AlgebraConfigRef, coefficients: HashMap<BasisIndex, f64>) -> Self {
        Self {
            config,
            coefficients,
        }
    }

    /// 从密集多向量转换
    pub fn from_dense(dense: &crate::multivector::dense::DenseMultiVector) -> Self {
        let mut coefficients = HashMap::new();

        for (i, &coef) in dense.coefficients.iter().enumerate() {
            if coef.abs() > 1e-15 {
                coefficients.insert(i as BasisIndex, coef);
            }
        }

        Self {
            config: dense.config.clone(),
            coefficients,
        }
    }

    /// 转换为密集多向量
    pub fn to_dense(&self) -> crate::multivector::dense::DenseMultiVector {
        use crate::multivector::dense::DenseMultiVector;

        let mut coefficients = vec![0.0f64; self.config.basis_count()];
        for (&idx, &coef) in &self.coefficients {
            coefficients[idx as usize] = coef;
        }

        DenseMultiVector::from_coefficients(self.config.clone(), coefficients)
    }

    /// 从标量创建多向量
    pub fn from_scalar(config: AlgebraConfigRef, scalar: f64) -> Self {
        let mut mv = Self::zeros(config);
        if scalar.abs() > 1e-15 {
            mv.coefficients.insert(0, scalar);
        }
        mv
    }

    /// 创建基向量 e_i
    pub fn basis_vector(config: AlgebraConfigRef, i: u32) -> Self {
        let mut mv = Self::zeros(config);
        mv.coefficients.insert(1 << i, 1.0);
        mv
    }

    /// 获取标量部分
    #[inline]
    pub fn scalar_part(&self) -> f64 {
        *self.coefficients.get(&0).unwrap_or(&0.0)
    }

    /// 获取第 i 个系数
    #[inline]
    pub fn get_coefficient(&self, index: BasisIndex) -> f64 {
        *self.coefficients.get(&index).unwrap_or(&0.0)
    }

    /// 设置第 i 个系数
    #[inline]
    pub fn set_coefficient(&mut self, index: BasisIndex, value: f64) {
        if value.abs() > 1e-15 {
            self.coefficients.insert(index, value);
        } else {
            self.coefficients.remove(&index);
        }
    }

    /// 获取非零项的数量
    #[inline]
    pub fn non_zero_count(&self) -> usize {
        self.coefficients.len()
    }

    /// 获取密度（非零项比例）
    pub fn density(&self) -> f64 {
        self.non_zero_count() as f64 / self.config.basis_count() as f64
    }

    /// 检查是否应该使用稀疏表示
    pub fn should_use_sparse(&self) -> bool {
        self.density() < SPARSE_THRESHOLD
    }

    /// 检查是否为零多向量
    pub fn is_zero(&self) -> bool {
        self.coefficients.is_empty()
    }

    /// 获取阶次 r 的部分
    pub fn grade_projection(&self, r: u32) -> Self {
        let mut result = Self::zeros(self.config.clone());

        for (&idx, &coef) in &self.coefficients {
            if self.config.grade_of_index(idx as usize) == r {
                result.coefficients.insert(idx, coef);
            }
        }

        result
    }

    /// 获取偶部
    pub fn even_part(&self) -> Self {
        let mut result = Self::zeros(self.config.clone());

        for (&idx, &coef) in &self.coefficients {
            if self.config.grade_of_index(idx as usize) % 2 == 0 {
                result.coefficients.insert(idx, coef);
            }
        }

        result
    }

    /// 获取奇部
    pub fn odd_part(&self) -> Self {
        let mut result = Self::zeros(self.config.clone());

        for (&idx, &coef) in &self.coefficients {
            if self.config.grade_of_index(idx as usize) % 2 == 1 {
                result.coefficients.insert(idx, coef);
            }
        }

        result
    }

    /// 加法
    pub fn add(&self, other: &Self) -> Self {
        let mut result = self.clone();

        for (&idx, &coef) in &other.coefficients {
            let current = result.coefficients.entry(idx).or_insert(0.0);
            *current += coef;
            if current.abs() < 1e-15 {
                result.coefficients.remove(&idx);
            }
        }

        result
    }

    /// 减法
    pub fn sub(&self, other: &Self) -> Self {
        let mut result = self.clone();

        for (&idx, &coef) in &other.coefficients {
            let current = result.coefficients.entry(idx).or_insert(0.0);
            *current -= coef;
            if current.abs() < 1e-15 {
                result.coefficients.remove(&idx);
            }
        }

        result
    }

    /// 标量乘法
    pub fn scale(&self, scalar: f64) -> Self {
        if scalar.abs() < 1e-15 {
            return Self::zeros(self.config.clone());
        }

        let mut result = Self::zeros(self.config.clone());

        for (&idx, &coef) in &self.coefficients {
            let new_coef = coef * scalar;
            if new_coef.abs() > 1e-15 {
                result.coefficients.insert(idx, new_coef);
            }
        }

        result
    }

    /// 取负
    pub fn neg(&self) -> Self {
        let mut result = Self::zeros(self.config.clone());

        for (&idx, &coef) in &self.coefficients {
            result.coefficients.insert(idx, -coef);
        }

        result
    }

    /// 获取代数配置引用
    #[inline]
    pub fn config(&self) -> &AlgebraConfig {
        &self.config
    }

    /// 迭代所有非零项
    pub fn iter(&self) -> impl Iterator<Item = (&BasisIndex, &f64)> {
        self.coefficients.iter()
    }
}

impl std::fmt::Display for SparseMultiVector {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.coefficients.is_empty() {
            return write!(f, "0");
        }

        let mut terms: Vec<(BasisIndex, f64)> =
            self.coefficients.iter().map(|(&k, &v)| (k, v)).collect();
        terms.sort_by_key(|(k, _)| *k);

        let term_strings: Vec<String> = terms
            .iter()
            .filter(|(_, c)| c.abs() > 1e-15)
            .map(|(i, c)| {
                let grade = self.config.grade_of_index(*i as usize);
                let basis = crate::basis::index::index_to_string(*i, self.config.dimension());
                if grade == 0 {
                    format!("{:.6}", c)
                } else {
                    format!("{:.6}{}", c, basis)
                }
            })
            .collect();

        if term_strings.is_empty() {
            write!(f, "0")
        } else {
            write!(f, "{}", term_strings.join(" + "))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra_config::AlgebraConfig;
    use std::sync::Arc;

    fn test_config() -> AlgebraConfigRef {
        Arc::new(AlgebraConfig::euclidean(3))
    }

    #[test]
    fn test_sparse_zeros() {
        let config = test_config();
        let mv = SparseMultiVector::zeros(config.clone());
        assert!(mv.is_zero());
        assert_eq!(mv.non_zero_count(), 0);
    }

    #[test]
    fn test_sparse_one() {
        let config = test_config();
        let mv = SparseMultiVector::one(config.clone());
        assert_eq!(mv.scalar_part(), 1.0);
        assert_eq!(mv.non_zero_count(), 1);
    }

    #[test]
    fn test_sparse_basis_vector() {
        let config = test_config();
        let e1 = SparseMultiVector::basis_vector(config.clone(), 0);
        assert_eq!(e1.get_coefficient(1), 1.0);
    }

    #[test]
    fn test_sparse_density() {
        let config = test_config();
        let mv = SparseMultiVector::zeros(config.clone());
        assert_eq!(mv.density(), 0.0);

        let mv = SparseMultiVector::one(config.clone());
        assert!((mv.density() - 1.0 / 8.0).abs() < 1e-10);
    }

    #[test]
    fn test_sparse_add_sub() {
        let config = test_config();
        let a = SparseMultiVector::basis_vector(config.clone(), 0);
        let b = SparseMultiVector::basis_vector(config.clone(), 1);

        let sum = a.add(&b);
        assert_eq!(sum.get_coefficient(1), 1.0);
        assert_eq!(sum.get_coefficient(2), 1.0);
        assert_eq!(sum.non_zero_count(), 2);
    }

    #[test]
    fn test_sparse_scale() {
        let config = test_config();
        let mv = SparseMultiVector::basis_vector(config.clone(), 0);
        let scaled = mv.scale(2.0);
        assert_eq!(scaled.get_coefficient(1), 2.0);
    }

    #[test]
    fn test_sparse_grade_projection() {
        let config = test_config();
        let mut mv = SparseMultiVector::zeros(config.clone());
        mv.set_coefficient(0, 1.0); // 标量
        mv.set_coefficient(1, 2.0); // e1
        mv.set_coefficient(3, 3.0); // e12

        let scalar_part = mv.grade_projection(0);
        assert_eq!(scalar_part.get_coefficient(0), 1.0);
        assert_eq!(scalar_part.non_zero_count(), 1);
    }

    #[test]
    fn test_sparse_dense_conversion() {
        use crate::multivector::dense::DenseMultiVector;

        let config = test_config();
        let sparse = SparseMultiVector::basis_vector(config.clone(), 0);
        let dense = sparse.to_dense();

        assert_eq!(dense.get_coefficient(1), 1.0);

        let back_to_sparse = SparseMultiVector::from_dense(&dense);
        assert_eq!(back_to_sparse.get_coefficient(1), 1.0);
    }

    #[test]
    fn test_sparse_from_map() {
        let config = test_config();
        let mut map = HashMap::new();
        map.insert(0, 5.0);
        map.insert(1, 3.0);
        map.insert(3, 2.0); // e12

        let mv = SparseMultiVector::from_map(config.clone(), map);
        assert_eq!(mv.scalar_part(), 5.0);
        assert_eq!(mv.get_coefficient(1), 3.0);
        assert_eq!(mv.get_coefficient(3), 2.0);
        assert_eq!(mv.non_zero_count(), 3);
    }

    #[test]
    fn test_sparse_from_scalar() {
        let config = test_config();

        // Non-zero scalar
        let mv = SparseMultiVector::from_scalar(config.clone(), 7.0);
        assert_eq!(mv.scalar_part(), 7.0);
        assert_eq!(mv.non_zero_count(), 1);

        // Zero scalar (should be empty)
        let mv_zero = SparseMultiVector::from_scalar(config.clone(), 0.0);
        assert!(mv_zero.is_zero());
        assert_eq!(mv_zero.scalar_part(), 0.0);
    }

    #[test]
    fn test_sparse_set_coefficient_zero_removal() {
        let config = test_config();
        let mut mv = SparseMultiVector::basis_vector(config.clone(), 0);
        assert_eq!(mv.non_zero_count(), 1);

        // Set to zero should remove the coefficient
        mv.set_coefficient(1, 0.0);
        assert!(mv.is_zero());
        assert_eq!(mv.non_zero_count(), 0);
    }

    #[test]
    fn test_sparse_should_use_sparse() {
        let config = test_config();

        // Empty (0% density) - should use sparse
        let empty = SparseMultiVector::zeros(config.clone());
        assert!(empty.should_use_sparse());

        // 1/8 = 12.5% density, which is ABOVE SPARSE_THRESHOLD (10%)
        // So should_use_sparse() returns false
        let one = SparseMultiVector::one(config.clone());
        assert!(!one.should_use_sparse());

        // Create a sparse vector in higher dimension where 1 coeff is < 10%
        let config_5d = Arc::new(AlgebraConfig::euclidean(5));
        let one_5d = SparseMultiVector::one(config_5d.clone());
        // 1/32 = 3.125% < 10%, so should use sparse
        assert!(one_5d.should_use_sparse());
    }

    #[test]
    fn test_sparse_grade_projection_all_grades() {
        let config = test_config();
        let mut mv = SparseMultiVector::zeros(config.clone());
        mv.set_coefficient(0, 1.0); // grade 0: scalar
        mv.set_coefficient(1, 2.0); // grade 1: e1
        mv.set_coefficient(2, 3.0); // grade 1: e2
        mv.set_coefficient(3, 4.0); // grade 2: e12
        mv.set_coefficient(7, 5.0); // grade 3: e123 (pseudoscalar)

        // Grade 0
        let g0 = mv.grade_projection(0);
        assert_eq!(g0.get_coefficient(0), 1.0);
        assert_eq!(g0.non_zero_count(), 1);

        // Grade 1
        let g1 = mv.grade_projection(1);
        assert_eq!(g1.get_coefficient(1), 2.0);
        assert_eq!(g1.get_coefficient(2), 3.0);
        assert_eq!(g1.non_zero_count(), 2);

        // Grade 2
        let g2 = mv.grade_projection(2);
        assert_eq!(g2.get_coefficient(3), 4.0);
        assert_eq!(g2.non_zero_count(), 1);

        // Grade 3
        let g3 = mv.grade_projection(3);
        assert_eq!(g3.get_coefficient(7), 5.0);
        assert_eq!(g3.non_zero_count(), 1);

        // Grade 4 (empty in 3D)
        let g4 = mv.grade_projection(4);
        assert!(g4.is_zero());
    }

    #[test]
    fn test_sparse_even_part() {
        let config = test_config();
        let mut mv = SparseMultiVector::zeros(config.clone());
        mv.set_coefficient(0, 1.0); // grade 0 (even): scalar
        mv.set_coefficient(1, 2.0); // grade 1 (odd): e1
        mv.set_coefficient(3, 4.0); // grade 2 (even): e12

        let even = mv.even_part();
        assert_eq!(even.get_coefficient(0), 1.0);
        assert_eq!(even.get_coefficient(3), 4.0);
        assert_eq!(even.non_zero_count(), 2);
        assert_eq!(even.get_coefficient(1), 0.0); // odd part excluded
    }

    #[test]
    fn test_sparse_odd_part() {
        let config = test_config();
        let mut mv = SparseMultiVector::zeros(config.clone());
        mv.set_coefficient(0, 1.0); // grade 0 (even): scalar
        mv.set_coefficient(1, 2.0); // grade 1 (odd): e1
        mv.set_coefficient(3, 4.0); // grade 2 (even): e12
        mv.set_coefficient(7, 6.0); // grade 3 (odd): e123

        let odd = mv.odd_part();
        assert_eq!(odd.get_coefficient(1), 2.0);
        assert_eq!(odd.get_coefficient(7), 6.0);
        assert_eq!(odd.non_zero_count(), 2);
        assert_eq!(odd.get_coefficient(0), 0.0); // even part excluded
        assert_eq!(odd.get_coefficient(3), 0.0); // even part excluded
    }

    #[test]
    fn test_sparse_sub() {
        let config = test_config();
        let a = SparseMultiVector::basis_vector(config.clone(), 0);
        let b = SparseMultiVector::basis_vector(config.clone(), 1);

        let diff = a.sub(&b);
        assert_eq!(diff.get_coefficient(1), 1.0);
        assert_eq!(diff.get_coefficient(2), -1.0);
        assert_eq!(diff.non_zero_count(), 2);
    }

    #[test]
    fn test_sparse_sub_cancellation() {
        let config = test_config();
        let a = SparseMultiVector::basis_vector(config.clone(), 0);
        let same = SparseMultiVector::basis_vector(config.clone(), 0);

        // a - a = 0
        let result = a.sub(&same);
        assert!(result.is_zero());
        assert_eq!(result.non_zero_count(), 0);
    }

    #[test]
    fn test_sparse_add_cancellation() {
        let config = test_config();
        let mut a = SparseMultiVector::basis_vector(config.clone(), 0);
        let mut neg_a = SparseMultiVector::basis_vector(config.clone(), 0);
        neg_a = neg_a.scale(-1.0);

        // a + (-a) = 0
        let result = a.add(&neg_a);
        assert!(result.is_zero());
    }

    #[test]
    fn test_sparse_scale_zero() {
        let config = test_config();
        let mv = SparseMultiVector::basis_vector(config.clone(), 0);

        // Scale by zero should return empty multivector
        let scaled = mv.scale(0.0);
        assert!(scaled.is_zero());
        assert_eq!(scaled.non_zero_count(), 0);
    }

    #[test]
    fn test_sparse_neg() {
        let config = test_config();
        let mut mv = SparseMultiVector::zeros(config.clone());
        mv.set_coefficient(0, 1.0);
        mv.set_coefficient(1, 2.0);
        mv.set_coefficient(3, -3.0);

        let negated = mv.neg();
        assert_eq!(negated.get_coefficient(0), -1.0);
        assert_eq!(negated.get_coefficient(1), -2.0);
        assert_eq!(negated.get_coefficient(3), 3.0);
        assert_eq!(negated.non_zero_count(), 3);
    }

    #[test]
    fn test_sparse_neg_empty() {
        let config = test_config();
        let mv = SparseMultiVector::zeros(config.clone());
        let negated = mv.neg();
        assert!(negated.is_zero());
    }

    #[test]
    fn test_sparse_config() {
        let config = test_config();
        let mv = SparseMultiVector::zeros(config.clone());
        assert_eq!(mv.config().dimension(), 3);
    }

    #[test]
    fn test_sparse_iter() {
        let config = test_config();
        let mut mv = SparseMultiVector::zeros(config.clone());
        mv.set_coefficient(0, 1.0);
        mv.set_coefficient(1, 2.0);
        mv.set_coefficient(3, 3.0);

        let items: Vec<_> = mv.iter().collect();
        assert_eq!(items.len(), 3);

        // Check that we can iterate over all coefficients
        let sum: f64 = items.iter().map(|(_, &v)| v).sum();
        assert!((sum - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_sparse_display_empty() {
        let config = test_config();
        let mv = SparseMultiVector::zeros(config.clone());
        let display = format!("{}", mv);
        assert_eq!(display, "0");
    }

    #[test]
    fn test_sparse_display_scalar() {
        let config = test_config();
        let mv = SparseMultiVector::from_scalar(config.clone(), 3.5);
        let display = format!("{}", mv);
        assert!(display.contains("3.5"));
    }

    #[test]
    fn test_sparse_display_vector() {
        let config = test_config();
        let mv = SparseMultiVector::basis_vector(config.clone(), 0);
        let display = format!("{}", mv);
        assert!(display.contains("1.0"));
        assert!(display.contains("e1"));
    }

    #[test]
    fn test_sparse_display_multiple_terms() {
        let config = test_config();
        let mut mv = SparseMultiVector::zeros(config.clone());
        mv.set_coefficient(0, 1.0);
        mv.set_coefficient(1, 2.0);
        mv.set_coefficient(3, 3.0);

        let display = format!("{}", mv);
        assert!(display.contains(" + "));
        assert!(display.contains("e1"));
        assert!(display.contains("e1∧e2"));
    }

    #[test]
    fn test_sparse_from_dense_all_zeros() {
        use crate::multivector::dense::DenseMultiVector;

        let config = test_config();
        let dense = DenseMultiVector::zeros(config.clone());
        let sparse = SparseMultiVector::from_dense(&dense);

        assert!(sparse.is_zero());
        assert_eq!(sparse.non_zero_count(), 0);
    }

    #[test]
    fn test_sparse_from_dense_multiple_nonzero() {
        use crate::multivector::dense::DenseMultiVector;

        let config = test_config();
        let mut dense = DenseMultiVector::zeros(config.clone());
        dense.set_coefficient(0, 1.0);
        dense.set_coefficient(1, 2.0);
        dense.set_coefficient(7, 5.0); // pseudoscalar

        let sparse = SparseMultiVector::from_dense(&dense);

        assert_eq!(sparse.scalar_part(), 1.0);
        assert_eq!(sparse.get_coefficient(1), 2.0);
        assert_eq!(sparse.get_coefficient(7), 5.0);
        assert_eq!(sparse.non_zero_count(), 3);
    }

    #[test]
    fn test_sparse_to_dense_multiple_nonzero() {
        let config = test_config();
        let mut sparse = SparseMultiVector::zeros(config.clone());
        sparse.set_coefficient(0, 1.0);
        sparse.set_coefficient(1, 2.0);
        sparse.set_coefficient(7, 5.0);

        let dense = sparse.to_dense();

        assert_eq!(dense.get_coefficient(0), 1.0);
        assert_eq!(dense.get_coefficient(1), 2.0);
        assert_eq!(dense.get_coefficient(7), 5.0);
        // Other coefficients should be zero
        assert_eq!(dense.get_coefficient(2), 0.0);
        assert_eq!(dense.get_coefficient(3), 0.0);
    }

    #[test]
    fn test_sparse_get_coefficient_missing() {
        let config = test_config();
        let mv = SparseMultiVector::zeros(config.clone());

        // Get coefficient that doesn't exist
        assert_eq!(mv.get_coefficient(999), 0.0);
    }

    #[test]
    fn test_sparse_scalar_part_missing() {
        let config = test_config();
        let mv = SparseMultiVector::basis_vector(config.clone(), 0);

        // No scalar part
        assert_eq!(mv.scalar_part(), 0.0);
    }
}
