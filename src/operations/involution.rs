//! 对合运算模块
//!
//! 实现三种对合运算：
//! - 阶次对合 (Grade Involution): A* = (-1)^r A_r
//! - 反转 (Reversion): A† = (-1)^{r(r-1)/2} A_r
//! - Clifford 共轭：A‡ = A*†

use crate::multivector::{DenseMultiVector, MultiVector, SparseMultiVector};

/// 密集多向量的阶次对合
pub fn grade_involution_dense(mv: &DenseMultiVector) -> DenseMultiVector {
    let mut result = mv.coefficients.clone();

    for (i, coef) in result.iter_mut().enumerate() {
        let grade = mv.config.grade_of_index(i);
        if grade % 2 == 1 {
            *coef = -*coef;
        }
    }

    DenseMultiVector {
        config: mv.config.clone(),
        coefficients: result,
    }
}

/// 稀疏多向量的阶次对合
pub fn grade_involution_sparse(mv: &SparseMultiVector) -> SparseMultiVector {
    let mut result = mv.coefficients.clone();

    for (&idx, coef) in result.iter_mut() {
        let grade = mv.config.grade_of_index(idx as usize);
        if grade % 2 == 1 {
            *coef = -*coef;
        }
    }

    SparseMultiVector {
        config: mv.config.clone(),
        coefficients: result,
    }
}

/// 密集多向量的反转
pub fn reversion_dense(mv: &DenseMultiVector) -> DenseMultiVector {
    let mut result = mv.coefficients.clone();

    for (i, coef) in result.iter_mut().enumerate() {
        let grade = mv.config.grade_of_index(i);
        // r(r-1)/2 的奇偶性
        let sign_change = if grade == 0 {
            0
        } else {
            (grade * (grade - 1) / 2) % 2
        };
        if sign_change == 1 {
            *coef = -*coef;
        }
    }

    DenseMultiVector {
        config: mv.config.clone(),
        coefficients: result,
    }
}

/// 稀疏多向量的反转
pub fn reversion_sparse(mv: &SparseMultiVector) -> SparseMultiVector {
    let mut result = mv.coefficients.clone();

    for (&idx, coef) in result.iter_mut() {
        let grade = mv.config.grade_of_index(idx as usize);
        let sign_change = if grade == 0 {
            0
        } else {
            (grade * (grade - 1) / 2) % 2
        };
        if sign_change == 1 {
            *coef = -*coef;
        }
    }

    SparseMultiVector {
        config: mv.config.clone(),
        coefficients: result,
    }
}

/// Clifford 共轭 = 阶次对合 + 反转
pub fn clifford_conjugate_dense(mv: &DenseMultiVector) -> DenseMultiVector {
    grade_involution_dense(&reversion_dense(mv))
}

pub fn clifford_conjugate_sparse(mv: &SparseMultiVector) -> SparseMultiVector {
    grade_involution_sparse(&reversion_sparse(mv))
}

impl MultiVector {
    /// 阶次对合 A*
    pub fn grade_involution(&self) -> Self {
        match self {
            MultiVector::Dense(d) => MultiVector::Dense(grade_involution_dense(d)),
            MultiVector::Sparse(s) => MultiVector::Sparse(grade_involution_sparse(s)),
        }
    }

    /// 反转 A†
    pub fn reversion(&self) -> Self {
        match self {
            MultiVector::Dense(d) => MultiVector::Dense(reversion_dense(d)),
            MultiVector::Sparse(s) => MultiVector::Sparse(reversion_sparse(s)),
        }
    }

    /// Clifford 共轭 A‡
    pub fn clifford_conjugate(&self) -> Self {
        match self {
            MultiVector::Dense(d) => MultiVector::Dense(clifford_conjugate_dense(d)),
            MultiVector::Sparse(s) => MultiVector::Sparse(clifford_conjugate_sparse(s)),
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
    fn test_grade_involution() {
        let config = test_config();
        let mut mv = MultiVector::zeros(config.clone());
        mv = mv.add(&MultiVector::from_scalar(config.clone(), 1.0)); // 标量 (偶)
        mv = mv.add(&MultiVector::basis_vector(config.clone(), 0)); // 向量 (奇)
        mv = mv.add(&MultiVector::basis_vector(config.clone(), 1)); // 向量 (奇)

        let result = mv.grade_involution();

        // 标量部分不变
        assert!((result.get_coefficient(0) - 1.0).abs() < 1e-10);
        // 向量部分变号
        assert!((result.get_coefficient(1) + 1.0).abs() < 1e-10);
        assert!((result.get_coefficient(2) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_reversion() {
        let config = test_config();

        // 向量反转不变：e1† = e1
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let result = e1.reversion();
        assert!((result.get_coefficient(1) - 1.0).abs() < 1e-10);

        // 二重向量反转变号：(e1∧e2)† = -e1∧e2
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);
        let result = e12.reversion();
        assert!((result.get_coefficient(3) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_clifford_conjugate() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        // e1‡ = -e1 (阶次对合变号，反转不变)
        let result = e1.clifford_conjugate();
        assert!((result.get_coefficient(1) + 1.0).abs() < 1e-10);
    }

    // ========== Sparse版本测试 ==========

    #[test]
    fn test_grade_involution_sparse() {
        let config = test_config();

        // 创建稀疏多向量：标量 + e1 + e2
        let sparse = SparseMultiVector::from_scalar(config.clone(), 1.0);
        let e1 = SparseMultiVector::basis_vector(config.clone(), 0);
        let e2 = SparseMultiVector::basis_vector(config.clone(), 1);

        let mut coefficients = sparse.coefficients.clone();
        for (&k, &v) in &e1.coefficients {
            coefficients.insert(k, v);
        }
        for (&k, &v) in &e2.coefficients {
            coefficients.insert(k, v);
        }
        let mv = SparseMultiVector::from_map(config.clone(), coefficients);

        let result = grade_involution_sparse(&mv);

        // 标量部分不变 (grade 0 is even)
        assert!((result.get_coefficient(0) - 1.0).abs() < 1e-10);
        // 向量部分变号 (grade 1 is odd)
        assert!((result.get_coefficient(1) + 1.0).abs() < 1e-10);
        assert!((result.get_coefficient(2) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_reversion_sparse() {
        let config = test_config();

        // 向量反转不变：e1† = e1
        let e1 = SparseMultiVector::basis_vector(config.clone(), 0);
        let result = reversion_sparse(&e1);
        assert!((result.get_coefficient(1) - 1.0).abs() < 1e-10);

        // 二重向量反转变号：(e1∧e2)† = -e1∧e2
        // 对于 grade 2: r(r-1)/2 = 2*1/2 = 1 (odd, sign changes)
        let e1_mv = MultiVector::basis_vector(config.clone(), 0);
        let e2_mv = MultiVector::basis_vector(config.clone(), 1);
        let e12_mv = e1_mv.outer_product(&e2_mv);
        let e12_sparse = match e12_mv {
            MultiVector::Dense(d) => SparseMultiVector::from_dense(&d),
            MultiVector::Sparse(s) => s,
        };

        let result = reversion_sparse(&e12_sparse);
        assert!((result.get_coefficient(3) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_clifford_conjugate_sparse() {
        let config = test_config();

        // e1‡ = -e1
        let e1 = SparseMultiVector::basis_vector(config.clone(), 0);
        let result = clifford_conjugate_sparse(&e1);
        assert!((result.get_coefficient(1) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_clifford_conjugate_dense() {
        let config = test_config();

        // e1‡ = -e1 (阶次对合变号，反转不变)
        let e1 = DenseMultiVector::basis_vector(config.clone(), 0);
        let result = clifford_conjugate_dense(&e1);
        assert!((result.get_coefficient(1) + 1.0).abs() < 1e-10);

        // 标量不变
        let scalar = DenseMultiVector::from_scalar(config.clone(), 5.0);
        let result = clifford_conjugate_dense(&scalar);
        assert!((result.scalar_part() - 5.0).abs() < 1e-10);
    }

    // ========== 高阶次测试 ==========

    #[test]
    fn test_grade_involution_trivector() {
        let config = test_config();

        // 三重向量 e1∧e2∧e3 (grade 3, odd, sign changes)
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);
        let e12 = e1.outer_product(&e2);
        let e123 = e12.outer_product(&e3);

        let result = e123.grade_involution();
        // grade 3 is odd, sign should change
        assert!((result.get_coefficient(7) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_reversion_trivector() {
        let config = test_config();

        // 三重向量 (e1∧e2∧e3)†
        // grade 3: r(r-1)/2 = 3*2/2 = 3 (odd, sign changes)
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);
        let e12 = e1.outer_product(&e2);
        let e123 = e12.outer_product(&e3);

        let result = e123.reversion();
        assert!((result.get_coefficient(7) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_clifford_conjugate_trivector() {
        let config = test_config();

        // 三重向量 e1∧e2∧e3
        // grade_involution: sign changes (grade 3 odd)
        // reversion: sign changes (3*2/2 = 3 odd)
        // Combined: sign*sign = sign (one change)
        // Actually clifford_conjugate = grade_involution(reversion(mv))
        // reversion changes sign, then grade_involution changes sign again
        // For grade 3: reversion sign = -1, then grade_involution on grade 3 = -1
        // Result: (-1)*(-1) = +1
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);
        let e12 = e1.outer_product(&e2);
        let e123 = e12.outer_product(&e3);

        let result = e123.clifford_conjugate();
        // grade 3: reversion(-1) then grade_involution(-1) = +1
        assert!((result.get_coefficient(7) - 1.0).abs() < 1e-10);
    }

    // ========== MultiVector枚举的Sparse分支测试 ==========

    #[test]
    fn test_multivector_sparse_grade_involution() {
        let config = test_config();

        // 创建一个稀疏多向量
        let sparse = SparseMultiVector::basis_vector(config.clone(), 0);
        let mv = MultiVector::Sparse(sparse);

        let result = mv.grade_involution();
        match result {
            MultiVector::Sparse(s) => {
                // grade 1 is odd, sign changes
                assert!((s.get_coefficient(1) + 1.0).abs() < 1e-10);
            }
            MultiVector::Dense(_) => panic!("Expected Sparse variant"),
        }
    }

    #[test]
    fn test_multivector_sparse_reversion() {
        let config = test_config();

        // 创建稀疏多向量：e1∧e2 (bivector)
        let e1_mv = MultiVector::basis_vector(config.clone(), 0);
        let e2_mv = MultiVector::basis_vector(config.clone(), 1);
        let e12_mv = e1_mv.outer_product(&e2_mv);
        let e12_sparse = match e12_mv {
            MultiVector::Dense(d) => SparseMultiVector::from_dense(&d),
            MultiVector::Sparse(s) => s,
        };
        let mv = MultiVector::Sparse(e12_sparse);

        let result = mv.reversion();
        match result {
            MultiVector::Sparse(s) => {
                // grade 2: r(r-1)/2 = 1 (odd, sign changes)
                assert!((s.get_coefficient(3) + 1.0).abs() < 1e-10);
            }
            MultiVector::Dense(_) => panic!("Expected Sparse variant"),
        }
    }

    #[test]
    fn test_multivector_sparse_clifford_conjugate() {
        let config = test_config();

        let sparse = SparseMultiVector::basis_vector(config.clone(), 0);
        let mv = MultiVector::Sparse(sparse);

        let result = mv.clifford_conjugate();
        match result {
            MultiVector::Sparse(s) => {
                // e1‡ = -e1
                assert!((s.get_coefficient(1) + 1.0).abs() < 1e-10);
            }
            MultiVector::Dense(_) => panic!("Expected Sparse variant"),
        }
    }

    // ========== 边界情况测试 ==========

    #[test]
    fn test_grade_involution_empty_sparse() {
        let config = test_config();

        // 空稀疏多向量
        let empty = SparseMultiVector::zeros(config.clone());
        let result = grade_involution_sparse(&empty);
        assert!(result.coefficients.is_empty());
    }

    #[test]
    fn test_reversion_empty_sparse() {
        let config = test_config();

        let empty = SparseMultiVector::zeros(config.clone());
        let result = reversion_sparse(&empty);
        assert!(result.coefficients.is_empty());
    }

    #[test]
    fn test_clifford_conjugate_empty_sparse() {
        let config = test_config();

        let empty = SparseMultiVector::zeros(config.clone());
        let result = clifford_conjugate_sparse(&empty);
        assert!(result.coefficients.is_empty());
    }

    #[test]
    fn test_grade_involution_scalar_only() {
        let config = test_config();

        // 仅标量
        let scalar = DenseMultiVector::from_scalar(config.clone(), 3.5);
        let result = grade_involution_dense(&scalar);
        assert!((result.scalar_part() - 3.5).abs() < 1e-10);
    }

    #[test]
    fn test_reversion_scalar_only() {
        let config = test_config();

        let scalar = DenseMultiVector::from_scalar(config.clone(), 7.0);
        let result = reversion_dense(&scalar);
        assert!((result.scalar_part() - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_mixed_grades() {
        let config = test_config();

        // 混合阶次：标量 + 向量 + 二重向量 + 三重向量
        let scalar = MultiVector::from_scalar(config.clone(), 1.0);
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);
        let e12 = e1.outer_product(&e2);
        let e123 = e12.outer_product(&e3);

        let mv = scalar.add(&e1).add(&e12).add(&e123);

        // 验证各阶次的阶次对合
        let result = mv.grade_involution();
        // 标量 (grade 0): 不变
        assert!((result.get_coefficient(0) - 1.0).abs() < 1e-10);
        // 向量 (grade 1): 变号
        assert!((result.get_coefficient(1) + 1.0).abs() < 1e-10);
        // 二重向量 (grade 2): 不变
        assert!((result.get_coefficient(3) - 1.0).abs() < 1e-10);
        // 三重向量 (grade 3): 变号
        assert!((result.get_coefficient(7) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_reversion_mixed_grades() {
        let config = test_config();

        // 混合阶次
        let scalar = MultiVector::from_scalar(config.clone(), 2.0);
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);
        let e12 = e1.outer_product(&e2);
        let e123 = e12.outer_product(&e3);

        let mv = scalar.add(&e1).add(&e12).add(&e123);

        let result = mv.reversion();
        // 标量 (grade 0): 不变
        assert!((result.get_coefficient(0) - 2.0).abs() < 1e-10);
        // 向量 (grade 1): r(r-1)/2 = 0, 不变
        assert!((result.get_coefficient(1) - 1.0).abs() < 1e-10);
        // 二重向量 (grade 2): r(r-1)/2 = 1, 变号
        assert!((result.get_coefficient(3) + 1.0).abs() < 1e-10);
        // 三重向量 (grade 3): r(r-1)/2 = 3, 变号
        assert!((result.get_coefficient(7) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_involution_double_application() {
        let config = test_config();

        // 双重应用阶次对合应恢复原值
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);

        let result1 = e12.grade_involution();
        let result2 = result1.grade_involution();

        // grade 2: 第一次不变，第二次也不变
        assert!((result2.get_coefficient(3) - 1.0).abs() < 1e-10);

        // 向量：第一次变号，第二次恢复
        let result1 = e1.grade_involution();
        let result2 = result1.grade_involution();
        assert!((result2.get_coefficient(1) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_reversion_double_application() {
        let config = test_config();

        // 双重应用反转会恢复原值
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);

        let result1 = e12.reversion();
        let result2 = result1.reversion();

        assert!((result2.get_coefficient(3) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_clifford_conjugate_double_application() {
        let config = test_config();

        // 双重应用Clifford共轭会恢复原值
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        let result1 = e1.clifford_conjugate();
        let result2 = result1.clifford_conjugate();

        assert!((result2.get_coefficient(1) - 1.0).abs() < 1e-10);
    }
}
