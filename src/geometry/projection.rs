//! 投影模块
//!
//! 实现：
//! - 正交投影：P_Ar(B) = (B⌋Ar)⌋Ar⁻¹
//! - 正交拒绝：R_Ar(B) = (B∧Ar)·Ar⁻¹ (使用几何积，而非左内积)

use crate::multivector::MultiVector;

/// 正交投影到子空间
///
/// 公式：P_Ar(B) = (B⌋Ar)⌋Ar⁻¹
pub fn project_to(b: &MultiVector, blade: &MultiVector) -> Result<MultiVector, String> {
    // 计算 blade 的逆
    let blade_inv = blade.inverse()?;

    // B⌋Ar
    let inner = b.left_inner(blade);

    // (B⌋Ar)⌋Ar⁻¹
    let result = inner.left_inner(&blade_inv);

    Ok(result)
}

/// 正交拒绝（投影到正交补）
///
/// 公式：R_Ar(B) = (B∧Ar)·Ar⁻¹
pub fn reject_from(b: &MultiVector, blade: &MultiVector) -> Result<MultiVector, String> {
    // 计算 blade 的逆
    let blade_inv = blade.inverse()?;

    // B∧Ar
    let outer = b.outer_product(blade);

    // (B∧Ar)·Ar⁻¹ (使用几何积)
    let result = outer.geometric_product(&blade_inv);

    Ok(result)
}

impl MultiVector {
    /// 正交投影到由 blade 表示的子空间
    pub fn project_to(&self, blade: &Self) -> Result<Self, String> {
        project_to(self, blade)
    }

    /// 正交拒绝（投影到正交补空间）
    pub fn reject_from(&self, blade: &Self) -> Result<Self, String> {
        reject_from(self, blade)
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
    fn test_project_vector_to_vector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        // e1 投影到 e1 上 = e1
        let proj = e1.project_to(&e1).unwrap();
        assert!((proj.get_coefficient(1) - 1.0).abs() < 1e-10);

        // e2 投影到 e1 上 = 0
        let proj = e2.project_to(&e1).unwrap();
        assert!(proj.is_zero());
    }

    #[test]
    fn test_project_vector_to_plane() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);
        let e12 = e1.outer_product(&e2);

        // e3 投影到 e12 平面 = 0
        let proj = e3.project_to(&e12).unwrap();
        assert!(proj.is_zero());

        // e1 投影到 e12 平面 = e1
        let proj = e1.project_to(&e12).unwrap();
        assert!((proj.get_coefficient(1) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_reject_vector_from_vector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        // e2 拒绝从 e1 = e2 (e2 已经正交于 e1)
        let rej = e2.reject_from(&e1).unwrap();
        assert!((rej.get_coefficient(2) - 1.0).abs() < 1e-10);

        // e1 拒绝从 e1 = 0
        let rej = e1.reject_from(&e1).unwrap();
        assert!(rej.is_zero());
    }

    #[test]
    fn test_project_plus_reject() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);

        // v = proj + reject
        let v = e1
            .add(&e2)
            .add(&MultiVector::basis_vector(config.clone(), 2));
        let proj = v.project_to(&e12).unwrap();
        let rej = v.reject_from(&e12).unwrap();

        let sum = proj.add(&rej);
        // 验证 v = proj + reject
        for i in 0..config.basis_count() {
            assert!((sum.get_coefficient(i) - v.get_coefficient(i)).abs() < 1e-10);
        }
    }
}
