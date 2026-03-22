//! 密集多向量模块
//!
//! 使用密集数组表示多向量，适用于低维代数或稠密数据

use crate::algebra_config::{AlgebraConfig, AlgebraConfigRef};
use ndarray::Array1;

#[cfg(feature = "pool")]
use crate::multivector::pool::BufferPool;

/// 密集多向量表示
///
/// 使用固定大小的数组存储所有 2^n 个系数
#[derive(Debug, Clone)]
pub struct DenseMultiVector {
    /// 代数配置（共享引用）
    pub config: AlgebraConfigRef,
    /// 系数数组（长度为 2^n）
    pub coefficients: Array1<f64>,
}

impl DenseMultiVector {
    /// 创建新的零多向量
    pub fn zeros(config: AlgebraConfigRef) -> Self {
        let size = config.basis_count();
        Self {
            config,
            coefficients: Array1::zeros(size),
        }
    }

    /// 创建单位多向量（标量部分为 1）
    pub fn one(config: AlgebraConfigRef) -> Self {
        let mut mv = Self::zeros(config.clone());
        mv.coefficients[0] = 1.0;
        mv
    }

    /// 从系数向量创建多向量
    pub fn from_coefficients(config: AlgebraConfigRef, coefficients: Vec<f64>) -> Self {
        assert_eq!(
            coefficients.len(),
            config.basis_count(),
            "Coefficient vector length {} must match algebra basis count {}",
            coefficients.len(),
            config.basis_count()
        );
        Self {
            config,
            coefficients: Array1::from_vec(coefficients),
        }
    }

    /// 从池化缓冲区创建多向量（需要 pool 特性）
    #[cfg(feature = "pool")]
    pub fn from_pooled_buffer(config: AlgebraConfigRef, buffer: Vec<f64>) -> Self {
        let expected_size = config.basis_count();
        assert_eq!(
            buffer.len(),
            expected_size,
            "Buffer size {} must match algebra basis count {}",
            buffer.len(),
            expected_size
        );
        Self {
            config,
            coefficients: Array1::from_vec(buffer),
        }
    }

    /// 获取适合池化的空缓冲区（需要 pool 特性）
    #[cfg(feature = "pool")]
    pub fn acquire_buffer(size: usize) -> Vec<f64> {
        BufferPool::get_buffer(size)
    }

    /// 归还缓冲区到池中（需要 pool 特性）
    #[cfg(feature = "pool")]
    pub fn release_buffer(buffer: Vec<f64>) {
        BufferPool::return_buffer(buffer);
    }

    /// 从标量创建多向量
    pub fn from_scalar(config: AlgebraConfigRef, scalar: f64) -> Self {
        let mut mv = Self::zeros(config);
        mv.coefficients[0] = scalar;
        mv
    }

    /// 创建基向量 e_i
    pub fn basis_vector(config: AlgebraConfigRef, i: u32) -> Self {
        assert!(
            i < config.dimension(),
            "Index {} must be less than dimension {}",
            i,
            config.dimension()
        );
        let mut mv = Self::zeros(config);
        mv.coefficients[1 << i] = 1.0;
        mv
    }

    /// 获取标量部分
    #[inline]
    pub fn scalar_part(&self) -> f64 {
        self.coefficients[0]
    }

    /// 设置标量部分
    #[inline]
    pub fn set_scalar_part(&mut self, value: f64) {
        self.coefficients[0] = value;
    }

    /// 获取第 i 个系数
    #[inline]
    pub fn get_coefficient(&self, index: usize) -> f64 {
        self.coefficients[index]
    }

    /// 设置第 i 个系数
    #[inline]
    pub fn set_coefficient(&mut self, index: usize, value: f64) {
        self.coefficients[index] = value;
    }

    /// 获取所有非零系数的索引和值
    pub fn non_zero_terms(&self) -> Vec<(usize, f64)> {
        self.coefficients
            .iter()
            .enumerate()
            .filter(|(_, &c)| c.abs() > 1e-15)
            .map(|(i, &c)| (i, c))
            .collect()
    }

    /// 检查是否为零多向量
    pub fn is_zero(&self) -> bool {
        self.coefficients.iter().all(|&c| c.abs() < 1e-15)
    }

    /// 获取阶次 r 的部分
    pub fn grade_projection(&self, r: u32) -> Self {
        let mut result = Self::zeros(self.config.clone());

        for (i, &coef) in self.coefficients.iter().enumerate() {
            if self.config.grade_of_index(i) == r {
                result.coefficients[i] = coef;
            }
        }

        result
    }

    /// 获取偶部
    pub fn even_part(&self) -> Self {
        let mut result = Self::zeros(self.config.clone());

        for (i, &coef) in self.coefficients.iter().enumerate() {
            if self.config.grade_of_index(i) % 2 == 0 {
                result.coefficients[i] = coef;
            }
        }

        result
    }

    /// 获取奇部
    pub fn odd_part(&self) -> Self {
        let mut result = Self::zeros(self.config.clone());

        for (i, &coef) in self.coefficients.iter().enumerate() {
            if self.config.grade_of_index(i) % 2 == 1 {
                result.coefficients[i] = coef;
            }
        }

        result
    }

    /// 加法
    pub fn add(&self, other: &Self) -> Self {
        assert_eq!(self.config.dimension(), other.config.dimension());
        Self {
            config: self.config.clone(),
            coefficients: &self.coefficients + &other.coefficients,
        }
    }

    /// 减法
    pub fn sub(&self, other: &Self) -> Self {
        assert_eq!(self.config.dimension(), other.config.dimension());
        Self {
            config: self.config.clone(),
            coefficients: &self.coefficients - &other.coefficients,
        }
    }

    /// 标量乘法
    pub fn scale(&self, scalar: f64) -> Self {
        Self {
            config: self.config.clone(),
            coefficients: &self.coefficients * scalar,
        }
    }

    /// 取负
    pub fn neg(&self) -> Self {
        Self {
            config: self.config.clone(),
            coefficients: -&self.coefficients,
        }
    }

    /// 获取代数配置引用
    #[inline]
    pub fn config(&self) -> &AlgebraConfig {
        &self.config
    }

    /// 获取系数切片
    #[inline]
    pub fn as_slice(&self) -> &[f64] {
        self.coefficients.as_slice().unwrap()
    }

    /// 获取可变系数切片
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        self.coefficients.as_slice_mut().unwrap()
    }

    /// 获取内部 ndarray 数组的引用
    #[inline]
    pub fn as_array(&self) -> &ndarray::Array1<f64> {
        &self.coefficients
    }
}

impl std::fmt::Display for DenseMultiVector {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let terms: Vec<String> = self
            .non_zero_terms()
            .iter()
            .map(|(i, c)| {
                let grade = self.config.grade_of_index(*i);
                let basis =
                    crate::basis::index::index_to_string(*i as u64, self.config.dimension());
                if grade == 0 {
                    format!("{:.6}", c)
                } else {
                    format!("{:.6}{}", c, basis)
                }
            })
            .collect();

        if terms.is_empty() {
            write!(f, "0")
        } else {
            write!(f, "{}", terms.join(" + "))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    fn test_config() -> AlgebraConfigRef {
        Arc::new(AlgebraConfig::euclidean(3))
    }

    #[test]
    fn test_zeros() {
        let config = test_config();
        let mv = DenseMultiVector::zeros(config.clone());
        assert!(mv.is_zero());
        assert_eq!(mv.coefficients.len(), 8);
    }

    #[test]
    fn test_one() {
        let config = test_config();
        let mv = DenseMultiVector::one(config.clone());
        assert_eq!(mv.scalar_part(), 1.0);
        assert!(mv.grade_projection(1).is_zero());
    }

    #[test]
    fn test_basis_vector() {
        let config = test_config();
        let e1 = DenseMultiVector::basis_vector(config.clone(), 0);
        assert_eq!(e1.get_coefficient(1), 1.0);
        assert_eq!(e1.config.grade_of_index(1), 1);
    }

    #[test]
    fn test_grade_projection() {
        let config = test_config();
        let mut mv = DenseMultiVector::zeros(config.clone());
        mv.set_coefficient(0, 1.0); // 标量
        mv.set_coefficient(1, 2.0); // e1
        mv.set_coefficient(3, 3.0); // e12

        let scalar_part = mv.grade_projection(0);
        assert_eq!(scalar_part.get_coefficient(0), 1.0);
        assert_eq!(
            scalar_part.coefficients.iter().skip(1).all(|&c| c == 0.0),
            true
        );

        let vector_part = mv.grade_projection(1);
        assert_eq!(vector_part.get_coefficient(1), 2.0);
    }

    #[test]
    fn test_even_odd_part() {
        let config = test_config();
        let mut mv = DenseMultiVector::zeros(config.clone());
        mv.set_coefficient(0, 1.0); // 标量 (偶)
        mv.set_coefficient(1, 2.0); // e1 (奇)
        mv.set_coefficient(3, 3.0); // e12 (偶)

        let even = mv.even_part();
        assert_eq!(even.get_coefficient(0), 1.0);
        assert_eq!(even.get_coefficient(3), 3.0);
        assert_eq!(even.get_coefficient(1), 0.0);

        let odd = mv.odd_part();
        assert_eq!(odd.get_coefficient(1), 2.0);
        assert_eq!(odd.get_coefficient(0), 0.0);
    }

    #[test]
    fn test_add_sub_scale() {
        let config = test_config();
        let a = DenseMultiVector::basis_vector(config.clone(), 0);
        let b = DenseMultiVector::basis_vector(config.clone(), 1);

        let sum = a.add(&b);
        assert_eq!(sum.get_coefficient(1), 1.0);
        assert_eq!(sum.get_coefficient(2), 1.0);

        let scaled = a.scale(2.0);
        assert_eq!(scaled.get_coefficient(1), 2.0);
    }
}
