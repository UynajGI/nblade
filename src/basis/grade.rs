//! 阶次计算模块
//!
//! 提供阶次相关的工具函数

use crate::basis::index::{grade_of_index, BasisIndex};

/// 检查索引是否为 r-阶次
#[inline]
pub fn is_grade(index: BasisIndex, r: u32) -> bool {
    grade_of_index(index) == r
}

/// 获取指定阶次的所有基索引
///
/// # 参数
/// * `dimension` - 向量空间维度
/// * `grade` - 目标阶次
///
/// # 返回
/// 包含所有指定阶次基索引的向量
pub fn get_grade_indices(dimension: u32, grade: u32) -> Vec<BasisIndex> {
    let mut indices = Vec::new();
    let max = 1u64 << dimension;

    for i in 0..max {
        if grade_of_index(i) == grade {
            indices.push(i);
        }
    }

    indices
}

/// 计算阶二项式系数 C(n, r)
pub fn binomial_coeff(n: u32, r: u32) -> usize {
    if r > n {
        return 0;
    }

    let mut result = 1usize;
    for i in 0..r {
        result = result * (n - i) as usize / (i + 1) as usize;
    }
    result
}

/// 获取指定阶次的基数量
///
/// 等于二项式系数 C(n, r)
#[inline]
pub fn grade_count(dimension: u32, grade: u32) -> usize {
    binomial_coeff(dimension, grade)
}

/// 提取偶阶次部分的索引
pub fn get_even_grade_indices(dimension: u32) -> Vec<BasisIndex> {
    let mut indices = Vec::new();
    let max = 1u64 << dimension;

    for i in 0..max {
        if grade_of_index(i) % 2 == 0 {
            indices.push(i);
        }
    }

    indices
}

/// 提取奇阶次部分的索引
pub fn get_odd_grade_indices(dimension: u32) -> Vec<BasisIndex> {
    let mut indices = Vec::new();
    let max = 1u64 << dimension;

    for i in 0..max {
        if grade_of_index(i) % 2 == 1 {
            indices.push(i);
        }
    }

    indices
}

/// 计算阶次投影的掩码
///
/// 返回一个位掩码，其中指定阶次的位为 1
pub fn grade_mask(dimension: u32, grade: u32) -> u64 {
    let mut mask = 0u64;
    let max = 1u64 << dimension;

    for i in 0..max {
        if grade_of_index(i) == grade {
            mask |= 1u64 << i;
        }
    }

    mask
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_grade() {
        assert!(is_grade(0, 0));
        assert!(is_grade(1, 1));
        assert!(is_grade(3, 2));
        assert!(!is_grade(3, 1));
    }

    #[test]
    fn test_get_grade_indices() {
        let indices = get_grade_indices(3, 0);
        assert_eq!(indices, vec![0]);

        let indices = get_grade_indices(3, 1);
        assert_eq!(indices, vec![1, 2, 4]);

        let indices = get_grade_indices(3, 2);
        assert_eq!(indices, vec![3, 5, 6]);

        let indices = get_grade_indices(3, 3);
        assert_eq!(indices, vec![7]);
    }

    #[test]
    fn test_binomial_coeff() {
        assert_eq!(binomial_coeff(3, 0), 1);
        assert_eq!(binomial_coeff(3, 1), 3);
        assert_eq!(binomial_coeff(3, 2), 3);
        assert_eq!(binomial_coeff(3, 3), 1);
    }

    #[test]
    fn test_grade_count() {
        assert_eq!(grade_count(3, 0), 1);
        assert_eq!(grade_count(3, 1), 3);
        assert_eq!(grade_count(3, 2), 3);
        assert_eq!(grade_count(3, 3), 1);
    }

    #[test]
    fn test_even_odd_indices() {
        let even = get_even_grade_indices(3);
        assert_eq!(even, vec![0, 3, 5, 6]);

        let odd = get_odd_grade_indices(3);
        assert_eq!(odd, vec![1, 2, 4, 7]);
    }
}
