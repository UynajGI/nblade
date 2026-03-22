//! 基索引模块
//!
//! 实现几何代数基向量的位运算表示和运算
//!
//! 对于 n 维空间，每个基向量用一个 n 位整数表示：
//! - 位模式：第 i 位为 1 表示包含 ei
//! - 例如 n=5: 0b00101 = e1∧e3 (第 0 和第 2 位为 1)

use crate::signature::Signature;

/// 基索引类型（支持最多 64 维空间）
pub type BasisIndex = u64;

/// 计算基索引的阶次（popcount = 1 的位数）
///
/// # 参数
/// * `index` - 基索引
///
/// # 返回
/// 阶次（1 的位数）
///
/// # 示例
/// ```
/// use nblade::basis::grade_of_index;
///
/// assert_eq!(grade_of_index(0), 0);    // 标量
/// assert_eq!(grade_of_index(1), 1);    // e1
/// assert_eq!(grade_of_index(3), 2);    // e1∧e2
/// assert_eq!(grade_of_index(7), 3);    // e1∧e2∧e3
/// ```
#[inline]
pub fn grade_of_index(index: BasisIndex) -> u32 {
    index.count_ones()
}

/// 计算两个基索引的几何积
///
/// # 参数
/// * `i` - 第一个基索引
/// * `j` - 第二个基索引
/// * `signature` - 度量签名
///
/// # 返回
/// (结果索引，符号，度量因子)
///
/// # 算法
/// 1. 结果索引 = i XOR j（对称差）
/// 2. 符号由交换次数决定
/// 3. 度量因子由公共位的基向量平方值决定
#[inline]
pub fn basis_geometric_product(
    i: BasisIndex,
    j: BasisIndex,
    signature: &Signature,
) -> (BasisIndex, i8, f64) {
    // 结果索引 = XOR（对称差）
    let result = i ^ j;

    // 计算符号（总是需要计算重排序符号）
    let sign = compute_geometric_sign(i, j);

    // 计算度量因子（处理非欧几里得签名）
    let common = i & j;
    let metric_factor = if common == 0 {
        1.0
    } else {
        compute_metric_factor(common, signature)
    };

    (result, sign, metric_factor)
}

/// 计算几何积的符号
///
/// 统计将 i 和 j 合并需要的交换次数
/// 符号 = (-1)^{pairs (i_elem, j_elem) where i_elem > j_elem}
fn compute_geometric_sign(i: BasisIndex, j: BasisIndex) -> i8 {
    // 优化算法：统计 j 中每个设置位对应的 i 中更高位的数量
    // 这等同于计算需要交换的次数
    let mut count = 0;
    let mut j_bits = j;

    while j_bits != 0 {
        let bit = j_bits.trailing_zeros();
        // 统计 i 中高于 bit 位置的位数（i > j 的对数）
        // 使用 !0 << (bit + 1) 来获取所有高于 bit 的位
        let mask = !0u64 << (bit + 1);
        count += (i & mask).count_ones();
        j_bits &= !(1u64 << bit);
    }

    if count % 2 == 0 {
        1
    } else {
        -1
    }
}

/// 计算度量因子
fn compute_metric_factor(common: BasisIndex, signature: &Signature) -> f64 {
    let mut factor = 1.0;
    let mut bits = common;

    while bits != 0 {
        let bit = bits.trailing_zeros();
        factor *= signature.basis_square_sign(bit);
        if factor == 0.0 {
            return 0.0;
        }
        bits &= !(1u64 << bit);
    }

    factor
}

/// 计算外积（只有当没有公共位时才非零）
///
/// # 返回
/// Some((结果索引，符号)) 或 None（如果有公共位）
#[inline]
pub fn basis_outer_product(i: BasisIndex, j: BasisIndex) -> Option<(BasisIndex, i8)> {
    if i & j == 0 {
        let result = i | j;
        let sign = compute_outer_sign(i, j);
        Some((result, sign))
    } else {
        None
    }
}

/// 计算外积的符号
/// 符号 = (-1)^{pairs (i_elem, j_elem) where i_elem > j_elem}
fn compute_outer_sign(i: BasisIndex, j: BasisIndex) -> i8 {
    // 外积符号与几何积相同：统计 i > j 的对数
    if j == 0 {
        return 1;
    }

    let mut count = 0;
    let mut j_bits = j;

    while j_bits != 0 {
        let bit = j_bits.trailing_zeros();
        let mask = !0u64 << (bit + 1);
        count += (i & mask).count_ones();
        j_bits &= !(1u64 << bit);
    }

    if count % 2 == 0 {
        1
    } else {
        -1
    }
}

/// 计算左内积的基运算
///
/// A⌋B 的基索引运算：只有当 i 是 j 的子集时才非零
/// 特殊情况：标量 ⌋ 任何元素 = 该元素（标量作为收缩的单位元）
#[inline]
pub fn basis_left_inner(i: BasisIndex, j: BasisIndex) -> Option<(BasisIndex, i8)> {
    // 标量 ⌋ 任何元素 = 该元素
    if i == 0 {
        return Some((j, 1));
    }
    // 左内积条件：i 必须是 j 的子集
    if i & j == i {
        let result = i ^ j; // j - i
        let sign = compute_left_inner_sign(i, j);
        Some((result, sign))
    } else {
        None
    }
}

/// 计算左内积的符号
///
/// 对于基向量 e_I ⌋ e_J，当 I ⊆ J 时，结果为 ±e_{J-I}
/// 符号由将 I 中元素移到 J 前面所需的交换次数决定
///
/// # 算法
/// 统计将 i 的位移动到 j 的最低位之前需要的交换次数
fn compute_left_inner_sign(i: BasisIndex, j: BasisIndex) -> i8 {
    let mut sign = 1i8;
    let mut i_bits = i;

    while i_bits != 0 {
        let bit_i = i_bits.trailing_zeros();
        // 统计 j 中低于 bit_i 位置的位数
        let mask = (1u64 << bit_i) - 1;
        let swaps = (j & mask).count_ones();
        if swaps % 2 == 1 {
            sign = -sign;
        }
        i_bits &= !(1u64 << bit_i);
    }

    sign
}

/// 计算右内积的基运算
///
/// A⌊B 的基索引运算：只有当 j 是 i 的子集时才非零
#[inline]
pub fn basis_right_inner(i: BasisIndex, j: BasisIndex) -> Option<(BasisIndex, i8)> {
    // 右内积条件：j 必须是 i 的子集
    if i & j == j {
        let result = i ^ j; // i - j
        let sign = compute_right_inner_sign(i, j);
        Some((result, sign))
    } else {
        None
    }
}

/// 计算右内积的符号
fn compute_right_inner_sign(i: BasisIndex, j: BasisIndex) -> i8 {
    let _grade_j = grade_of_index(j);
    let mut count = 0;
    let mut j_bits = j;

    while j_bits != 0 {
        let bit = j_bits.trailing_zeros();
        let mask = (1u64 << bit) - 1;
        count += (i & mask).count_ones();
        j_bits &= !(1u64 << bit);
    }

    if count % 2 == 0 {
        1
    } else {
        -1
    }
}

/// 获取基索引的字符串表示
pub fn index_to_string(index: BasisIndex, dimension: u32) -> String {
    if index == 0 {
        return "1".to_string(); // 标量
    }

    let mut parts = Vec::new();
    for i in 0..dimension {
        if index & (1u64 << i) != 0 {
            parts.push(format!("e{}", i + 1));
        }
    }

    parts.join("∧")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_grade_of_index() {
        assert_eq!(grade_of_index(0), 0);
        assert_eq!(grade_of_index(1), 1);
        assert_eq!(grade_of_index(3), 2);
        assert_eq!(grade_of_index(7), 3);
        assert_eq!(grade_of_index(15), 4);
    }

    #[test]
    fn test_geometric_product_euclidean() {
        let sig = Signature::euclidean(3);

        // e1 * e1 = 1
        let (result, sign, metric) = basis_geometric_product(1, 1, &sig);
        assert_eq!(result, 0);
        assert_eq!(sign, 1);
        assert_eq!(metric, 1.0);

        // e1 * e2 = e12
        let (result, sign, metric) = basis_geometric_product(1, 2, &sig);
        assert_eq!(result, 3);
        assert_eq!(sign, 1);
        assert_eq!(metric, 1.0);

        // e2 * e1 = -e12
        let (result, sign, metric) = basis_geometric_product(2, 1, &sig);
        assert_eq!(result, 3);
        assert_eq!(sign, -1);
        assert_eq!(metric, 1.0);
    }

    #[test]
    fn test_geometric_product_spacetime() {
        let sig = Signature::new(1, 1, 0);

        // e0 * e0 = 1
        let (result, sign, metric) = basis_geometric_product(1, 1, &sig);
        assert_eq!(result, 0);
        assert_eq!(sign, 1);
        assert_eq!(metric, 1.0);

        // e1 * e1 = -1
        let (result, sign, metric) = basis_geometric_product(2, 2, &sig);
        assert_eq!(result, 0);
        assert_eq!(sign, 1);
        assert_eq!(metric, -1.0);
    }

    #[test]
    fn test_outer_product() {
        // e1 ∧ e2 = e12
        let result = basis_outer_product(1, 2);
        assert!(result.is_some());
        let (idx, sign) = result.unwrap();
        assert_eq!(idx, 3);
        assert_eq!(sign, 1);

        // e1 ∧ e1 = 0 (有公共位)
        let result = basis_outer_product(1, 1);
        assert!(result.is_none());
    }

    #[test]
    fn test_left_inner() {
        // e1 ⌋ (e1∧e2) = e2
        let result = basis_left_inner(1, 3);
        assert!(result.is_some());
        let (idx, sign) = result.unwrap();
        assert_eq!(idx, 2);
        assert_eq!(sign, 1);

        // e2 ⌋ (e1∧e2) = -e1
        let result = basis_left_inner(2, 3);
        assert!(result.is_some());
        let (idx, sign) = result.unwrap();
        assert_eq!(idx, 1);
        assert_eq!(sign, -1);
    }

    #[test]
    fn test_index_to_string() {
        assert_eq!(index_to_string(0, 3), "1");
        assert_eq!(index_to_string(1, 3), "e1");
        assert_eq!(index_to_string(3, 3), "e1∧e2");
        assert_eq!(index_to_string(7, 3), "e1∧e2∧e3");
    }
}
