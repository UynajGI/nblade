//! 度量签名模块
//! 
//! 定义几何代数的度量签名 (p, q, r)，其中：
//! - p: 正平方基向量数量 (e_i² = +1)
//! - q: 负平方基向量数量 (e_i² = -1)
//! - r: 零平方基向量数量 (e_i² = 0)

use std::fmt;

/// 几何代数的度量签名
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Signature {
    /// 正平方基向量数量 (e_i² = +1)
    pub positive: u32,
    /// 负平方基向量数量 (e_i² = -1)
    pub negative: u32,
    /// 零平方基向量数量 (e_i² = 0)
    pub null: u32,
}

impl Signature {
    /// 创建新的签名
    /// 
    /// # 参数
    /// * `p` - 正平方基向量数量
    /// * `q` - 负平方基向量数量
    /// * `r` - 零平方基向量数量
    /// 
    /// # 示例
    /// ```
    /// use nblade::signature::Signature;
    /// 
    /// // 欧几里得空间 G(3,0,0)
    /// let euclidean = Signature::new(3, 0, 0);
    /// 
    /// // 时空代数 G(1,3,0)
    /// let spacetime = Signature::new(1, 3, 0);
    /// 
    /// // 共形几何代数 G(4,1,0)
    /// let cga = Signature::new(4, 1, 0);
    /// ```
    pub fn new(p: u32, q: u32, r: u32) -> Self {
        Self {
            positive: p,
            negative: q,
            null: r,
        }
    }

    /// 创建欧几里得签名 (n, 0, 0)
    pub fn euclidean(n: u32) -> Self {
        Self::new(n, 0, 0)
    }

    /// 创建时空签名 (1, n-1, 0)
    pub fn spacetime(n: u32) -> Self {
        Self::new(1, n - 1, 0)
    }

    /// 获取向量空间维度 n = p + q + r
    #[inline]
    pub fn dimension(&self) -> u32 {
        self.positive + self.negative + self.null
    }

    /// 获取第 i 个基向量的平方值
    /// 
    /// # 参数
    /// * `i` - 基向量索引 (0-based)
    /// 
    /// # 返回
    /// * `1.0` - 如果 i < p (正平方)
    /// * `-1.0` - 如果 p <= i < p+q (负平方)
    /// * `0.0` - 如果 p+q <= i (零平方)
    #[inline]
    pub fn basis_square(&self, i: u32) -> f64 {
        if i >= self.dimension() {
            return 1.0;
        }
        if i < self.positive {
            1.0
        } else if i < self.positive + self.negative {
            -1.0
        } else {
            0.0
        }
    }

    /// 获取第 i 个基向量的平方符号（用于位运算优化）
    #[inline]
    pub fn basis_square_sign(&self, i: u32) -> f64 {
        self.basis_square(i)
    }

    /// 检查是否为欧几里得签名
    pub fn is_euclidean(&self) -> bool {
        self.negative == 0 && self.null == 0
    }

    /// 检查是否为非退化签名（无零向量）
    pub fn is_non_degenerate(&self) -> bool {
        self.null == 0
    }
}

impl fmt::Display for Signature {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "G({},{},{})", self.positive, self.negative, self.null)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_signature_creation() {
        let sig = Signature::new(3, 1, 0);
        assert_eq!(sig.positive, 3);
        assert_eq!(sig.negative, 1);
        assert_eq!(sig.null, 0);
        assert_eq!(sig.dimension(), 4);
    }

    #[test]
    fn test_euclidean_signature() {
        let sig = Signature::euclidean(3);
        assert_eq!(sig.positive, 3);
        assert_eq!(sig.negative, 0);
        assert_eq!(sig.null, 0);
        assert!(sig.is_euclidean());
    }

    #[test]
    fn test_spacetime_signature() {
        let sig = Signature::spacetime(4);
        assert_eq!(sig.positive, 1);
        assert_eq!(sig.negative, 3);
        assert_eq!(sig.null, 0);
    }

    #[test]
    fn test_basis_square() {
        let sig = Signature::new(2, 1, 1);
        assert_eq!(sig.basis_square(0), 1.0);   // e0² = +1
        assert_eq!(sig.basis_square(1), 1.0);   // e1² = +1
        assert_eq!(sig.basis_square(2), -1.0);  // e2² = -1
        assert_eq!(sig.basis_square(3), 0.0);   // e3² = 0
    }

    #[test]
    fn test_display() {
        let sig = Signature::new(3, 1, 0);
        assert_eq!(format!("{}", sig), "G(3,1,0)");
    }
}
