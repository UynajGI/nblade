//! 预计算乘法表模块
//!
//! 使用并行计算预先生成几何积、内积、外积的查找表
//! 加速后续的多向量运算

use crate::algebra_config::AlgebraConfig;
use crate::basis::index::{
    basis_geometric_product, basis_left_inner, basis_outer_product, basis_right_inner, BasisIndex,
};
use rayon::prelude::*;

/// 乘法表条目
#[derive(Debug, Clone, Copy)]
pub struct MultiplicationEntry {
    /// 结果基索引
    pub result_index: BasisIndex,
    /// 符号 (+1 或 -1)
    pub sign: i8,
    /// 度量因子（处理非欧几里得签名）
    pub metric_factor: f64,
}

impl MultiplicationEntry {
    /// 创建零条目
    pub fn zero() -> Self {
        Self {
            result_index: 0,
            sign: 1,
            metric_factor: 0.0,
        }
    }

    /// 检查是否为零条目
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.metric_factor == 0.0
    }
}

/// 几何积乘法表
#[derive(Debug, Clone)]
pub struct GeometricProductTable {
    /// 向量空间维度
    dimension: u32,
    /// 表大小 (2^n × 2^n)
    size: usize,
    /// 乘法表数据
    table: Vec<MultiplicationEntry>,
}

impl GeometricProductTable {
    /// 创建新的乘法表（并行计算）
    pub fn new(config: &AlgebraConfig) -> Self {
        let n = config.dimension() as usize;
        let size = (1 << n) * (1 << n);
        let signature = &config.signature;

        // 并行计算乘法表
        let table: Vec<MultiplicationEntry> = (0..size)
            .into_par_iter()
            .map(|idx| {
                let i = (idx / (1 << n)) as BasisIndex;
                let j = (idx % (1 << n)) as BasisIndex;
                let (result, sign, metric) = basis_geometric_product(i, j, signature);
                MultiplicationEntry {
                    result_index: result,
                    sign,
                    metric_factor: metric,
                }
            })
            .collect();

        Self {
            dimension: config.dimension(),
            size,
            table,
        }
    }

    /// 获取乘法表条目
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> &MultiplicationEntry {
        let idx = i * (1 << self.dimension) + j;
        &self.table[idx]
    }

    /// 获取表大小
    #[inline]
    pub fn size(&self) -> usize {
        self.size
    }

    /// 获取维度
    #[inline]
    pub fn dimension(&self) -> u32 {
        self.dimension
    }
}

/// 外积乘法表
#[derive(Debug, Clone)]
pub struct OuterProductTable {
    dimension: u32,
    #[allow(dead_code)]
    size: usize,
    table: Vec<Option<MultiplicationEntry>>,
}

impl OuterProductTable {
    /// 创建外积表
    pub fn new(config: &AlgebraConfig) -> Self {
        let n = config.dimension() as usize;
        let size = (1 << n) * (1 << n);

        let table: Vec<Option<MultiplicationEntry>> = (0..size)
            .into_par_iter()
            .map(|idx| {
                let i = (idx / (1 << n)) as BasisIndex;
                let j = (idx % (1 << n)) as BasisIndex;
                basis_outer_product(i, j).map(|(result, sign)| MultiplicationEntry {
                    result_index: result,
                    sign,
                    metric_factor: 1.0,
                })
            })
            .collect();

        Self {
            dimension: config.dimension(),
            size,
            table,
        }
    }

    /// 获取外积条目
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> Option<&MultiplicationEntry> {
        let idx = i * (1 << self.dimension) + j;
        self.table[idx].as_ref()
    }
}

/// 左内积乘法表
#[derive(Debug, Clone)]
pub struct LeftInnerProductTable {
    dimension: u32,
    #[allow(dead_code)]
    size: usize,
    table: Vec<Option<MultiplicationEntry>>,
}

impl LeftInnerProductTable {
    /// 创建左内积表
    pub fn new(config: &AlgebraConfig) -> Self {
        let n = config.dimension() as usize;
        let size = (1 << n) * (1 << n);

        let table: Vec<Option<MultiplicationEntry>> = (0..size)
            .into_par_iter()
            .map(|idx| {
                let i = (idx / (1 << n)) as BasisIndex;
                let j = (idx % (1 << n)) as BasisIndex;
                basis_left_inner(i, j).map(|(result, sign)| MultiplicationEntry {
                    result_index: result,
                    sign,
                    metric_factor: 1.0,
                })
            })
            .collect();

        Self {
            dimension: config.dimension(),
            size,
            table,
        }
    }

    /// 获取左内积条目
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> Option<&MultiplicationEntry> {
        let idx = i * (1 << self.dimension) + j;
        self.table[idx].as_ref()
    }
}

/// 右内积乘法表
#[derive(Debug, Clone)]
pub struct RightInnerProductTable {
    dimension: u32,
    #[allow(dead_code)]
    size: usize,
    table: Vec<Option<MultiplicationEntry>>,
}

impl RightInnerProductTable {
    /// 创建右内积表
    pub fn new(config: &AlgebraConfig) -> Self {
        let n = config.dimension() as usize;
        let size = (1 << n) * (1 << n);

        let table: Vec<Option<MultiplicationEntry>> = (0..size)
            .into_par_iter()
            .map(|idx| {
                let i = (idx / (1 << n)) as BasisIndex;
                let j = (idx % (1 << n)) as BasisIndex;
                basis_right_inner(i, j).map(|(result, sign)| MultiplicationEntry {
                    result_index: result,
                    sign,
                    metric_factor: 1.0,
                })
            })
            .collect();

        Self {
            dimension: config.dimension(),
            size,
            table,
        }
    }

    /// 获取右内积条目
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> Option<&MultiplicationEntry> {
        let idx = i * (1 << self.dimension) + j;
        self.table[idx].as_ref()
    }
}

/// 完整的乘法表集合
#[derive(Debug, Clone)]
pub struct MultiplicationTables {
    pub geometric: GeometricProductTable,
    pub outer: OuterProductTable,
    pub left_inner: LeftInnerProductTable,
    pub right_inner: RightInnerProductTable,
}

impl MultiplicationTables {
    /// 创建所有乘法表
    pub fn new(config: &AlgebraConfig) -> Self {
        Self {
            geometric: GeometricProductTable::new(config),
            outer: OuterProductTable::new(config),
            left_inner: LeftInnerProductTable::new(config),
            right_inner: RightInnerProductTable::new(config),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra_config::AlgebraConfig;
    use crate::signature::Signature;

    #[test]
    fn test_geometric_table_euclidean() {
        let config = AlgebraConfig::euclidean(2);
        let table = GeometricProductTable::new(&config);

        // e1 * e1 = 1
        let entry = table.get(1, 1);
        assert_eq!(entry.result_index, 0);
        assert_eq!(entry.sign, 1);
        assert_eq!(entry.metric_factor, 1.0);

        // e1 * e2 = e12
        let entry = table.get(1, 2);
        assert_eq!(entry.result_index, 3);
        assert_eq!(entry.sign, 1);

        // e2 * e1 = -e12
        let entry = table.get(2, 1);
        assert_eq!(entry.result_index, 3);
        assert_eq!(entry.sign, -1);
    }

    #[test]
    fn test_geometric_table_spacetime() {
        let config = AlgebraConfig::new(2, Signature::new(1, 1, 0));
        let table = GeometricProductTable::new(&config);

        // e0 * e0 = 1
        let entry = table.get(1, 1);
        assert_eq!(entry.result_index, 0);
        assert_eq!(entry.metric_factor, 1.0);

        // e1 * e1 = -1
        let entry = table.get(2, 2);
        assert_eq!(entry.result_index, 0);
        assert_eq!(entry.metric_factor, -1.0);
    }

    #[test]
    fn test_outer_table() {
        let config = AlgebraConfig::euclidean(2);
        let table = OuterProductTable::new(&config);

        // e1 ∧ e2 = e12
        let entry = table.get(1, 2);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 3);

        // e1 ∧ e1 = 0
        let entry = table.get(1, 1);
        assert!(entry.is_none());
    }

    #[test]
    fn test_left_inner_table() {
        let config = AlgebraConfig::euclidean(2);
        let table = LeftInnerProductTable::new(&config);

        // e1 ⌋ (e1∧e2) = e2
        let entry = table.get(1, 3);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 2);
    }

    #[test]
    fn test_right_inner_table() {
        let config = AlgebraConfig::euclidean(2);
        let table = RightInnerProductTable::new(&config);

        let entry = table.get(3, 1);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 2);

        let entry = table.get(3, 2);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 1);

        let entry = table.get(1, 1);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 0);
    }

    #[test]
    fn test_multiplication_entry_zero() {
        let zero_entry = MultiplicationEntry::zero();
        assert_eq!(zero_entry.result_index, 0);
        assert_eq!(zero_entry.sign, 1);
        assert_eq!(zero_entry.metric_factor, 0.0);
        assert!(zero_entry.is_zero());
    }

    #[test]
    fn test_multiplication_entry_is_zero() {
        let zero_entry = MultiplicationEntry::zero();
        assert!(zero_entry.is_zero());

        let non_zero = MultiplicationEntry {
            result_index: 1,
            sign: 1,
            metric_factor: 1.0,
        };
        assert!(!non_zero.is_zero());

        let negative_metric = MultiplicationEntry {
            result_index: 0,
            sign: -1,
            metric_factor: -1.0,
        };
        assert!(!negative_metric.is_zero());
    }

    #[test]
    fn test_geometric_table_size_and_dimension() {
        let config_2d = AlgebraConfig::euclidean(2);
        let table_2d = GeometricProductTable::new(&config_2d);
        assert_eq!(table_2d.dimension(), 2);
        assert_eq!(table_2d.size(), 16);

        let config_3d = AlgebraConfig::euclidean(3);
        let table_3d = GeometricProductTable::new(&config_3d);
        assert_eq!(table_3d.dimension(), 3);
        assert_eq!(table_3d.size(), 64);

        let config_4d = AlgebraConfig::euclidean(4);
        let table_4d = GeometricProductTable::new(&config_4d);
        assert_eq!(table_4d.dimension(), 4);
        assert_eq!(table_4d.size(), 256);
    }

    #[test]
    fn test_geometric_table_3d() {
        let config = AlgebraConfig::euclidean(3);
        let table = GeometricProductTable::new(&config);

        let entry = table.get(1, 2);
        assert_eq!(entry.result_index, 3);
        assert_eq!(entry.sign, 1);

        let entry = table.get(2, 1);
        assert_eq!(entry.result_index, 3);
        assert_eq!(entry.sign, -1);

        let entry = table.get(1, 4);
        assert_eq!(entry.result_index, 5);
        assert_eq!(entry.sign, 1);

        let entry = table.get(2, 4);
        assert_eq!(entry.result_index, 6);
        assert_eq!(entry.sign, 1);

        let entry = table.get(3, 4);
        assert_eq!(entry.result_index, 7);
        assert_eq!(entry.sign, 1);

        let entry = table.get(4, 3);
        assert_eq!(entry.result_index, 7);
    }

    #[test]
    fn test_outer_table_3d() {
        let config = AlgebraConfig::euclidean(3);
        let table = OuterProductTable::new(&config);

        let entry = table.get(1, 2);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 3);

        let entry = table.get(2, 1);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 3);
        assert_eq!(entry.unwrap().sign, -1);

        let entry = table.get(3, 4);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 7);

        let entry = table.get(1, 6);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 7);

        let entry = table.get(1, 1);
        assert!(entry.is_none());
    }

    #[test]
    fn test_left_inner_table_3d() {
        let config = AlgebraConfig::euclidean(3);
        let table = LeftInnerProductTable::new(&config);

        let entry = table.get(1, 3);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 2);

        let entry = table.get(1, 7);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 6);

        let entry = table.get(2, 7);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 5);
        assert_eq!(entry.unwrap().sign, -1);

        let entry = table.get(4, 7);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 3);
    }

    #[test]
    fn test_right_inner_table_3d() {
        let config = AlgebraConfig::euclidean(3);
        let table = RightInnerProductTable::new(&config);

        let entry = table.get(3, 1);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 2);

        let entry = table.get(7, 1);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 6);

        let entry = table.get(7, 2);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 5);
        assert_eq!(entry.unwrap().sign, -1);

        let entry = table.get(7, 4);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 3);
    }

    #[test]
    fn test_multiplication_tables() {
        let config = AlgebraConfig::euclidean(2);
        let tables = MultiplicationTables::new(&config);

        let geo_entry = tables.geometric.get(1, 1);
        assert_eq!(geo_entry.result_index, 0);

        let outer_entry = tables.outer.get(1, 2);
        assert!(outer_entry.is_some());
        assert_eq!(outer_entry.unwrap().result_index, 3);

        let left_entry = tables.left_inner.get(1, 3);
        assert!(left_entry.is_some());
        assert_eq!(left_entry.unwrap().result_index, 2);

        let right_entry = tables.right_inner.get(3, 1);
        assert!(right_entry.is_some());
        assert_eq!(right_entry.unwrap().result_index, 2);
    }

    #[test]
    fn test_outer_table_edge_cases() {
        let config = AlgebraConfig::euclidean(3);
        let table = OuterProductTable::new(&config);

        let entry = table.get(0, 1);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 1);
        assert_eq!(entry.unwrap().sign, 1);

        let entry = table.get(2, 0);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 2);
        assert_eq!(entry.unwrap().sign, 1);

        let entry = table.get(0, 0);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().result_index, 0);
    }

    #[test]
    fn test_geometric_table_minkowski() {
        let config = AlgebraConfig::new(2, Signature::new(1, 1, 0));
        let table = GeometricProductTable::new(&config);

        let entry = table.get(1, 1);
        assert_eq!(entry.result_index, 0);
        assert_eq!(entry.sign, 1);
        assert_eq!(entry.metric_factor, 1.0);

        let entry = table.get(2, 2);
        assert_eq!(entry.result_index, 0);
        assert_eq!(entry.sign, 1);
        assert_eq!(entry.metric_factor, -1.0);

        let entry = table.get(1, 2);
        assert_eq!(entry.result_index, 3);
        assert_eq!(entry.sign, 1);

        let entry = table.get(2, 1);
        assert_eq!(entry.result_index, 3);
        assert_eq!(entry.sign, -1);
    }

    #[test]
    fn test_multiplication_entry_debug_clone_copy() {
        let entry = MultiplicationEntry {
            result_index: 5,
            sign: -1,
            metric_factor: 2.0,
        };

        let cloned = entry.clone();
        assert_eq!(cloned.result_index, entry.result_index);
        assert_eq!(cloned.sign, entry.sign);
        assert_eq!(cloned.metric_factor, entry.metric_factor);

        let copied = entry;
        assert_eq!(copied.result_index, entry.result_index);

        let debug_str = format!("{:?}", entry);
        assert!(debug_str.contains("MultiplicationEntry"));
    }
}
