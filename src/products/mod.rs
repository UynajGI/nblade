//! 乘积运算模块 / Product Operations Module
//!
//! 实现几何代数的各种乘积运算，支持自适应并行调度：
//! - 低维（< 6维）：顺序执行，避免并行开销
//! - 高维（>= 6维）：并行执行，充分利用多核
//!
//! Implements geometric algebra product operations with adaptive parallel dispatch:
//! - Low-dimensional (< 6D): Sequential execution to avoid parallel overhead
//! - High-dimensional (>= 6D): Parallel execution to leverage multiple cores

pub mod geometric;
pub mod geometric_simd;
pub mod inner;
pub mod outer;

/// 并行调度阈值：系数数量阈值
/// 当 2^n >= PARALLEL_THRESHOLD 时启用并行计算
///
/// Parallel dispatch threshold: coefficient count threshold
/// Enables parallel computation when 2^n >= PARALLEL_THRESHOLD
///
/// - 64 对应 6 维代数 (2^6 = 64)
/// - 经验值：低于此值时，并行开销大于收益
pub const PARALLEL_THRESHOLD: usize = 64;

/// 检查是否应使用并行计算
/// Check if parallel computation should be used
#[inline]
pub fn should_use_parallel(basis_count: usize) -> bool {
    basis_count >= PARALLEL_THRESHOLD
}
