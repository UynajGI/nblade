//! 多向量缓冲区对象池 / Buffer pool for MultiVector allocations
//!
//! 减少密集多向量分配开销，特别适用于并行计算场景。
//! Reduces allocation overhead for dense multivectors, especially in parallel computation.
//!
//! # 设计 / Design
//!
//! - 线程本地存储，无锁访问 / Thread-local storage for lock-free access
//! - 按大小分桶存储 / Size-bucketed storage
//! - 自动限制池大小防止内存膨胀 / Automatic pool size limit to prevent memory bloat
//!
//! # 使用 / Usage
//!
//! ```ignore
//! use nblade::multivector::pool::BufferPool;
//!
//! // 获取缓冲区 / Acquire buffer
//! let mut buffer = BufferPool::get_buffer(8);
//!
//! // 使用缓冲区 / Use buffer
//! buffer[0] = 1.0;
//!
//! // 归还缓冲区 / Return buffer
//! BufferPool::return_buffer(buffer);
//! ```

use std::cell::RefCell;
use std::collections::HashMap;

/// 每个大小桶的最大缓冲区数量 / Maximum buffers per size bucket
const MAX_BUFFERS_PER_SIZE: usize = 32;

/// 支持池化的缓冲区大小（2的幂次，最高 256 = 8 维） / Poolable buffer sizes (powers of 2, max 256 = 8D)
const POOLABLE_SIZES: [usize; 7] = [4, 8, 16, 32, 64, 128, 256];

/// 检查大小是否可以被池化
fn is_poolable_size(size: usize) -> bool {
    POOLABLE_SIZES.contains(&size)
}

thread_local! {
    /// 线程本地缓冲区池 / Thread-local buffer pool
    ///
    /// 键：缓冲区大小，值：可用缓冲区列表
    /// Key: buffer size, Value: list of available buffers
    static BUFFER_POOL: RefCell<HashMap<usize, Vec<Vec<f64>>>> = RefCell::new(HashMap::new());
}

/// 缓冲区池 / Buffer Pool
///
/// 提供线程安全的缓冲区复用机制。
/// Provides thread-safe buffer reuse mechanism.
pub struct BufferPool;

impl BufferPool {
    /// 获取指定大小的缓冲区 / Acquire a buffer of the given size
    ///
    /// 优先从池中获取，池空则创建新缓冲区。
    /// Preferentially gets from pool, creates new if pool is empty.
    ///
    /// # 参数 / Arguments
    /// * `size` - 缓冲区大小 / Buffer size
    ///
    /// # 返回 / Returns
    /// 初始化为零的缓冲区 / Zero-initialized buffer
    pub fn get_buffer(size: usize) -> Vec<f64> {
        BUFFER_POOL.with(|pool| {
            let mut pool = pool.borrow_mut();

            if let Some(buffers) = pool.get_mut(&size) {
                if let Some(mut buffer) = buffers.pop() {
                    buffer.fill(0.0);
                    return buffer;
                }
            }

            vec![0.0; size]
        })
    }

    /// 归还缓冲区到池中 / Return a buffer to the pool
    ///
    /// 缓冲区会被清零并存储以供复用。
    /// Buffer is zeroed and stored for reuse.
    ///
    /// # 参数 / Arguments
    /// * `buffer` - 要归还的缓冲区 / Buffer to return
    pub fn return_buffer(mut buffer: Vec<f64>) {
        let size = buffer.len();

        if !is_poolable_size(size) {
            return;
        }

        BUFFER_POOL.with(|pool| {
            let mut pool = pool.borrow_mut();

            let buffers = pool.entry(size).or_insert_with(Vec::new);

            if buffers.len() < MAX_BUFFERS_PER_SIZE {
                buffer.fill(0.0);
                buffers.push(buffer);
            }
        });
    }

    /// 获取池统计信息 / Get pool statistics
    ///
    /// # 返回 / Returns
    /// (总缓冲区数, 总内存字节) / (total buffers, total memory bytes)
    pub fn stats() -> (usize, usize) {
        BUFFER_POOL.with(|pool| {
            let pool = pool.borrow();
            let mut total_buffers = 0;
            let mut total_bytes = 0;

            for (size, buffers) in pool.iter() {
                total_buffers += buffers.len();
                total_bytes += buffers.len() * size * std::mem::size_of::<f64>();
            }

            (total_buffers, total_bytes)
        })
    }

    /// 清空池 / Clear the pool
    ///
    /// 释放所有缓冲区，用于测试或内存压力场景。
    /// Releases all buffers, useful for testing or memory pressure scenarios.
    pub fn clear() {
        BUFFER_POOL.with(|pool| {
            pool.borrow_mut().clear();
        });
    }
}

/// 池化的缓冲区守卫 / Pooled buffer guard
///
/// RAII 守卫，自动归还缓冲区。
/// RAII guard that automatically returns buffer.
pub struct PooledBuffer {
    buffer: Option<Vec<f64>>,
}

impl PooledBuffer {
    /// 创建新的池化缓冲区 / Create new pooled buffer
    pub fn new(size: usize) -> Self {
        Self {
            buffer: Some(BufferPool::get_buffer(size)),
        }
    }

    /// 获取缓冲区引用 / Get buffer reference
    pub fn as_slice(&self) -> &[f64] {
        self.buffer.as_ref().unwrap()
    }

    /// 获取可变缓冲区引用 / Get mutable buffer reference
    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        self.buffer.as_mut().unwrap()
    }
}

impl Drop for PooledBuffer {
    fn drop(&mut self) {
        if let Some(buffer) = self.buffer.take() {
            BufferPool::return_buffer(buffer);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_and_return_buffer() {
        BufferPool::clear();

        // 获取缓冲区 / Acquire buffer
        let buffer = BufferPool::get_buffer(8);
        assert_eq!(buffer.len(), 8);
        assert!(buffer.iter().all(|&x| x == 0.0));

        // 归还缓冲区 / Return buffer
        BufferPool::return_buffer(buffer);

        // 再次获取应该复用 / Should reuse
        let buffer2 = BufferPool::get_buffer(8);
        assert_eq!(buffer2.len(), 8);

        BufferPool::return_buffer(buffer2);
    }

    #[test]
    fn test_poolable_sizes() {
        assert!(is_poolable_size(4));
        assert!(is_poolable_size(8));
        assert!(is_poolable_size(16));
        assert!(is_poolable_size(32));
        assert!(is_poolable_size(64));
        assert!(is_poolable_size(128));
        assert!(is_poolable_size(256));

        assert!(!is_poolable_size(3));
        assert!(!is_poolable_size(512));
        assert!(!is_poolable_size(1024));
    }

    #[test]
    fn test_non_poolable_size() {
        BufferPool::clear();

        let buffer = BufferPool::get_buffer(512);
        assert_eq!(buffer.len(), 512);

        BufferPool::return_buffer(buffer);

        let (buffers, _) = BufferPool::stats();
        assert_eq!(buffers, 0);

        BufferPool::clear();
    }

    #[test]
    fn test_pool_stats() {
        BufferPool::clear();

        let (buffers, bytes) = BufferPool::stats();
        assert_eq!(buffers, 0);
        assert_eq!(bytes, 0);

        // 添加一些缓冲区 / Add some buffers
        let b1 = BufferPool::get_buffer(8);
        let b2 = BufferPool::get_buffer(16);

        BufferPool::return_buffer(b1);
        BufferPool::return_buffer(b2);

        let (buffers, bytes) = BufferPool::stats();
        assert_eq!(buffers, 2);
        assert_eq!(bytes, (8 + 16) * std::mem::size_of::<f64>());

        BufferPool::clear();
    }

    #[test]
    fn test_pooled_buffer_guard() {
        BufferPool::clear();

        {
            let mut guard = PooledBuffer::new(8);
            guard.as_mut_slice()[0] = 1.0;
            assert_eq!(guard.as_slice()[0], 1.0);
        } // 自动归还 / Auto return

        let (buffers, _) = BufferPool::stats();
        assert!(buffers >= 1);

        BufferPool::clear();
    }

    #[test]
    fn test_max_pool_size() {
        BufferPool::clear();

        // 创建超过限制的缓冲区 / Create more buffers than limit
        for _ in 0..MAX_BUFFERS_PER_SIZE + 10 {
            let buffer = BufferPool::get_buffer(8);
            BufferPool::return_buffer(buffer);
        }

        // 池大小应该受限 / Pool size should be limited
        let (buffers, _) = BufferPool::stats();
        assert!(buffers <= MAX_BUFFERS_PER_SIZE);

        BufferPool::clear();
    }
}
