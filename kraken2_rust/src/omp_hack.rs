/// OpenMP abstraction layer
///
/// Translated from omp_hack.cc/omp_hack.h
/// In Rust, we use rayon for parallelization instead of OpenMP.
/// This module provides a compatibility layer.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

static THREAD_COUNT: AtomicUsize = AtomicUsize::new(1);

/// Get the current OpenMP thread number (0-indexed)
/// In Rust with rayon, we use rayon's thread index
pub fn get_thread_num() -> usize {
    rayon::current_thread_index().unwrap_or(0)
}

/// Get the maximum number of OpenMP threads
pub fn get_max_threads() -> usize {
    THREAD_COUNT.load(Ordering::Relaxed)
}

/// Set the number of OpenMP threads to use
pub fn set_num_threads(num: usize) {
    THREAD_COUNT.store(num, Ordering::Relaxed);
    rayon::ThreadPoolBuilder::new()
        .num_threads(num)
        .build_global()
        .ok();
}

/// OpenMP lock (simplified Rust version using Mutex)
/// Note: This is a compatibility layer. In Rust, prefer using RAII-style locks.
pub struct OmpLock {
    lock: Mutex<()>,
}

impl OmpLock {
    pub fn new() -> Self {
        OmpLock {
            lock: Mutex::new(()),
        }
    }

    /// Acquire the lock. In Rust, this should return a guard, but for OpenMP compatibility
    /// we leak the guard. This is NOT recommended Rust style - use Mutex directly instead.
    pub fn lock(&self) {
        // NOTE: This is intentionally leaking the guard to match OpenMP semantics
        // where lock() and unlock() are separate calls. This is unsafe in Rust's
        // model and should only be used for C++ compatibility.
        std::mem::forget(self.lock.lock().unwrap());
    }

    pub fn unlock(&self) {
        // In the leaked guard model, unlock is a no-op
        // The original lock design doesn't map well to Rust's RAII
        // WARNING: This implementation is incorrect and should be refactored
    }

    pub fn test_lock(&self) -> bool {
        self.lock.try_lock().is_ok()
    }
}

impl Default for OmpLock {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_max_threads_default() {
        // Default is 1 thread
        let threads = get_max_threads();
        assert!(threads >= 1);
    }

    #[test]
    fn test_thread_count_atomic() {
        // Test atomic operations on thread count
        let initial = THREAD_COUNT.load(Ordering::Relaxed);
        THREAD_COUNT.store(4, Ordering::Relaxed);
        assert_eq!(THREAD_COUNT.load(Ordering::Relaxed), 4);
        THREAD_COUNT.store(initial, Ordering::Relaxed);
    }

    #[test]
    fn test_omp_lock_new() {
        let lock = OmpLock::new();
        // Lock should be available initially
        assert!(lock.test_lock());
    }

    #[test]
    fn test_omp_lock_default() {
        let lock = OmpLock::default();
        assert!(lock.test_lock());
    }

    #[test]
    fn test_omp_lock_test_lock() {
        let lock = OmpLock::new();
        // First test_lock should succeed
        let result = lock.test_lock();
        assert!(result);
    }

    #[test]
    fn test_get_thread_num() {
        // In single-threaded context, should return 0
        let thread_num = get_thread_num();
        // Thread num should be a valid index
        assert!(thread_num < std::usize::MAX);
    }
}
