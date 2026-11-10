//! Keep R's default CPU use bounded without changing Rayon's global pool.
use extendr_api::{error::Result, prelude::Error};
use rayon::{ThreadPool, ThreadPoolBuilder};
use std::sync::OnceLock;

/// Worker count from `WITHINR_NUM_THREADS`, defaulting to one.
fn requested_threads() -> Result<usize> {
    match std::env::var("WITHINR_NUM_THREADS") {
        Ok(value) => value
            .trim()
            .parse::<usize>()
            .ok()
            .filter(|&n| n > 0)
            .ok_or_else(|| {
                Error::Other(format!(
                    "WITHINR_NUM_THREADS must be a positive integer, got \"{value}\""
                ))
            }),
        // Adaptive can create a second, equally sized isolation pool.
        // One worker here bounds the combined pools to CRAN's two cores.
        Err(_) => Ok(1),
    }
}

/// Run native work on the private pool, built on first use. A failed build
/// (such as an invalid `WITHINR_NUM_THREADS`) is not cached, so correcting
/// the variable takes effect on the next call.
pub(crate) fn run<T: Send>(f: impl FnOnce() -> T + Send) -> Result<T> {
    static POOL: OnceLock<ThreadPool> = OnceLock::new();
    let pool = match POOL.get() {
        Some(pool) => pool,
        None => {
            let pool = ThreadPoolBuilder::new()
                .num_threads(requested_threads()?)
                .build()
                .map_err(|e| Error::Other(format!("thread pool build failed: {e}")))?;
            // Calls arrive on the single R thread, so no other build can race this one.
            POOL.get_or_init(|| pool)
        }
    };
    // Only native Rust work runs in the pool. R allocations and conversions
    // must remain on the calling R thread.
    Ok(pool.install(f))
}
