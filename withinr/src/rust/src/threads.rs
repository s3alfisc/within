//! Keep R's default CPU use bounded without changing Rayon's global pool.
use extendr_api::{error::Result, prelude::Error};
use rayon::{ThreadPool, ThreadPoolBuilder};
use std::sync::OnceLock;

pub(crate) fn run<T: Send>(f: impl FnOnce() -> T + Send) -> Result<T> {
    static POOL: OnceLock<std::result::Result<ThreadPool, String>> = OnceLock::new();
    let pool = POOL.get_or_init(|| {
        let n = match std::env::var("WITHINR_NUM_THREADS") {
            Ok(value) => value
                .parse::<usize>()
                .ok()
                .filter(|&n| n > 0)
                .ok_or("WITHINR_NUM_THREADS must be a positive integer")?,
            // Adaptive can create a second, equally sized isolation pool.
            // One worker here bounds the combined pools to CRAN's two cores.
            Err(_) => 1,
        };
        ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .map_err(|e| e.to_string())
    });
    let pool = pool.as_ref().map_err(|e| Error::Other(e.clone()))?;
    // Only native Rust work runs in the pool. R allocations and conversions
    // must remain on the calling R thread.
    Ok(pool.install(f))
}
