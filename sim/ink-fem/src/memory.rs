//! A solve's memory, estimated from its cells before it allocates them,
//! and a warning on stderr where the estimate is more than is available:
//! the solve goes on, but may be killed by the system for want of memory,
//! which says nothing of why.

use crate::occupancy::Classified;
use crate::surface::split_cells;
use crate::{Physics, Setup};

/// What a solve's memory grows with, counted on its cells before it
/// allocates for them: its full cells, its cells with a matrix of their
/// own (cut by the surface), and the surface its reads mesh, a cut grid
/// cell's at half its cells, or at half the finest child split within one
/// cell of it (`part_surface`), so 4^level times as many facets.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Cells {
    pub full: f64,
    pub cut: f64,
    pub surface: f64,
}

impl Cells {
    /// Classified cells (children of split grid cells among them) on a grid
    /// of `dims`.
    pub(crate) fn of(cells: &[Classified], dims: [i64; 3]) -> Cells {
        let split = split_cells(cells, dims);
        let mut cut: Vec<[i64; 3]> = cells.iter().filter(|c| c.moments.is_some()).map(|c| c.idx).collect();
        let n = cut.len() as f64;
        cut.sort_unstable();
        cut.dedup();
        let level = |i: [i64; 3]| split.binary_search_by(|c| c.0.cmp(&i)).map_or(0, |k| split[k].1);
        let surface = cut.iter().map(|&i| 4f64.powi(i32::from(level(i)))).sum();
        Cells { full: cells.len() as f64 - n, cut: n, surface }
    }
}

/// Bytes a solve holds at its peak per full cell, cut cell and unit of
/// surface (`Cells`), by what it solves, and per cell for each mode it
/// solves for: fitted to the peak resident set of 63 solves of a quarter
/// gigabyte or more (benchmarks/memory: solid cubes, plates, a ring, a
/// fin and gyroid lattices at several cell sizes, and every study of
/// benchmarks, examples and the rc car's analyses).
///
/// It is a rough estimate. Measured, it is 0.43 to 1.35 times what a solve
/// takes: blocks and plates high; real parts low (a heat sink, a compliant
/// guide, a shelf near 0.6), and lowest on coarse cells split around a
/// thin feature (a 0.3 mm fin on 2.8 mm cells, 0.43 to 0.52), where the
/// surface the reads mesh grows with more than the cells' counts.
fn per_cell(s: &Setup) -> [f64; 3] {
    match (s.physics, s.modes > 0, s.transient.is_some()) {
        (Physics::Elastic { .. }, true, _) => MODES,
        (Physics::Elastic { .. }, false, _) => ELASTIC,
        (Physics::Heat, _, true) => TRANSIENT,
        (Physics::Heat, _, false) => HEAT,
    }
}

// The surface's (the facets the reads mesh) is the same for every kind of
// solve; a cut cell's is its own matrix with what it reads by; a full
// cell's is its share of the system, the multigrid's levels and the
// fields read from them (elastic: with a buckling solve).
const ELASTIC: [f64; 3] = [3800.0, 8900.0, 4900.0];
const MODES: [f64; 3] = [2600.0, 1750.0, 4900.0];
const HEAT: [f64; 3] = [1560.0, 0.0, 4900.0];
const TRANSIENT: [f64; 3] = [1330.0, 0.0, 4900.0];
/// A mode's vectors in the eigensolve, per cell (16 modes against 4 on
/// one cube).
const PER_MODE: f64 = 500.0;

/// The bytes a solve of these cells holds at its peak, estimated.
pub(crate) fn needed(s: &Setup, c: Cells) -> f64 {
    let [f, k, a] = per_cell(s);
    f * c.full + k * c.cut + a * c.surface + PER_MODE * s.modes as f64 * (c.full + c.cut)
}

/// The bytes this process may still take: the machine's available memory
/// and swap (MemAvailable, SwapFree), or less where a control group it
/// runs in caps them (memory.max less the memory it holds that cannot be
/// reclaimed, and memory.swap.max less the swap it holds). `INK_MEMORY`
/// sets it: bytes, or a number with K, M or G (8G, 512MB, 1.5GiB); one that
/// is no size is said so on stderr and not used. None where unknown (not
/// Linux).
pub(crate) fn available() -> Option<f64> {
    if let Ok(v) = std::env::var("INK_MEMORY") {
        match bytes(&v) {
            Some(b) => return Some(b),
            None => eprintln!("ink: INK_MEMORY '{v}' is not a size (bytes, or a number with K, M or G: 8G, 512MB); it is not used"),
        }
    }
    let meminfo = std::fs::read_to_string("/proc/meminfo").ok()?;
    let kb = |key: &str| meminfo.lines().find_map(|l| l.strip_prefix(key)).and_then(|v| v.trim().trim_end_matches("kB").trim().parse::<f64>().ok()).map(|k| k * 1024.0);
    let swap = kb("SwapFree:").unwrap_or(0.0);
    let mut least = kb("MemAvailable:")? + swap;
    // cgroup v2: "0::/path"; each group from it up to the root may cap it.
    let cgroup = std::fs::read_to_string("/proc/self/cgroup").unwrap_or_default();
    if let Some(path) = cgroup.lines().find_map(|l| l.strip_prefix("0::")) {
        let mut dir = std::path::PathBuf::from("/sys/fs/cgroup").join(path.trim_start_matches('/'));
        loop {
            let read = |f: &str| std::fs::read_to_string(dir.join(f)).ok().map(|s| s.trim().to_string());
            let num = |f: &str| read(f).and_then(|v| v.parse::<f64>().ok());
            if let Some(max) = num("memory.max") {
                let anon = read("memory.stat").and_then(|st| st.lines().find_map(|l| l.strip_prefix("anon ")).and_then(|v| v.trim().parse::<f64>().ok()));
                let held = anon.or_else(|| num("memory.current")).unwrap_or(0.0);
                // "max" (or no file): its swap is not capped here.
                let swap_left = num("memory.swap.max").map_or(swap, |m| (m - num("memory.swap.current").unwrap_or(0.0)).clamp(0.0, swap));
                least = least.min((max - held).max(0.0) + swap_left);
            }
            if !dir.pop() || !dir.starts_with("/sys/fs/cgroup") || dir.as_os_str() == "/sys/fs/cgroup" {
                break;
            }
        }
    }
    Some(least)
}

/// "123", "4G", "4GB", "4 GiB", "512m", "1.5M": bytes (K, M, G, T binary).
fn bytes(v: &str) -> Option<f64> {
    let t = v.trim().to_ascii_uppercase();
    let t = t.strip_suffix("IB").or_else(|| t.strip_suffix('B')).unwrap_or(&t).trim_end();
    let (num, scale) = match t.chars().last()? {
        'K' => (&t[..t.len() - 1], 1024.0),
        'M' => (&t[..t.len() - 1], 1024.0 * 1024.0),
        'G' => (&t[..t.len() - 1], 1024.0 * 1024.0 * 1024.0),
        'T' => (&t[..t.len() - 1], 1024.0f64.powi(4)),
        _ => (t, 1.0),
    };
    num.trim().parse::<f64>().ok().filter(|x| x.is_finite() && *x >= 0.0).map(|x| x * scale)
}

/// Bytes to two significant figures, in GB, MB or kB.
fn size(b: f64) -> String {
    if b < 1024.0 {
        return format!("{b:.0} bytes");
    }
    let (v, unit) = match b {
        b if b >= 1024.0 * 1024.0 * 1024.0 => (b / (1024.0 * 1024.0 * 1024.0), "GB"),
        b if b >= 1024.0 * 1024.0 => (b / (1024.0 * 1024.0), "MB"),
        b => (b / 1024.0, "kB"),
    };
    let places = 1usize.saturating_sub(v.log10().floor() as usize);
    format!("{v:.places$} {unit}")
}

/// Says on stderr, before a solve on cells of `h` allocates for these
/// cells, where their estimate (`needed`) is more than is `available`,
/// unless `said` already; whether it has said so.
pub(crate) fn warn(s: &Setup, h: f64, c: Cells, said: bool) -> bool {
    if said {
        return true;
    }
    let need = needed(s, c);
    let Some(have) = available().filter(|&have| need > have) else { return false };
    let mm = format!("{h:.6}");
    eprintln!(
        "ink: a solve on cells of {} mm ({:.0} cells) needs about {} of memory, {} available; it may run out of memory (the estimate is rough: 0.43 to 1.35 times what solves take)",
        mm.trim_end_matches('0').trim_end_matches('.'),
        c.full + c.cut,
        size(need),
        size(have)
    );
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sizes_read_as_written() {
        let g = 1024.0 * 1024.0 * 1024.0;
        for (v, b) in [("123", 123.0), ("4G", 4.0 * g), ("8GB", 8.0 * g), ("8 GiB", 8.0 * g), ("1.5M", 1.5 * 1024.0 * 1024.0), ("512mb", 512.0 * 1024.0 * 1024.0)] {
            assert_eq!(bytes(v), Some(b), "{v}");
        }
        for v in ["x", "8GX", "", "-1G"] {
            assert_eq!(bytes(v), None, "{v}");
        }
        assert_eq!(size(378.4 * 1024.0 * 1024.0), "378 MB");
        assert_eq!(size(1.034 * 1024.0 * 1024.0 * 1024.0), "1.0 GB");
        assert_eq!(size(16.24 * 1024.0 * 1024.0 * 1024.0), "16 GB");
    }
}
