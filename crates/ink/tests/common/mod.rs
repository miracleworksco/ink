//! What the integration tests share.
// A file whose tests all need a package this build lacks uses none of it.
#![allow(dead_code)]

use std::path::{Path, PathBuf};

/// A directory of a test's own under the system's temporary directory
/// (`ink-TAG-PID`), empty when made and removed with what is in it when
/// dropped, a test that panics too: tests leave nothing behind.
pub struct Scratch(PathBuf);

impl Scratch {
    pub fn new(tag: &str) -> Scratch {
        let dir = std::env::temp_dir().join(format!("ink-{tag}-{}", std::process::id()));
        _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        Scratch(dir)
    }
}

impl std::ops::Deref for Scratch {
    type Target = Path;
    fn deref(&self) -> &Path {
        &self.0
    }
}

impl AsRef<Path> for Scratch {
    fn as_ref(&self) -> &Path {
        &self.0
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        _ = std::fs::remove_dir_all(&self.0);
    }
}
