//! Progress, for a caller that may stop a long run before its answer.

/// Results as they are decided, for a caller that may stop a long run before
/// its answer (`ink mcp` at its time limit reports what finished): one JSON
/// line each, appended to the file the environment's `INK_PROGRESS` names,
/// if set. Lines: `{"total": n, "of": "asserts", "kept": true}` first (`of`
/// the answer's list the rows go in; `kept`: each is kept in the cache as it
/// is decided, so a run again takes it from there), `{"running": name}` as
/// one starts (several may be running at once: interfere's moving pairs; a
/// claim's `ranged`: the ranged params it reads, each with its span), `{"sweep": {...}}` as steps over
/// some of them start and about each second they go on (motion's
/// `step_over`: how many params it steps, the boxes checked of the most it
/// may make, the share of their box proven), `{"done": name, "status": s, "ms": t,
/// "row": r}` as one is decided (`row`: its row of the answer, where the
/// caller gives one), and `{"answer": {...}}`, the answer's other
/// fields (its units, scope, method). `ink interfere` reports its pairs and
/// a package's command its items so (an optimizer its designs). Each line
/// says when it was written (`at`, Unix seconds): a caller that stops the
/// run knows how long the one running had run.
pub fn progress(mut line: serde_json::Value) {
    use std::io::Write;
    static SINK: std::sync::OnceLock<Option<std::sync::Mutex<std::fs::File>>> = std::sync::OnceLock::new();
    let sink = SINK.get_or_init(|| {
        let path = std::env::var_os("INK_PROGRESS")?;
        std::fs::OpenOptions::new().create(true).append(true).open(path).ok().map(std::sync::Mutex::new)
    });
    if let Some(f) = sink
        && let Ok(mut f) = f.lock()
    {
        if let (Some(o), Ok(now)) = (line.as_object_mut(), std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)) {
            o.insert("at".into(), serde_json::json!(now.as_secs_f64()));
        }
        // One write, so a reader (or a kill) never meets half a line:
        // `writeln!` on a File writes each formatted piece on its own.
        let _ = f.write_all(format!("{line}\n").as_bytes());
    }
}
