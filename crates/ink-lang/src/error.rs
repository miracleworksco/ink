//! The language's error: a message, where in which file, a hint, and the
//! calls that led there.

#[derive(Clone, Debug)]
pub struct InkError {
    pub message: String,
    pub file: String,
    pub line: u32,
    pub col: u32,
    pub hint: Option<String>,
    /// For errors inside functions (std or the user's): the calls that led
    /// there, innermost first, as (function, file, line, col) of each call.
    pub called_from: Vec<(String, String, u32, u32)>,
}

impl InkError {
    pub fn msg(m: impl Into<String>) -> InkError {
        InkError { message: m.into(), file: String::new(), line: 0, col: 0, hint: None, called_from: Vec::new() }
    }

    pub fn to_json(&self) -> serde_json::Value {
        let mut v = serde_json::json!({ "message": self.message });
        if self.line > 0 {
            v["at"] = serde_json::json!({ "file": self.file, "line": self.line, "col": self.col });
        }
        if let Some(h) = &self.hint {
            v["hint"] = serde_json::json!(h);
        }
        if !self.called_from.is_empty() {
            // A frame repeated in a row (a recursion) is listed once, with
            // how many times.
            let mut frames: Vec<(serde_json::Value, usize)> = Vec::new();
            for (f, file, line, col) in &self.called_from {
                let j = serde_json::json!({"fn": f, "file": file, "line": line, "col": col});
                match frames.last_mut() {
                    Some((last, n)) if *last == j => *n += 1,
                    _ => frames.push((j, 1)),
                }
            }
            v["called_from"] = frames
                .into_iter()
                .map(|(mut j, n)| {
                    if n > 1 {
                        j["times"] = serde_json::json!(n);
                    }
                    j
                })
                .collect();
        }
        v
    }
}

impl From<String> for InkError {
    fn from(m: String) -> Self {
        InkError::msg(m)
    }
}
