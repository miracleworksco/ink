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
    /// Where it says why something asked has no answer: which kind of
    /// reason, if the code that refused knows (`Open`).
    pub open: Option<Open>,
}

/// Why something asked has no answer, by kind, said by the code that
/// refuses: what an answer's reader does next differs. With neither, a
/// source tried and could not decide.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Open {
    /// The design lacks something its author can add (a moving body's
    /// mass, a bound on a shape).
    Lacks,
    /// Outside what the method covers: asking again will not help (a model
    /// that leaves a body out, one value where a range is asked).
    Outside,
}

/// A refusal in words, with its kind where the refuser knows it: what a
/// solve gives where it sets nothing (`Provider::provide`). A `String` is
/// one of no stated kind.
#[derive(Clone, Debug, PartialEq)]
pub struct Refused {
    pub why: String,
    pub open: Option<Open>,
}

impl Refused {
    /// The design lacks something: `why` says what.
    pub fn lacks(why: impl Into<String>) -> Refused {
        Refused { why: why.into(), open: Some(Open::Lacks) }
    }
    /// Outside what the method covers: `why` says how.
    pub fn outside(why: impl Into<String>) -> Refused {
        Refused { why: why.into(), open: Some(Open::Outside) }
    }
}

impl From<String> for Refused {
    fn from(why: String) -> Refused {
        Refused { why, open: None }
    }
}

impl From<&str> for Refused {
    fn from(why: &str) -> Refused {
        Refused { why: why.into(), open: None }
    }
}

impl std::fmt::Display for Refused {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str(&self.why)
    }
}

impl InkError {
    pub fn msg(m: impl Into<String>) -> InkError {
        InkError { message: m.into(), file: String::new(), line: 0, col: 0, hint: None, called_from: Vec::new(), open: None }
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

/// Its words, where only words are kept.
impl From<Refused> for String {
    fn from(refused: Refused) -> String {
        refused.why
    }
}
