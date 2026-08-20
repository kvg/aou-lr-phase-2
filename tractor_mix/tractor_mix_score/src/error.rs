use thiserror::Error;

pub type Result<T> = std::result::Result<T, ScoreError>;

#[derive(Debug, Error)]
pub enum ScoreError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("JSON error: {0}")]
    Json(#[from] serde_json::Error),
    #[error("{0}")]
    Msg(String),
}

impl ScoreError {
    pub fn msg(s: impl Into<String>) -> Self {
        Self::Msg(s.into())
    }
}
