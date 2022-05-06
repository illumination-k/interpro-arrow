/*
This parser is referred to https://github.com/Smittyvb/ttw/blob/f77fa34e62739b0225847317d243fc1a4ab29b96/taglogic/src/bool.rs#L187
*/

pub mod ast;
pub mod expr;
pub mod lex;

pub use expr::{Expr, MAX_RECURSION};
