use std::{collections::VecDeque, error::Error};

use super::ast::{Node, ParseError};
use super::lex::{lex, Token};

#[derive(Debug, Clone, PartialEq, Eq)]
enum ExprData {
    Empty,
    HasNodes(Node),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Expr(ExprData);

pub const MAX_RECURSION: u16 = 20;

impl Expr {
    pub fn from_string(s: &str) -> Result<Self, Box<dyn Error>> {
        // lex and convert to a deque
        let mut tokens: VecDeque<Token> = VecDeque::from(lex(s)?);
        if tokens.is_empty() {
            // no tokens
            return Ok(Self(ExprData::Empty));
        }

        let ast = Node::munch_tokens(&mut tokens, MAX_RECURSION)?;
        if !tokens.is_empty() {
            return Err(Box::new(ParseError::new(
                "expected EOF, found extra tokens",
            )));
        }

        Ok(Self(ExprData::HasNodes(ast)))
    }

    pub fn matches(&self, tags: &[&str]) -> Result<bool, Box<dyn Error>> {
        match &self.0 {
            ExprData::Empty => Ok(true),
            ExprData::HasNodes(node) => node.matches(tags),
        }
    }

    // pub fn matches_domains(&self, gene_record: &GeneRecord) -> Result<bool, Box<dyn Error>> {
    //     let tags: Vec<&str> = gene_record
    //         .iter_domains()
    //         .map(|domain| domain.domain_name.as_str())
    //         .collect();

    //     self.matches(&tags)
    // }
}

#[cfg(test)]
mod test_expr {
    use super::*;
    #[test]
    fn or_alias() {
        assert_eq!(
            Expr::from_string("a | b").unwrap().0,
            Expr::from_string("a , b").unwrap().0,
        );
    }

    #[test]
    fn simple_and() {
        let expr = Expr::from_string("a & b").unwrap();
        assert_eq!(
            expr.to_owned().0,
            ExprData::HasNodes(Node::And {
                lhs: Box::new(Node::Name("a".to_string())),
                rhs: Box::new(Node::Name("b".to_string())),
            })
        );

        assert!(expr.matches(&["a", "b"]).unwrap());
        assert!(!expr.matches(&["a"]).unwrap());
        assert!(!expr.matches(&["c"]).unwrap());
    }

    #[test]
    fn test_simple_count_and() {
        let expr = Expr::from_string("a$2 & b").unwrap();
        assert!(expr.matches(&["a", "a", "b"]).unwrap());
    }

    #[test]
    fn simple_inversion() {
        let expr = Expr::from_string("!a & b").unwrap();
        assert_eq!(
            expr.to_owned().0,
            ExprData::HasNodes(Node::And {
                lhs: Box::new(Node::Invert(Box::new(Node::Name("a".to_string())))),
                rhs: Box::new(Node::Name("b".to_string())),
            })
        );

        assert!(expr.matches(&["b"]).unwrap());
        assert!(!expr.matches(&["a", "b"]).unwrap());
        assert!(!expr.matches(&["a"]).unwrap());
        assert!(!expr.matches(&["c"]).unwrap());
    }

    #[test]
    fn simple_and_matching() {
        assert!(Expr::from_string("a & b & c")
            .unwrap()
            .matches(&["a", "b", "c"])
            .unwrap());
        assert!(!Expr::from_string("a & b & c")
            .unwrap()
            .matches(&["a", "c"])
            .unwrap());
        assert!(!Expr::from_string("a & b & c")
            .unwrap()
            .matches(&["a", "b"])
            .unwrap());
        assert!(!Expr::from_string("a & b & c")
            .unwrap()
            .matches(&["c", "b"])
            .unwrap());
        assert!(Expr::from_string("a & b & c")
            .unwrap()
            .matches(&["a", "b", "c", "d"])
            .unwrap());
        assert!(!Expr::from_string("a & b & c")
            .unwrap()
            .matches(&["a", "c", "d"])
            .unwrap());
        assert!(!Expr::from_string("a & b & c")
            .unwrap()
            .matches(&["a", "b", "d"])
            .unwrap());
        assert!(!Expr::from_string("a & b & c")
            .unwrap()
            .matches(&["c", "b", "d"])
            .unwrap());
    }

    #[test]
    fn simple_or_matching() {
        assert!(Expr::from_string("a | b | c")
            .unwrap()
            .matches(&["a", "b", "c"])
            .unwrap());
        assert!(Expr::from_string("a | b | c")
            .unwrap()
            .matches(&["a", "c"])
            .unwrap());
        assert!(Expr::from_string("a | b | c")
            .unwrap()
            .matches(&["a", "b"])
            .unwrap());
        assert!(Expr::from_string("a | b | c")
            .unwrap()
            .matches(&["c", "b"])
            .unwrap());
        assert!(Expr::from_string("a | b | c")
            .unwrap()
            .matches(&["c"])
            .unwrap());
        assert!(Expr::from_string("a | b | c")
            .unwrap()
            .matches(&["b"])
            .unwrap());
        assert!(Expr::from_string("a | b | c")
            .unwrap()
            .matches(&["a"])
            .unwrap());
        assert!(Expr::from_string("a | b | c")
            .unwrap()
            .matches(&["a", "b", "c", "d"])
            .unwrap());
        assert!(Expr::from_string("a | b | c")
            .unwrap()
            .matches(&["a", "c", "d"])
            .unwrap());
        assert!(Expr::from_string("a | b | c")
            .unwrap()
            .matches(&["a", "b", "d"])
            .unwrap());
        assert!(Expr::from_string("a | b | c")
            .unwrap()
            .matches(&["c", "b", "d"])
            .unwrap());
        assert!(Expr::from_string("a | b | c")
            .unwrap()
            .matches(&["c", "d"])
            .unwrap());
        assert!(Expr::from_string("a | b | c")
            .unwrap()
            .matches(&["b", "d"])
            .unwrap());
        assert!(Expr::from_string("a | b | c")
            .unwrap()
            .matches(&["a", "d"])
            .unwrap());
        assert!(!Expr::from_string("a | b | c")
            .unwrap()
            .matches(&["ddwf"])
            .unwrap());
        assert!(!Expr::from_string("a | b | c")
            .unwrap()
            .matches(&["d"])
            .unwrap());
        assert!(!Expr::from_string("a | b | c")
            .unwrap()
            .matches(&["hdwf", "dtw"])
            .unwrap());
    }

    #[test]
    fn lone_name() {
        assert!(Expr::from_string("a").unwrap().matches(&["a"]).unwrap());
        assert!(Expr::from_string("a")
            .unwrap()
            .matches(&["a", "b"])
            .unwrap());
        assert!(!Expr::from_string("a").unwrap().matches(&["b"]).unwrap());
    }

    #[test]
    fn lone_inverted_name() {
        assert!(!Expr::from_string("!a").unwrap().matches(&["a"]).unwrap());
        assert!(!Expr::from_string("!a")
            .unwrap()
            .matches(&["a", "b"])
            .unwrap());
        assert!(Expr::from_string("!a").unwrap().matches(&["b"]).unwrap());
    }

    #[test]
    fn lone_inverted_bracketed_name() {
        assert!(!Expr::from_string("!(a)").unwrap().matches(&["a"]).unwrap());
        assert!(!Expr::from_string("!(a)")
            .unwrap()
            .matches(&["a", "b"])
            .unwrap());
        assert!(Expr::from_string("!(a)").unwrap().matches(&["b"]).unwrap());
    }

    #[test]
    fn check_expr() {
        let s = "a | b | c";
        let expr = Expr::from_string(s).unwrap();
        assert!(expr.matches(&["a"]).unwrap());
        let s = "(a | b | c) | (d & e & c)";
        let expr = Expr::from_string(s).unwrap();
        assert!(expr.matches(&["a"]).unwrap());
        assert!(expr.matches(&["d", "e", "c"]).unwrap());
        assert!(!expr.matches(&["d"]).unwrap());
    }
}
