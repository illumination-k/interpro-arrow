use anyhow::Result;

#[derive(Debug, Clone, PartialEq, PartialOrd, Eq)]
pub enum Token {
    OpenBracket,
    CloseBracket,
    Invert,
    And,
    Or,
    Name(String),
}

impl Token {
    fn op_from_char(c: char) -> Option<Self> {
        match c {
            '(' => Some(Token::OpenBracket),
            ')' => Some(Token::CloseBracket),
            '|' | ',' => Some(Token::Or),
            '&' => Some(Token::And),
            '!' => Some(Token::Invert),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum ParseState {
    Ready,
    InName,
}

pub fn lex(s: &str) -> Result<Vec<Token>> {
    let mut state = ParseState::Ready;
    let mut tokens = vec![];

    let mut cur_name = String::new();

    for c in s.chars() {
        let op_token = Token::op_from_char(c);
        match state {
            ParseState::InName => {
                if let Some(op) = op_token {
                    tokens.push(Token::Name(cur_name.to_owned()));

                    tokens.push(op);

                    state = ParseState::Ready;
                    cur_name = String::new();
                } else if c.is_whitespace() {
                    tokens.push(Token::Name(cur_name.to_owned()));
                    state = ParseState::Ready;
                    cur_name = String::new();
                } else {
                    cur_name.push(c)
                }
            }
            ParseState::Ready => {
                if let Some(op) = op_token {
                    tokens.push(op);
                } else if !c.is_whitespace() {
                    cur_name.push(c);
                    state = ParseState::InName
                }
            }
        }
    }

    if !cur_name.is_empty() {
        tokens.push(Token::Name(cur_name.to_owned()));
    }

    Ok(tokens)
}

#[cfg(test)]
mod test_lex {
    use super::*;

    #[test]
    fn test_simple() {
        let s = "a & b";
        let tokens = lex(s).unwrap();
        assert_eq!(
            vec![
                Token::Name('a'.to_string()),
                Token::And,
                Token::Name('b'.to_string())
            ],
            tokens
        );
    }

    #[test]
    fn test_or_alias() {
        assert_eq!(lex("a | b").unwrap(), lex("a,b").unwrap());
    }

    #[test]
    fn test_with_invert() {
        let s = "!a & b";
        let tokens = lex(s).unwrap();
        assert_eq!(
            vec![
                Token::Invert,
                Token::Name("a".to_string()),
                Token::And,
                Token::Name("b".to_string())
            ],
            tokens
        )
    }

    #[test]
    fn test_with_bracket() {
        let s = "!(a & b) | c";
        let tokens = lex(s).unwrap();
        assert_eq!(
            vec![
                Token::Invert,
                Token::OpenBracket,
                Token::Name("a".to_string()),
                Token::And,
                Token::Name("b".to_string()),
                Token::CloseBracket,
                Token::Or,
                Token::Name("c".to_string()),
            ],
            tokens
        )
    }
}
