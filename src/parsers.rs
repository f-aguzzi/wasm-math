use std::collections::HashMap;
use crate::{SquareMatrix, Matrix};

#[derive(Clone, Copy)]
enum SupportedTypes {
    SquareMatrix,
    Matrix,
    Number,
    Constant
}

struct Data {
    square_matrices: HashMap<String, SquareMatrix>,
    matrices: HashMap<String, Matrix>,
    numbers: HashMap<String, f64>,
    constants: HashMap<String, f64>,
    all_variables: HashMap<String, SupportedTypes>
}

impl Data {
    pub fn new() -> Data {
        Data {
            square_matrices: HashMap::new(),
            matrices: HashMap::new(),
            numbers: HashMap::new(),
            constants: HashMap::new(),
            all_variables: HashMap::new(),
        }
    }

    pub fn add_matrix(&mut self, m: Matrix, name: String) {

        if self.all_variables.contains_key(&name) {
            let var_type: SupportedTypes = *self.all_variables.get(&name).unwrap();

            match var_type {
                SupportedTypes::SquareMatrix => {
                    self.square_matrices.remove(&name);
                }
                SupportedTypes::Matrix => {
                    self.matrices.remove(&name);
                }
                SupportedTypes::Number => {
                    self.numbers.remove(&name);
                }
                SupportedTypes::Constant => ()
            }
        }

        self.all_variables.insert(name.to_owned(), SupportedTypes::Matrix);
        self.matrices.insert(name, m);
    }
}

fn parse_value(value: String, name: String) {

}

fn parse(text: String) {
    let mut tokens = text.split_whitespace().peekable();

    while tokens.peek().is_some() {
        let cur: String = tokens.next().unwrap().to_owned();
        if cur == "let" && tokens.peek().is_some() {
            let name: String = tokens.next().unwrap().to_owned();

            if name != "=" && tokens.peek().is_some() && tokens.peek().unwrap() == &"=" {
                tokens.next();
                if tokens.peek().is_some() {
                    let value = tokens.next().unwrap().to_owned();
                    parse_value(value, name);
                }
            } else {
                panic!("Syntax Error");
            }
        }

    }
}