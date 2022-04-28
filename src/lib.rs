use serde::Deserialize;
use serde::Serialize;

fn sqrt(n: f64) -> f64 {
    f64::sqrt(n)
}



#[derive(Debug, Deserialize, Serialize)]
#[serde(rename_all = "PascalCase")]
struct Matrix {
    size: usize,
    matrix: Vec<f64>
}

impl Matrix {

    fn new(s: usize) -> Matrix {
        let mut mat: Vec<f64> = Vec::new();

        for _x in 0..(s*s) {
			mat.push(0.0);
		}

		Matrix {
			matrix: mat,
			size: s
		}
    }

    fn set(&mut self, row: usize, col:usize, val:f64) {
		self.matrix[((row-1) * self.size + (col - 1)) as usize] = val
	}

    fn get(&mut self, row: usize, col:usize) -> f64 {
        self.matrix[((row-1) * self.size + (col - 1))]
    }

    fn compl(&mut self, row: usize, col:usize) -> Matrix {
		let mut mat: Matrix = Matrix::new(self.size - 1);

		let mut row_1: usize;
		let mut col_1: usize;

		for x in 1..=mat.size {
			for y in 1..=mat.size {
				row_1 = x;
				col_1 = y;

				if x >= row {
					row_1 += 1;
				}

				if y >= col {
					col_1 += 1;
				}

				mat.set(x, y, self.get(row_1, col_1))
			}
		}

		mat
	}

    // Laplacian expansion algorithm (recursive)
	pub fn deter(&mut self) ->f64 {
		let mut deter: f64 = 0.0;

		// basic case
		if self.size == 2 {
			deter += self.get(1,1) * self.get(2,2) - self.get(1,2) * self.get(2,1)
		}

		// recursion
		else {
			for x in 1..=self.size {
                match (x+1)%2 {
                    0 => deter += self.get(1,x) * self.compl(1,x).deter(),
                    1 => deter += (-1.00) * self.get(1,x) * self.compl(1,x).deter(),
                    _ => break
                }
			}
		}

		deter

	}
}

fn determinant(m: String) -> f64 {
    let mut mat: Matrix = serde_json::from_str(&m).expect("The JSON was corrupted");
    return mat.deter()
}



#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn float_sqrt_test() {
        assert_eq!(16.0, sqrt(256.0));
    }

    #[test]
    fn determinant_test_1() {
        let mut test_object = Matrix::new(2);
        test_object.set(1, 1, 1.0);
        test_object.set(1, 2, 2.0);
        test_object.set(2, 1, 3.0);
        test_object.set(2, 2, 4.0);

        let string_result = serde_json::to_string(&test_object);
        let string = string_result.expect("Err");
        
        assert_eq!(determinant(string), -2.0);
    }

    #[test]
    fn determinant_test_2() {
        let mut test_object = Matrix::new(3);
        test_object.set(1, 1, 1.0);
        test_object.set(1, 2, 2.0);
        test_object.set(1, 3, 3.0);
        test_object.set(2, 1, 4.0);
        test_object.set(2, 2, 5.0);
        test_object.set(2, 3, 6.0);
        test_object.set(3, 1, 7.0);
        test_object.set(3, 2, 8.0);
        test_object.set(3, 3, 9.0);

        let string_result = serde_json::to_string(&test_object);
        let string = string_result.expect("Err");
        
        assert_eq!(determinant(string), 0.0);
    }
}