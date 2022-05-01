use serde::Deserialize;
use serde::Serialize;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub fn sqrt(n: f64) -> f64 {
    f64::sqrt(n)
}

trait MatrixTraits {
    fn set(&mut self, row: usize, col:usize, val:f64);
    fn get(&self, row: usize, col:usize) -> f64;
    fn transpose(&self) -> Self;
}

#[derive(Debug, Deserialize, Serialize, PartialEq)]
#[serde(rename_all = "PascalCase")]
struct SquareMatrix {
    size: usize,
    matrix: Vec<f64>
}

impl SquareMatrix {

    pub fn new(s: usize) -> SquareMatrix {
        let mut mat: Vec<f64> = Vec::new();

        for _x in 0..(s*s) {
			mat.push(0.0);
		}

		SquareMatrix {
			matrix: mat,
			size: s
		}
    }

    fn compl(&self, row: usize, col:usize) -> f64 {
		let mut mat: SquareMatrix = SquareMatrix::new(self.size - 1);

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

        let mut deter = 0.0;

        match (row + col) % 2 {
            0 => deter = mat.deter(),
            1 => deter = (-1.00) * mat.deter(),
            _ => ()
        }

		deter
	}

    // Laplacian expansion algorithm (recursive)
	pub fn deter(&self) ->f64 {
		let mut deter: f64 = 0.0;

		// basic case
		if self.size == 2 {
			deter += self.get(1,1) * self.get(2,2) - self.get(1,2) * self.get(2,1)
		}

		// recursion
		else {
			for x in 1..=self.size {
                deter += self.get(1, x) * self.compl(1, x)
			}
		}

		deter

	}

    // Matrix inversion
    pub fn invert(self) -> SquareMatrix {
        let d = 1.0 / self.deter();
        let transposed_mat = self.transpose();
        let mut inv_mat = SquareMatrix::new(self.size);

        for x in 1..=self.size {
            for y in 1..=self.size {
                inv_mat.set(x, y, d * transposed_mat.compl(x, y))
            }
        }

        inv_mat
    }
}

impl MatrixTraits for SquareMatrix {
    fn set(&mut self, row: usize, col:usize, val:f64) {
		self.matrix[((row-1) * self.size + (col - 1)) as usize] = val
	}

    fn get(&self, row: usize, col:usize) -> f64 {
        self.matrix[((row-1) * self.size + (col - 1))]
    }

    fn transpose(&self) -> SquareMatrix {
        let mut transposed_mat = SquareMatrix::new(self.size);
        for x in 1..=self.size {
            for y in 1..=self.size {
                transposed_mat.set(x, y, self.get(y, x));
            }
        }

        transposed_mat
    }
}

#[wasm_bindgen]
pub fn determinant(m: String) -> f64 {
    let mat: SquareMatrix = serde_json::from_str(&m).expect("The JSON was corrupted");
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
        let mut test_object = SquareMatrix::new(2);
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
        let mut test_object = SquareMatrix::new(3);
        test_object.set(1, 1, 1.0);
        test_object.set(1, 2, -1.0);
        test_object.set(1, 3, 0.0);
        test_object.set(2, 1, 1.0);
        test_object.set(2, 2, 0.0);
        test_object.set(2, 3, -1.0);
        test_object.set(3, 1, 2.0);
        test_object.set(3, 2, 3.0);
        test_object.set(3, 3, -4.0);

        let string_result = serde_json::to_string(&test_object);
        let string = string_result.expect("Err");
        
        assert_eq!(determinant(string), 1.0);
    }

    #[test]
    fn transposition_test() {
        let mut original_matrix = SquareMatrix::new(3);
        let mut transposed_matrix = SquareMatrix::new(3);

        original_matrix.set(1, 1, 1.0);
        original_matrix.set(1, 2, -1.0);
        original_matrix.set(1, 3, 0.0);
        original_matrix.set(2, 1, 1.0);
        original_matrix.set(2, 2, 0.0);
        original_matrix.set(2, 3, -1.0);
        original_matrix.set(3, 1, 2.0);
        original_matrix.set(3, 2, 3.0);
        original_matrix.set(3, 3, -4.0);

        transposed_matrix.set(1, 1, 1.0);
        transposed_matrix.set(1, 2, 1.0);
        transposed_matrix.set(1, 3, 2.0);
        transposed_matrix.set(2, 1, -1.0);
        transposed_matrix.set(2, 2, 0.0);
        transposed_matrix.set(2, 3, 3.0);
        transposed_matrix.set(3, 1, 0.0);
        transposed_matrix.set(3, 2, -1.0);
        transposed_matrix.set(3, 3, -4.0);

        assert_eq!(original_matrix.transpose(), transposed_matrix);
    }

    #[test]
    fn inversion_test() {
        let mut original_matrix = SquareMatrix::new(3);
        let mut inverse_matrix = SquareMatrix::new(3);

        original_matrix.set(1, 1, 1.0);
        original_matrix.set(1, 2, -1.0);
        original_matrix.set(1, 3, 0.0);
        original_matrix.set(2, 1, 1.0);
        original_matrix.set(2, 2, 0.0);
        original_matrix.set(2, 3, -1.0);
        original_matrix.set(3, 1, 2.0);
        original_matrix.set(3, 2, 3.0);
        original_matrix.set(3, 3, -4.0);

        inverse_matrix.set(1, 1, 3.0);
        inverse_matrix.set(1, 2, -4.0);
        inverse_matrix.set(1, 3, 1.0);
        inverse_matrix.set(2, 1, 2.0);
        inverse_matrix.set(2, 2, -4.0);
        inverse_matrix.set(2, 3, 1.0);
        inverse_matrix.set(3, 1, 3.0);
        inverse_matrix.set(3, 2, -5.0);
        inverse_matrix.set(3, 3, 1.0);

        assert_eq!(original_matrix, inverse_matrix.invert());
    }
}