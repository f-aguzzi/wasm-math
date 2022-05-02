use serde::Deserialize;
use serde::Serialize;
use wasm_bindgen::prelude::*;
use std::io::{Error, ErrorKind};

mod statistics;

#[wasm_bindgen]
pub fn sqrt(n: f64) -> f64 {
    f64::sqrt(n)
}

trait MatrixTraits {
    type SuperMatrix;

    fn set(&mut self, row: usize, col:usize, val:f64);
    fn get(&self, row: usize, col:usize) -> f64;
    fn transpose(&self) -> Self;
    fn sum(&self, mat2: Self) -> Result<Self::SuperMatrix, Error>;
}

#[derive(Debug, Deserialize, Serialize, PartialEq, Clone)]
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
    pub fn invert(self) -> Result<SquareMatrix, Error> {
        let d = self.deter();

        if d == 0.0 {
            return Err(Error::new(ErrorKind::Other, "Determinant is 0: can't be inverted"));
        } 
            
        let transposed_mat = self.transpose();
        let mut inv_mat = SquareMatrix::new(self.size);

        for x in 1..=self.size {
            for y in 1..=self.size {
                inv_mat.set(x, y, d * transposed_mat.compl(x, y))
            }
        }

        Ok(inv_mat)
        
    }
}

impl MatrixTraits for SquareMatrix {

    type SuperMatrix = SquareMatrix;

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

    fn sum(&self, mat2: SquareMatrix) -> Result<SquareMatrix, std::io::Error> {
        
        let mut mat;

        if self.size == mat2.size {
            mat = SquareMatrix::new(self.size);
            for x in 1..=self.size {
                for y in 1..=self.size {
                    mat.set(x, y, self.get(x,y) + mat2.get(x,y))
                }
            }

            return Ok(mat);
        } else {
            return Err(Error::new(ErrorKind::InvalidInput, "Mismatched matrix dimensions"));
        }


    }
}

#[wasm_bindgen]
pub fn determinant(m: String) -> f64 {
    let mat: SquareMatrix = serde_json::from_str(&m).expect("The JSON was corrupted");
    return mat.deter()
}


#[derive(Debug, Deserialize, Serialize, PartialEq, Clone)]
#[serde(rename_all = "PascalCase")]
struct Matrix {
    sizex: usize,
    sizey: usize,
    matrix: Vec<f64>
}

impl Matrix {
    pub fn new(sizex: usize, sizey: usize) -> Matrix {
        let mut mat: Vec<f64> = Vec::new();

        for _x in 0..(sizex * sizey) {
			mat.push(0.0);
		}

        Matrix {
            sizex: sizex,
            sizey: sizey,
            matrix: mat
        }
    }
}

impl MatrixTraits for Matrix {

    type SuperMatrix = Matrix;

    fn set(&mut self, row: usize, col:usize, val:f64) {
		self.matrix[((row-1) * self.sizex + (col - 1)) as usize] = val
	}

    fn get(&self, row: usize, col:usize) -> f64 {
        self.matrix[((row-1) * self.sizex + (col - 1))]
    }

    fn transpose(&self) -> Matrix {
        let mut transposed_mat = Matrix::new(self.sizey, self.sizex);
        for x in 1..=self.sizex {
            for y in 1..=self.sizey {
                transposed_mat.set(x, y, self.get(y, x));
            }
        }

        transposed_mat
    }

    fn sum(&self, mat2: Matrix) -> Result<Matrix, Error> {

        if self.sizex == mat2.sizex && self.sizey == mat2.sizey {
            let mut mat = Matrix::new(self.sizex, self.sizey);

            for x in 1..=self.sizey {
                for y in 1..=self.sizex {
                    println!("Riga {}, colonna {}", x, y);
                    mat.set(x, y, self.get(x,y) + mat2.get(x,y))
                }
            }

            return Ok(mat);

        } else {
            return Err(Error::new(ErrorKind::InvalidInput, "Mismatched matrix dimensions"));
        }

        
    }
}


#[cfg(test)]
mod tests {

    use super::*;

    // Function tests

    #[test]
    fn float_sqrt_test() {
        assert_eq!(16.0, sqrt(256.0));
    }

    // SquareMatrix tests

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

        // Reference matrix
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

        // Check if the determinant is computed correctly
        assert_eq!(test_object.deter(), 1.0);

        // Testing the passing of a matrix object as a stringified JSON
        let string = serde_json::to_string(&test_object).expect("Err");
        let recomposed_matrix = serde_json::from_str(&string).expect("Err");
        assert_eq!(test_object, recomposed_matrix);
        
        // Check if the determinant of the stringified matrix is computed correctly
        assert_eq!(determinant(string), 1.0);

        // Ckeck if the determinant of the recomposed matrix is correct
        assert_eq!(recomposed_matrix.deter(), 1.0);
    }

    #[test]
    fn square_transposition_test() {

        // Reference matrix
        let mut original_matrix = SquareMatrix::new(3);
        original_matrix.set(1, 1, 1.0);
        original_matrix.set(1, 2, -1.0);
        original_matrix.set(1, 3, 0.0);
        original_matrix.set(2, 1, 1.0);
        original_matrix.set(2, 2, 0.0);
        original_matrix.set(2, 3, -1.0);
        original_matrix.set(3, 1, 2.0);
        original_matrix.set(3, 2, 3.0);
        original_matrix.set(3, 3, -4.0);

        // Hand-transposed matrix
        let mut transposed_matrix = SquareMatrix::new(3);
        transposed_matrix.set(1, 1, 1.0);
        transposed_matrix.set(1, 2, 1.0);
        transposed_matrix.set(1, 3, 2.0);
        transposed_matrix.set(2, 1, -1.0);
        transposed_matrix.set(2, 2, 0.0);
        transposed_matrix.set(2, 3, 3.0);
        transposed_matrix.set(3, 1, 0.0);
        transposed_matrix.set(3, 2, -1.0);
        transposed_matrix.set(3, 3, -4.0);

        // Check if the function-transposed and the hand-transposed matrix are equal
        assert_eq!(original_matrix.transpose(), transposed_matrix);
    }

    #[test]
    fn inversion_test() {
        let mut original_matrix = SquareMatrix::new(3);

        original_matrix.set(1, 1, 1.0);
        original_matrix.set(1, 2, -1.0);
        original_matrix.set(1, 3, 0.0);
        original_matrix.set(2, 1, 1.0);
        original_matrix.set(2, 2, 0.0);
        original_matrix.set(2, 3, -1.0);
        original_matrix.set(3, 1, 2.0);
        original_matrix.set(3, 2, 3.0);
        original_matrix.set(3, 3, -4.0);

        let mut inverse_matrix = SquareMatrix::new(3);

        inverse_matrix.set(1, 1, 3.0);
        inverse_matrix.set(1, 2, -4.0);
        inverse_matrix.set(1, 3, 1.0);
        inverse_matrix.set(2, 1, 2.0);
        inverse_matrix.set(2, 2, -4.0);
        inverse_matrix.set(2, 3, 1.0);
        inverse_matrix.set(3, 1, 3.0);
        inverse_matrix.set(3, 2, -5.0);
        inverse_matrix.set(3, 3, 1.0);

        // Check if the matrix gets inverted correctly
        assert_eq!(original_matrix, inverse_matrix.invert().unwrap());

        let mut det_0_matrix = SquareMatrix::new(3);
        
        det_0_matrix.set(1, 1, 1.0);
        det_0_matrix.set(1, 2, 2.0);
        det_0_matrix.set(2, 1, 2.0);
        det_0_matrix.set(2, 2, 4.0);

        // Check for error on non-invertible matrices
        assert_eq!(det_0_matrix.invert().is_err(), true);
    }

    
    #[test]
    fn square_sum_test() {
        let mut original_matrix = SquareMatrix::new(3);

        original_matrix.set(1, 1, 1.0);
        original_matrix.set(1, 2, -1.0);
        original_matrix.set(1, 3, 0.0);
        original_matrix.set(2, 1, 1.0);
        original_matrix.set(2, 2, 0.0);
        original_matrix.set(2, 3, -1.0);
        original_matrix.set(3, 1, 2.0);
        original_matrix.set(3, 2, 3.0);
        original_matrix.set(3, 3, -4.0);

        let mut second_matrix = SquareMatrix::new(2);

        second_matrix.set(1, 1, 1.0);
        second_matrix.set(1, 2, -1.0);
        second_matrix.set(2, 1, 1.0);
        second_matrix.set(2, 2, 0.0);

        assert_eq!(original_matrix.sum(second_matrix).is_err(), true);
        assert_eq!(original_matrix.sum(original_matrix.to_owned()).is_err(), false);

        let mut third_matrix = SquareMatrix::new(3);

        third_matrix.set(1, 1, 2.0);
        third_matrix.set(1, 2, -2.0);
        third_matrix.set(1, 3, 0.0);
        third_matrix.set(2, 1, 2.0);
        third_matrix.set(2, 2, 0.0);
        third_matrix.set(2, 3, -2.0);
        third_matrix.set(3, 1, 4.0);
        third_matrix.set(3, 2, 6.0);
        third_matrix.set(3, 3, -8.0);

        assert_eq!(original_matrix.sum(original_matrix.to_owned()).unwrap(), third_matrix);
    }


    // Matrix tests

    #[test]
    fn transposition_test() {

        // Reference matrix
        let mut original_matrix = Matrix::new(3,2);
        original_matrix.set(1, 1, 1.0);
        original_matrix.set(1, 2, -1.0);
        original_matrix.set(1, 3, 0.0);
        original_matrix.set(2, 1, 1.0);
        original_matrix.set(2, 2, 0.0);
        original_matrix.set(2, 3, -1.0);

        // Matrix transposed by hand
        let mut transpose_matrix = Matrix::new(2,3);
        transpose_matrix.set(1, 1, 1.0);
        transpose_matrix.set(1, 2, 1.0);
        transpose_matrix.set(2, 1, -1.0);
        transpose_matrix.set(2, 2, 0.0);
        transpose_matrix.set(3, 1, 0.0);
        transpose_matrix.set(3, 2, -1.0);

        // Check if the function-transposed and the hand-transposed matrix are equal
        assert_eq!(original_matrix.transpose(), transpose_matrix);
    }

    #[test]
    fn sum_test() {
        let mut original_matrix = Matrix::new(3,2);

        original_matrix.set(1, 1, 1.0);
        original_matrix.set(1, 2, -1.0);
        original_matrix.set(1, 3, 0.0);
        original_matrix.set(2, 1, 1.0);
        original_matrix.set(2, 2, 0.0);
        original_matrix.set(2, 3, -1.0);

        let mut second_matrix = Matrix::new(2,3);
        
        second_matrix.set(1, 1, 1.0);
        second_matrix.set(1, 2, 1.0);
        second_matrix.set(2, 1, -1.0);
        second_matrix.set(2, 2, 0.0);
        second_matrix.set(3, 1, 0.0);
        second_matrix.set(3, 2, -1.0);

        // Check whether non-invertible matrices are detected correctly
        assert_eq!(original_matrix.sum(second_matrix).is_err(), true);
        assert_eq!(original_matrix.sum(original_matrix.to_owned()).is_err(), false);

        let mut third_matrix = Matrix::new(3,2);

        third_matrix.set(1, 1, 2.0);
        third_matrix.set(1, 2, -2.0);
        third_matrix.set(1, 3, 0.0);
        third_matrix.set(2, 1, 2.0);
        third_matrix.set(2, 2, 0.0);
        third_matrix.set(2, 3, -2.0);

        // Check if the matrix gets inverted correctly
        assert_eq!(original_matrix.sum(original_matrix.to_owned()).unwrap(), third_matrix);
    }
}