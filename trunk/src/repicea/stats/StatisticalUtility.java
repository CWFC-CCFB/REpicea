/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2012 Mathieu Fortin for Rouge-Epicea
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * This library is distributed with the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * Please see the license at http://www.gnu.org/copyleft/lesser.html.
 */
package repicea.stats;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import repicea.math.Matrix;

/**
 * This class contains static methods that are useful for statistical regressions.
 * @author Mathieu Fortin - August 2012
 */
public class StatisticalUtility {

	private static Random random;

	public static enum TypeMatrixR {LINEAR, LINEAR_LOG, COMPOUND_SYMMETRY, POWER, ARMA}
	
	
	/**
	 * Construct a within-subject correlation matrix using a variance parameter, a correlation parameter and a column vector of coordinates.
	 * @param coordinates a column vector of coordinates from which the distances are calculated
	 * @param varianceParameter the variance parameter
	 * @param covarianceParameter the covariance parameter
	 * @param type the type of correlation
	 * @return the resulting matrix
	 */
	public static Matrix constructRMatrix(Matrix coordinates, double varianceParameter, double covarianceParameter, TypeMatrixR type) {
		if (!coordinates.isColumnVector()) {
			throw new UnsupportedOperationException("Matrix.constructRMatrix() : The coordinates matrix is not a column vector");
		} else {
			int nrow = coordinates.m_iRows;
			Matrix matrixR = new Matrix(nrow,nrow);
			for (int i = 0; i < nrow; i++) {
				for (int j = i; j < nrow; j++) {
					double corr = 0d;
					switch(type) {
					case LINEAR:					// linear case
						corr = 1 - covarianceParameter * Math.abs(coordinates.m_afData[i][0] - coordinates.m_afData[j][0]);
						if (corr >= 0) {
							matrixR.m_afData[i][j] = varianceParameter * corr;
							matrixR.m_afData[j][i] = varianceParameter * corr;
						}
						break;
					case LINEAR_LOG:				// linear log case
						if (Math.abs(coordinates.m_afData[i][0] - coordinates.m_afData[j][0]) == 0) {
							corr = 1d;
						} else {
							corr = 1 - covarianceParameter*Math.log(Math.abs(coordinates.m_afData[i][0] - coordinates.m_afData[j][0]));
						}
						if (corr >= 0) {
							matrixR.m_afData[i][j] = varianceParameter * corr;
							matrixR.m_afData[j][i] = varianceParameter * corr;
						}
						break;
					case COMPOUND_SYMMETRY:
						if (i == j) {
							matrixR.m_afData[i][j] = varianceParameter + covarianceParameter;
						} else {
							matrixR.m_afData[i][j] = covarianceParameter;
							matrixR.m_afData[j][i] = covarianceParameter;
						}
						break;
					case POWER:                  // power case
                        if (Math.abs(coordinates.m_afData[i][0] - coordinates.m_afData[j][0]) == 0) {
                              corr = 1d;
                        } else {
                              corr = Math.pow (covarianceParameter, (Math.abs (coordinates.m_afData[i][0] - coordinates.m_afData[j][0])));
                        }
                        if (corr >= 0) {
                              matrixR.m_afData[i][j] = varianceParameter * corr;
                              matrixR.m_afData[j][i] = varianceParameter * corr;
                        }
       break;   

					default:
						throw new UnsupportedOperationException("Matrix.ConstructRMatrix() : This type of correlation structure is not supported in this function");
					}
				}
			}
			return matrixR;
		}
	}
	
	/**
	 * Construct a within-subject correlation matrix using a rho parameter, a gamma parameter, a residual parameter and a column vector of coordinates.
	 * @param coordinates a column vector of coordinates from which the distances are calculated
	 * @param rho the rho parameter
	 * @param gamma the gamma parameter
	 * @param residual the residual parameter
	 * @param type the type of correlation
	 * @return the resulting matrix
	 */
	public static Matrix constructRMatrix(Matrix coordinates, double rho, double gamma, double residual, TypeMatrixR type) {
		if (!coordinates.isColumnVector()) {
			throw new UnsupportedOperationException("Matrix.constructRMatrix() : The coordinates matrix is not a column vector");
		} else {
			int nrow = coordinates.m_iRows;
			Matrix matrixR = new Matrix(nrow,nrow);
			switch(type) {
			case ARMA:		
				double corr = 0d;	
				for (int i = 0; i < nrow; i++) {
					for (int j = i+1; j < nrow; j++) {							
						corr =Math.abs(i - j)-1;						
						double powCol = java.lang.Math.pow(rho, corr);
						matrixR.m_afData[i][j] = residual * gamma * powCol;
						matrixR.m_afData[j][i] = matrixR.m_afData[i][j] ;						

					}
					matrixR.m_afData[i][i]  = residual;
				}
				break;
			default:
				throw new UnsupportedOperationException("Matrix.ConstructRMatrix() : This type of correlation structure is not supported in this function");
			}
			return matrixR;
		}
	}
	
    /**
     * This method generates a random vector
     * @param nrow the number of elements to be generated
     * @param type the distribution type (a Distribution.Type enum variable)
     * @return a Matrix instance
     */
	public static Matrix drawRandomVector(int nrow, Distribution.Type type) {
		return StatisticalUtility.drawRandomVector(nrow, type, StatisticalUtility.getRandom());
	}

	
	/**
	 * This method returns a Random generator.
	 * @return a Random instance
	 */
	public static Random getRandom() {
		if (random == null) {
			random = new Random();
		}
		return random;
	}
	
    /**
     * This method generates a random vector
     * @param nrow the number of elements to be generated
     * @param type the distribution type (a Distribution enum variable)
     * @param random a Random instance
     * @return a Matrix instance
     */
	public static Matrix drawRandomVector(int nrow, Distribution.Type type, Random random) {
		try {
			boolean valid = true;
			Matrix matrix = new Matrix(nrow,1);
			for (int i=0; i<nrow; i++) {
				double number = 0.0f;
				switch (type) {
				case GAUSSIAN:		// Gaussian random number ~ N(0,1)
					number = random.nextGaussian();
					break;
				case UNIFORM:		// Uniform random number [0,1]
					number = random.nextDouble();
					break;
				default:
					i = nrow;
					System.out.println("Matrix.RandomVector() : The specified distribution is not supported in the function");
					valid = false;
					break;
				}
				if (valid) matrix.m_afData[i][0] = number;
			}
			if (valid) return matrix;
			else return null;
		} catch (Exception e) {
			System.out.println("Matrix.RandomVector() : Error while computing the random vector");
			return null;		
		}
	}

	/**
	 * This method performs a special addition in which only the elements different from 0 and 1 
	 * are involved. NOTE: this method is used with SAS output. 
	 * @param originalMatrix the matrix of parameters
	 * @param matrixToAdd the matrix of parameter deviates
	 * @return the new parameters in a new Matrix instance
	 */
	public static Matrix performSpecialAdd(Matrix originalMatrix, Matrix matrixToAdd) {
		Matrix oMat = originalMatrix.getDeepClone();
		List<Integer> oVector = new ArrayList<Integer>();
		oVector.clear();
		
		for (int i = 0; i < originalMatrix.m_iRows; i++) {
			if (oMat.m_afData[i][0] != 0.0 && oMat.m_afData[i][0] != 1.0) { 
				oVector.add(i);
			}
		}

		if (oVector.size() != matrixToAdd.m_iRows) {
			throw new InvalidParameterException("The number of rows do not match!");
		} else {
			for (int j = 0; j < oVector.size(); j++) {
				oMat.m_afData[oVector.get(j)][0] += matrixToAdd.m_afData[j][0];
			}
		}
		return oMat;
	}
	
	/**
	 * This method combines two row vectors of dummy variables. Useful for regressions.
	 * @param mat1 the first row vector
	 * @param mat2 the second row vector
	 * @return the resulting matrix
	 */
	public static Matrix combineMatrices(Matrix mat1, Matrix mat2) {
		if (mat1.m_iRows == mat2.m_iRows) {
			int nbCols = mat1.m_iCols * mat2.m_iCols;
			Matrix oMat = new Matrix(mat1.m_iRows, nbCols);
			for (int i = 0; i < mat1.m_iRows; i++) {
				for (int j = 0; j < mat1.m_iCols; j++) {
					for (int j_prime = 0; j_prime < mat2.m_iCols; j_prime++) {
						oMat.m_afData[i][j*mat2.m_iCols+j_prime] = mat1.m_afData[i][j] * mat2.m_afData[i][j_prime];
					}
				}
			}
			return oMat;
		} else {
			throw new UnsupportedOperationException("The two matrices do not have the same number of rows!");
		}
//		if (mat1.isRowVector() && mat2.isRowVector()) {
//			int nbCols = mat1.m_iCols * mat2.m_iCols;
//			Matrix oMat = new Matrix(1,nbCols);
//			for (int i = 0; i < mat1.m_iCols; i++)
//				for (int j = 0; j < mat2.m_iCols; j++)
//					oMat.m_afData[0][i*mat2.m_iCols+j] = mat1.m_afData[0][i] * mat2.m_afData[0][j];
//			return oMat;
//		} else {
//			throw new UnsupportedOperationException("One of the two matrices is not a row vector!");
//		}
	}

}