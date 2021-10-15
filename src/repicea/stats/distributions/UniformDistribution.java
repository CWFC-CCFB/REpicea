/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2019 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.distributions;

import java.security.InvalidParameterException;

import repicea.math.Matrix;
import repicea.stats.Distribution;
import repicea.stats.StatisticalUtility;

@SuppressWarnings("serial")
public class UniformDistribution implements Distribution, BoundedDistribution {

	private final BasicBound upperBound;
	private final BasicBound lowerBound;

	/**
	 * Constructor. </br>
	 * </br>
	 * The two arguments must column vectors of the same size.
	 * @param lowerBound a Matrix instance that represents a column vector.
	 * @param upperBound a Matrix instance that represents a column vector.
	 */
	public UniformDistribution(Matrix lowerBound, Matrix upperBound) {
		if (lowerBound == null || upperBound == null) {
			throw new InvalidParameterException("The lowerBound and upperBound matrices must be column vectors!");
		}
		if (!lowerBound.isColumnVector() || !upperBound.isColumnVector()) {
			throw new InvalidParameterException("The lowerBound and upperBound matrices must be column vectors!");
		} 
		if (lowerBound.m_iCols != upperBound.m_iCols) {
			throw new InvalidParameterException("The lowerBound and upperBound matrices must be column vectors of the same dimension!");
		}
		if (upperBound.subtract(lowerBound).anyElementSmallerOrEqualTo(0d)) {	// means the upper bound is equal or smaller than the lower bound
			throw new InvalidParameterException("The elements in the upperBound matrix must be larger than their corresponding element in the lowerBound matrix!");
		}
		this.upperBound = new BasicBound(true);
		this.upperBound.setBoundValue(upperBound);
		this.lowerBound = new BasicBound(false);
		this.lowerBound.setBoundValue(lowerBound);
	}
	
	
	
	@Override
	public Matrix getMean() {
		return upperBound.getBoundValue().add(lowerBound.getBoundValue()).scalarMultiply(.5);
	}

	@Override
	public Matrix getVariance() {
		Matrix diagonalDifference = upperBound.getBoundValue().subtract(lowerBound.getBoundValue()).matrixDiagonal();
		return diagonalDifference.multiply(diagonalDifference).scalarMultiply(1d/12);
	}

	@Override
	public boolean isParametric() {return true;}

	@Override
	public boolean isMultivariate() {return getMean().m_iRows > 1;}

	@Override
	public Type getType() {return Type.UNIFORM;}

	@Override
	public Matrix getRandomRealization() {
		Matrix diagonalDifference = upperBound.getBoundValue().subtract(lowerBound.getBoundValue()).matrixDiagonal();
		Matrix deviates = StatisticalUtility.drawRandomVector(getMean().m_iRows, Type.UNIFORM);
		return lowerBound.getBoundValue().add(diagonalDifference.multiply(deviates));
	}

	/**
	 * Return the probability density of a column vector.
	 * @param values a Matrix instance that stands for a column vector.
	 * @return the composite probability density
	 */
	public double getProbabilityDensity(Matrix values) {
		if (values == null) {
			throw new InvalidParameterException("The values parameter must be non null!");
		}
		if (values.m_iCols != upperBound.getBoundValue().m_iCols) {
			throw new InvalidParameterException("The number of elements in the values parameter is inconsistent with the bounds!");
		}
		Matrix densities = upperBound.getBoundValue().subtract(lowerBound.getBoundValue()).elementWisePower(-1d);
		double product = 1d;
		for (int i = 0; i < densities.m_iCols; i++) {
			product *= densities.getValueAt(i, 0);
		}
		return product;
	}



	@Override
	public void setLowerBoundValue(Matrix lowerBoundValue) {
		lowerBound.setBoundValue(lowerBoundValue);
	}

	@Override
	public void setUpperBoundValue(Matrix upperBoundValue) {
		upperBound.setBoundValue(upperBoundValue);
	}
	
	public BasicBound getLowerBound() {return lowerBound;}
	public BasicBound getUpperBound() {return upperBound;}

}
