/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2021 Mathieu Fortin for Rouge Epicea.
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
package repicea.simulation.metamodel;

import java.util.ArrayList;

import repicea.math.Matrix;
import repicea.stats.data.StatisticalDataException;
import repicea.stats.distributions.GaussianDistribution;
import repicea.stats.distributions.UniformDistribution;

/**
 * An implementation of the derivative form of the Chapman-Richards model including random effects.
 * @author Mathieu Fortin - October 2021
 */
class ChapmanRichardsDerivativeModelWithRandomEffectImplementation extends AbstractMixedModelImplementation {
	
	ChapmanRichardsDerivativeModelWithRandomEffectImplementation(String outputType, MetaModel model) throws StatisticalDataException {
		super(outputType, model);
	}

	
	@Override
	GaussianDistribution getStartingParmEst(double coefVar) {
		Matrix parmEst = new Matrix(5,1);
		parmEst.setValueAt(0, 0, 1000d);
		parmEst.setValueAt(1, 0, 0.02);
		parmEst.setValueAt(2, 0, 2d);
		parmEst.setValueAt(3, 0, 1000d);
		parmEst.setValueAt(4, 0, .92);
		
		fixedEffectsParameterIndices = new ArrayList<Integer>();
		fixedEffectsParameterIndices.add(0);
		fixedEffectsParameterIndices.add(1);
		fixedEffectsParameterIndices.add(2);
		
		this.indexRandomEffectVariance = 3;
		this.indexCorrelationParameter = 4;
		
		Matrix varianceDiag = new Matrix(parmEst.m_iRows,1);
		for (int i = 0; i < varianceDiag.m_iRows; i++) {
			varianceDiag.setValueAt(i, 0, Math.pow(parmEst.getValueAt(i, 0) * coefVar, 2d));
		}
		
		GaussianDistribution gd = new GaussianDistribution(parmEst, varianceDiag.matrixDiagonal());
		
		
		Matrix lowerBound = new Matrix(5,1);
		Matrix upperBound = new Matrix(5,1);
		lowerBound.setValueAt(0, 0, 0);
		upperBound.setValueAt(0, 0, 2000);
		
		lowerBound.setValueAt(1, 0, 0.00001);
		upperBound.setValueAt(1, 0, 0.05);
		
		lowerBound.setValueAt(2, 0, 1);
		upperBound.setValueAt(2, 0, 6);

		lowerBound.setValueAt(3, 0, 0);
		upperBound.setValueAt(3, 0, 2000);

		lowerBound.setValueAt(4, 0, .90);
		upperBound.setValueAt(4, 0, .99);

		priors = new UniformDistribution(lowerBound, upperBound);

		return gd;
	}


	@Override
	double getPrediction(double ageYr, double timeSinceBeginning, double r1) {
		double b1 = getParameters().getValueAt(0, 0);
		double b2 = getParameters().getValueAt(1, 0);
		double b3 = getParameters().getValueAt(2, 0);
		double pred = (b1 + r1) * Math.exp(-b2 * ageYr) * Math.pow(1 - Math.exp(-b2 * ageYr), b3);
		return pred;
	}


	@Override
	Matrix getFirstDerivative(double ageYr, double timeSinceBeginning, double r1) {
		double b1 = getParameters().getValueAt(0, 0);
		double b2 = getParameters().getValueAt(1, 0);
		double b3 = getParameters().getValueAt(2, 0);
		
		double exp = Math.exp(-b2 * ageYr);
		double root = 1 - exp;
		
		Matrix derivatives = new Matrix(3,1);
		derivatives.setValueAt(0, 0, exp * Math.pow(root, b3));
		derivatives.setValueAt(1, 0, - ageYr * b1 * exp * Math.pow(root, b3) + 
				b1 * exp * b3 * Math.pow(root, b3 - 1) * exp * ageYr);
		derivatives.setValueAt(2, 0, b1 * exp * Math.pow(root, b3) * Math.log(root));
		return derivatives;
	}

	
}
