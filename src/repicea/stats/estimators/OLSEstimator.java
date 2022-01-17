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
package repicea.stats.estimators;

import java.util.logging.Level;

import repicea.math.Matrix;
import repicea.stats.data.StatisticalDataStructure;
import repicea.stats.estimates.Estimate;
import repicea.stats.estimates.GaussianEstimate;
import repicea.stats.estimates.VarianceEstimate;
import repicea.stats.model.AbstractStatisticalModel;
import repicea.stats.model.StatisticalModel;
import repicea.stats.model.lm.LinearModel;
import repicea.util.REpiceaLogManager;

/**
 * The OLSOptimizer implements the Ordinary Least Squares estimator.
 * @author Mathieu Fortin - November 2012
 */
public class OLSEstimator implements Estimator {

	private VarianceEstimate residualVariance;
	private boolean hasConverged;
	private Estimate<?> betaVector;
	
	@Override
	public boolean doEstimation(StatisticalModel<? extends StatisticalDataStructure> model)	throws EstimatorException {
		if (!(model instanceof LinearModel)) {
			throw new EstimatorException("The OLS optimizer is designed to work with instances of LinearModel only!");
		}
		Matrix matrixX = model.getDataStructure().getMatrixX();
		Matrix matrixY = model.getDataStructure().getVectorY();
		Matrix matrixXT = matrixX.transpose();
		betaVector = new GaussianEstimate();
		Matrix inverseProduct = matrixXT.multiply(matrixX).getInverseMatrix();
		((GaussianEstimate) betaVector).setMean(inverseProduct.multiply(matrixX.transpose()).multiply(matrixY));
		model.setParameters(betaVector.getMean());
		Matrix residual = model.getResiduals();
		int degreesOfFreedom = model.getDataStructure().getNumberOfObservations() - betaVector.getMean().m_iRows;
		double sse = residual.transpose().multiply(residual).getValueAt(0, 0);
		REpiceaLogManager.logMessage(AbstractStatisticalModel.LoggerName,
				Level.FINE,
				null, 
				"Sum of squared errors = " + sse);
		double resVar = sse / degreesOfFreedom;
		residualVariance = new VarianceEstimate(degreesOfFreedom, resVar);
		((GaussianEstimate) betaVector).setVariance(inverseProduct.scalarMultiply(resVar));
		return true;
	}

	/**
	 * This method returns the residual variance of the OLS algorithm.
	 * @return a VarianceEstimate instance
	 */
	public VarianceEstimate getResidualVariance() {
		return residualVariance;
	}

	@Override
	public boolean isConvergenceAchieved() {return hasConverged;}

	@Override
	public Estimate<?> getParameterEstimates() {return betaVector;}


}
