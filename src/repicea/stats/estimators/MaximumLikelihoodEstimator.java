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

import java.util.ArrayList;
import java.util.List;

import repicea.math.Matrix;
import repicea.math.optimizer.AbstractOptimizer.OptimizationException;
import repicea.stats.data.StatisticalDataStructure;
import repicea.stats.estimates.Estimate;
import repicea.stats.estimates.GaussianEstimate;
import repicea.stats.model.CompositeLogLikelihood;
import repicea.stats.model.StatisticalModel;

/**
 * This interface specifies if the optimizer is based on a maximum likelihood theory.
 * @author Mathieu Fortin - August 2011
 */
public class MaximumLikelihoodEstimator implements Estimator {
	
	protected boolean hasConverged;
	protected GaussianEstimate parameterEstimate;
	
	protected final repicea.math.optimizer.NewtonRaphsonOptimizer nro;

	public MaximumLikelihoodEstimator() {
		nro = new repicea.math.optimizer.NewtonRaphsonOptimizer();
	}
	
	
	@Override
	public boolean doEstimation(StatisticalModel<? extends StatisticalDataStructure> model) throws EstimatorException {
		nro.setConvergenceCriterion(model.getConvergenceCriterion());
		
		CompositeLogLikelihood llk = model.getCompleteLogLikelihood();
		List<Integer> indices = new ArrayList<Integer>();
		for (int i = 0; i < model.getParameters().m_iRows; i++) {
			indices.add(i);
		}
		try {
			nro.optimize(llk, indices);
		} catch (OptimizationException e) {
			System.out.println("Newton-Raphson optimisation failed.");
			parameterEstimate = null;
			return false;
		}
		if (nro.isConvergenceAchieved()) {
			Matrix varCov = nro.getHessianAtMaximum().getInverseMatrix().scalarMultiply(-1d);
			parameterEstimate = new GaussianEstimate(nro.getParametersAtMaximum(), varCov);
			return true;
		} else {
			parameterEstimate = null;
			return false;
		}
		
	}

	/**
	 * Returns the maximum log likelihood value after convergence.
	 * @return a double
	 */
	public double getMaximumLogLikelihood() {return nro.getOptimalValue();}


	@Override
	public boolean isConvergenceAchieved() {return nro.isConvergenceAchieved();}

	@Override
	public Estimate<?> getParameterEstimates() {
		return parameterEstimate;
	}
	
	@Override
	public String toString() {return "Maximum likelihood estimator";}


	@Override
	public void setVerboseEnabled(boolean bool) {
		nro.setVerboseEnabled(bool);
	}
	
	

}
