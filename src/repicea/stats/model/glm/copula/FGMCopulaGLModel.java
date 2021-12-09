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
package repicea.stats.model.glm.copula;

import java.security.InvalidParameterException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;

import repicea.math.Matrix;
import repicea.stats.data.DataSet;
import repicea.stats.data.GenericHierarchicalSpatialDataStructure;
import repicea.stats.data.HierarchicalStatisticalDataStructure;
import repicea.stats.data.StatisticalDataException;
import repicea.stats.estimates.CorrelationEstimate;
import repicea.stats.model.glm.GeneralizedLinearModel;

/**
 * The FGMCopulaGLModel class implements a Farlie-Gumbel-Morgenstein (FGM) copula for binary outcomes. The method is described in 
 * Bhat and Sener 2009.
 * @see <a href=http://www.springerlink.com/content/t33474476r231w60/> 
 * Bhat, C.R., and Sener, I. N. 2009. A copula-based closed form binary logit choice model for accomodating spatial correlation across observational units. Journal of Geographical Systems 11(3): 243-272.
 * </a>
 * @author Mathieu Fortin - June 2011
 */
public class FGMCopulaGLModel extends GeneralizedLinearModel {
	
	@SuppressWarnings("rawtypes")
	protected static class LikelihoodValue implements Comparable {

		private double llk;
		private Matrix beta;
		
		protected LikelihoodValue(Matrix beta, double llk) {
			this.beta = beta.getDeepClone();
			this.llk = llk;
		}
		
		@Override
		public int compareTo(Object arg0) {
			double reference = ((LikelihoodValue) arg0).llk;
			if (this.llk < reference) {
				return 1;
			} else if (this.llk == reference) {
				return 0;
			} else {
				return -1;
			}
		}
		
		protected Matrix getParameters() {return beta;}
	}
	
	
	
	
	private CopulaExpression copula;
	
	/**
	 * Constructor for this class
	 * @param glm a GeneralizedLinearModel instance
	 * @param copula a CopulaExpression instance
	 * @throws StatisticalDataException if the hierarchical level specification in the copula is not found in the data set
	 */
	public FGMCopulaGLModel(GeneralizedLinearModel glm, CopulaExpression copula) throws StatisticalDataException {
		super(glm);
		if (!glm.getEstimator().isConvergenceAchieved()) {
			glm.doEstimation();
		}
		glm.setParameters(glm.getEstimator().getParameterEstimates().getMean());
		this.copula = copula;
		this.copula.initialize(this, getDataStructure());
		setCompleteLLK();
	}
	
	@Override
	public Matrix getParameters() {return individualLLK.getBeta().matrixStack(copula.getBeta(), true);}
	

	@Override
	public void setParameters(Matrix beta) {
		if (beta == null) {
			individualLLK.setBeta(new Matrix(matrixX.m_iCols, 1));		// default starting parameters at 0
		} else {
			individualLLK.setBeta(beta.getSubMatrix(0, beta.m_iRows - copula.getNumberOfParameters() - 1, 0, 0));
			copula.setBeta(beta.getSubMatrix(beta.m_iRows - copula.getNumberOfParameters(), beta.m_iRows - 1, 0, 0));
		}
	}

	/**
	 * This method scans the log likelihood function within a range of values for a particular parameter.
	 * @param parameterName the index of the parameter
	 * @param start the starting value
	 * @param end the ending value
	 * @param step the step between these two values.
	 */
	@SuppressWarnings("unchecked")
	public void gridSearch(int parameterName, double start, double end, double step) {
		if (start >= end) {
			throw new InvalidParameterException("The argument end must be greater than the start argument!");
		}
		System.out.println("Initializing grid search...");
		ArrayList<LikelihoodValue> likelihoodValues = new ArrayList<LikelihoodValue>();
		Matrix originalParameters = getParameters();
		double llk;
		for (double value = start; value < end + step; value+=step) {
			Matrix beta = originalParameters.getDeepClone();
			beta.setValueAt(parameterName, 0, value);
			setParameters(beta);
			((FGMCompositeLogLikelihood) getCompleteLogLikelihood()).reset();
			llk = getCompleteLogLikelihood().getValue();
			likelihoodValues.add(new LikelihoodValue(beta, llk));
			System.out.println("Parameter value : " + value + "; Log-likelihood : " + llk);
		}
		
		Collections.sort(likelihoodValues);
		LikelihoodValue lk;
		Matrix bestFittingParameters = null;
		for (int i = 0; i < likelihoodValues.size(); i++) {
			lk = likelihoodValues.get(i);
			if (!Double.isNaN(lk.llk)) {
				bestFittingParameters = lk.getParameters();
				break;
			}
		}
		if (bestFittingParameters == null) {
			throw new InvalidParameterException("All the likelihoods of the grid are NaN!");
		} else {
			setParameters(bestFittingParameters);
		}
	}

	@Override
	protected void setCompleteLLK() {completeLLK = new FGMCompositeLogLikelihood(individualLLK,	matrixX, y,	getDataStructure(),	copula);}
	
	@Override
	public String toString() {
		return "Generalized linear model based on FGM copula";
	}

	@Override
	public HierarchicalStatisticalDataStructure getDataStructure() {
		return (HierarchicalStatisticalDataStructure) super.getDataStructure();
	}
	
	@Override
	public void getSummary() {
		super.getSummary();
		if (getEstimator().isConvergenceAchieved()) {
			NumberFormat formatter = NumberFormat.getInstance();
			formatter.setMaximumFractionDigits(4);
			SpearmanCorrelationCoefficient scc = new SpearmanCorrelationCoefficient();
			CorrelationEstimate[] spearmanCorr = scc.getSpearmanCorrelationCoefficient(this);
			System.out.println("Spearman's correlation coefficients");
			for (int i = 0; i < spearmanCorr.length; i++) {
				CorrelationEstimate estimate = spearmanCorr[i];
				System.out.println(i + " : " + "Estimate : " + formatter.format(estimate.getMean().getValueAt(0, 0)) +
						"; n : " + estimate.getSampleSize() + 
						"; t-value : " + formatter.format(estimate.getStudentT()));
			}
		}
	}
	
	protected CopulaExpression getCopula() {return copula;}
	
	
	@Override
	protected HierarchicalStatisticalDataStructure getDataStructureFromDataSet(	DataSet dataSet) {
		return new GenericHierarchicalSpatialDataStructure(dataSet);
	}

	
	
}
