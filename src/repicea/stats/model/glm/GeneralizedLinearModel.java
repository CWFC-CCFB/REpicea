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
package repicea.stats.model.glm;


import repicea.math.Matrix;
import repicea.stats.LinearStatisticalExpression;
import repicea.stats.data.DataSet;
import repicea.stats.data.GenericHierarchicalStatisticalDataStructure;
import repicea.stats.data.HierarchicalStatisticalDataStructure;
import repicea.stats.data.StatisticalDataException;
import repicea.stats.estimators.Estimator;
import repicea.stats.estimators.MaximumLikelihoodEstimator;
import repicea.stats.model.AbstractStatisticalModel;
import repicea.stats.model.CompositeLogLikelihood;
import repicea.stats.model.LogLikelihood;
import repicea.stats.model.glm.LinkFunction.LFParameter;
import repicea.stats.model.glm.LinkFunction.Type;

/**
 * This class implements generalized linear models. 
 * @author Mathieu Fortin - August 2011
 */
public class GeneralizedLinearModel extends AbstractStatisticalModel<HierarchicalStatisticalDataStructure> {

	
//	/**
//	 * This class defines the log likelihood function of the model.
//	 * @author Mathieu Fortin - August 2011
//	 */
//	@SuppressWarnings("serial")
//	protected static class OverallLogLikelihood extends LogLikelihood {
//
//		private GeneralizedLinearModel model;
//		
//		protected OverallLogLikelihood(GeneralizedLinearModel generalizedLinearModel) {
//			super(generalizedLinearModel.individualLLK);
//			this.model = generalizedLinearModel;
//			init();
//		}
//
//		protected void init() {
//			
//		}
//		
//		@Override
//		public Matrix getGradient() {
//			int numberParameters = getModel().le.getNumberOfParameters();
//			Matrix gradient = new Matrix(numberParameters, 1);
//			for (int i = 0; i < getModel().matrixX.m_iRows; i++) {
//				getModel().setObservationInLogLikelihoodFunction(i);
//				Matrix tmp = getModel().individualLLK.getGradient();
//				MatrixUtility.add(gradient, tmp);		// report the result in Matrix gradient
//			}
//			return gradient;
//		}
//
//		@Override
//		public Matrix getHessian() {
//			int numberParameters = getModel().le.getNumberOfParameters();
//			Matrix hessian = new Matrix(numberParameters, numberParameters);
//			for (int i = 0; i < getModel().matrixX.m_iRows; i++) {
//				getModel().setObservationInLogLikelihoodFunction(i);
//				MatrixUtility.add(hessian, getModel().individualLLK.getHessian());		// report the result in Matrix hessian
//			}
//			return hessian;
//		}
//
//		@Override
//		public Double getValue() {
//			double loglikelihood = 0;
//			for (int i = 0; i < getModel().matrixX.m_iRows; i++) {
//				getModel().setObservationInLogLikelihoodFunction(i);
//				loglikelihood += getModel().individualLLK.getValue();
//			}
//			return loglikelihood;
//		}
//		
//		protected GeneralizedLinearModel getModel() {return model;}
//		
//	}
	
	protected LinkFunction.Type linkFunctionType;
	protected Matrix matrixX;		// reference
	protected Matrix y;				// reference

	protected LogLikelihood individualLLK;
	protected LinearStatisticalExpression le;
	

	/**
	 * General constructor
	 * @param dataSet the fitting data
	 * @param linkFunctionType the type of ling function (Logit, CLogLog, ...)
	 * @param modelDefinition a String that defines the dependent variable and the effects of the model
	 * @param startingBeta the starting values of the parameters
	 */
	public GeneralizedLinearModel(DataSet dataSet, Type linkFunctionType, String modelDefinition, Matrix startingBeta) {
		super(dataSet);
		this.linkFunctionType = linkFunctionType;
		// first initialize the model structure
		initializeLinkFunction(linkFunctionType);

		// then define the model effects and retrieve matrix X and vector y
		try {
			setModelDefinition(modelDefinition);
		} catch (StatisticalDataException e) {
			System.out.println("Unable to define this model : " + modelDefinition);
			e.printStackTrace();
		}
		setParameters(startingBeta);
	}

	/**
	 * Constructor using a vector of 0s as starting values for the parameters
	 * @param dataSet the fitting data
	 * @param linkFunctionType the type of ling function (Logit, CLogLog, ...)
	 * @param modelDefinition a String that defines the dependent variable and the effects of the model
	 */
	public GeneralizedLinearModel(DataSet dataSet, Type linkFunctionType, String modelDefinition) {
		this(dataSet, linkFunctionType, modelDefinition, null);
	}

	/**
	 * This method returns the type of the link function.
	 * @return a LinkFunction.Type enum variable
	 */
	public LinkFunction.Type getLinkFunctionType() {
		return linkFunctionType;
	}
	
	/**
	 * Constructor for derived class.
	 */
	protected GeneralizedLinearModel(GeneralizedLinearModel glm) {
		this(glm.getDataStructure().getDataSet(), glm.getLinkFunctionType(), glm.getModelDefinition());
		this.matrixX = glm.matrixX;
		this.y = glm.y;

		this.individualLLK = glm.individualLLK;
		this.le = glm.le;
		
		setOverallLLK();
	}
	
	
	
	
	
	
//	/**
//	 * This method sets the variable and parameter in the individual log-likelihood function to the observation i.
//	 * @param i the index of the observation
//	 */
//	protected void setObservationInLogLikelihoodFunction(int i) {
//		le.setX(matrixX.getSubMatrix(i, i, 0, matrixX.m_iCols - 1));
//		((LikelihoodGLM) individualLLK.getOriginalFunction()).setObservedValue(y.m_afData[i][0] == 1d);
//	}

	@Override
	protected void setModelDefinition(String modelDefinition) throws StatisticalDataException {
		super.setModelDefinition(modelDefinition);
		matrixX = getDataStructure().getMatrixX();
		y = getDataStructure().getVectorY();
	}
	
	protected void initializeLinkFunction(Type linkFunctionType) {
		le = new LinearStatisticalExpression();
		LinkFunction lf = new LinkFunction(linkFunctionType);
		lf.setParameterValue(LFParameter.Eta, le);
		
		individualLLK = new LogLikelihood(new LikelihoodGLM(lf));
//		individualLLK.setParameterValue(LLKGLMParameter.Predicted, lf);
		setOverallLLK();
	}

	
	@Override
	public void setParameters(Matrix beta) {
		if (beta == null) {
			le.setBeta(new Matrix(matrixX.m_iCols, 1));		// default starting parameters at 0
		} else {
			le.setBeta(beta);
		}
	}
		
	@Override
	public Matrix getParameters() {return le.getBeta();}

	@Override
	protected Estimator instantiateDefaultOptimizer() {
		return new MaximumLikelihoodEstimator();
	}

	@Override
	protected void setOverallLLK() {overallLLK = new CompositeLogLikelihood(this);}

	@Override
	public String toString() {
		return "Generalized linear model";
	}

	
	public Matrix getLinearPredictions() {
		if (getOptimizer().isConvergenceAchieved()) {
			int nbObs = getDataStructure().getNumberOfObservations();
			Matrix pred = new Matrix(getDataStructure().getNumberOfObservations(),2);
			Matrix beta = getOptimizer().getParameterEstimates().getMean().getSubMatrix(0, matrixX.m_iCols - 1, 0, 0);
			Matrix linearPred = matrixX.multiply(beta);
			Matrix omega = getOptimizer().getParameterEstimates().getVariance().getSubMatrix(0, matrixX.m_iCols - 1, 0, matrixX.m_iCols - 1);
			for (int i = 0; i < nbObs; i++) {
				pred.m_afData[i][0] = linearPred.m_afData[i][0];
				Matrix x_i = matrixX.getSubMatrix(i, i, 0, matrixX.m_iCols - 1);
				pred.m_afData[i][1] = x_i.multiply(omega).multiply(x_i.transpose()).m_afData[0][0];
			}
			return pred;
		} else {
			return null;
		}
	}

	
	@Override
	public Matrix getPredicted() {
		if (getOptimizer().isConvergenceAchieved()) {
			int nbObs = getDataStructure().getNumberOfObservations();
			Matrix pred = new Matrix(getDataStructure().getNumberOfObservations(),1);
			setParameters(getOptimizer().getParameterEstimates().getMean());
			for (int i = 0; i < nbObs; i++) {
				setObservationInLogLikelihoodFunction(i);
				pred.m_afData[i][0] = ((LikelihoodGLM) individualLLK.getOriginalFunction()).getPrediction();
			}
			return pred;
		} else {
			return null;
		}
	}

	@Override
	public Matrix getResiduals() {
		if (getOptimizer().isConvergenceAchieved()) {
			return y.subtract(getPredicted());
		} else {
			return null;
		}
	}

	@Override
	protected HierarchicalStatisticalDataStructure getDataStructureFromDataSet(DataSet dataSet) {
		return new GenericHierarchicalStatisticalDataStructure(dataSet);
	}

}
