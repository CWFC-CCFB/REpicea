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
import repicea.stats.data.DataSet;
import repicea.stats.data.GenericHierarchicalStatisticalDataStructure;
import repicea.stats.data.HierarchicalStatisticalDataStructure;
import repicea.stats.data.StatisticalDataException;
import repicea.stats.estimators.Estimator;
import repicea.stats.estimators.MaximumLikelihoodEstimator;
import repicea.stats.model.AbstractStatisticalModel;
import repicea.stats.model.CompositeLogLikelihood;
import repicea.stats.model.IndividualLogLikelihood;
import repicea.stats.model.glm.LinkFunction.Type;

/**
 * This class implements generalized linear models. 
 * @author Mathieu Fortin - August 2011
 */
public class GeneralizedLinearModel extends AbstractStatisticalModel<HierarchicalStatisticalDataStructure> {
	
	protected LinkFunction.Type linkFunctionType;
	protected Matrix matrixX;		// reference
	protected Matrix y;				// reference

	protected IndividualLogLikelihood individualLLK;
	

	/**
	 * General constructor
	 * @param dataSet the fitting data
	 * @param linkFunctionType the type of ling function (Logit, CLogLog, ...)
	 * @param modelDefinition a String that defines the dependent variable and the effects of the model
	 * @param startingBeta the starting values of the parameters
	 */
	public GeneralizedLinearModel(DataSet dataSet, Type linkFunctionType, String modelDefinition, Matrix startingBeta) {
		super(dataSet);

		// then define the model effects and retrieve matrix X and vector y
		try {
			setModelDefinition(modelDefinition);
		} catch (StatisticalDataException e) {
			System.out.println("Unable to define this model : " + modelDefinition);
			e.printStackTrace();
		}
		initializeLinkFunction(linkFunctionType);
		this.linkFunctionType = linkFunctionType;
		// first initialize the model structure
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
	 * Constructor for derived class.
	 */
	protected GeneralizedLinearModel(GeneralizedLinearModel glm) {
		this(glm.getDataStructure().getDataSet(), glm.getLinkFunctionType(), glm.getModelDefinition());
		this.matrixX = glm.matrixX;
		this.y = glm.y;

		this.individualLLK = glm.individualLLK;
	}

	/**
	 * This method returns the type of the link function.
	 * @return a LinkFunction.Type enum variable
	 */
	public LinkFunction.Type getLinkFunctionType() {
		return linkFunctionType;
	}
	

	@Override
	protected void setModelDefinition(String modelDefinition) throws StatisticalDataException {
		super.setModelDefinition(modelDefinition);
		matrixX = getDataStructure().getMatrixX();
		y = getDataStructure().getVectorY();
	}
	
	protected void initializeLinkFunction(Type linkFunctionType) {
//		LinearStatisticalExpression le = new LinearStatisticalExpression();
		LinkFunction lf = new LinkFunction(linkFunctionType);
//		lf.setParameterValue(LFParameter.Eta, le);
		
		individualLLK = new IndividualLogLikelihood(new LikelihoodGLM(lf));
		setCompleteLLK();
	}

	
	@Override
	public void setParameters(Matrix beta) {
		if (beta == null) {
			individualLLK.setBeta(new Matrix(matrixX.m_iCols, 1));		// default starting parameters at 0
		} else {
			individualLLK.setBeta(beta);
		}
	}
		
	@Override
	public Matrix getParameters() {return individualLLK.getBeta();}

	@Override
	protected Estimator instantiateDefaultEstimator() {
		return new MaximumLikelihoodEstimator();
	}

	@Override
	public String toString() {
		return "Generalized linear model";
	}

	
	public Matrix getLinearPredictions() {
		if (getEstimator().isConvergenceAchieved()) {
			int nbObs = getDataStructure().getNumberOfObservations();
			Matrix pred = new Matrix(getDataStructure().getNumberOfObservations(),2);
			Matrix beta = getEstimator().getParameterEstimates().getMean().getSubMatrix(0, matrixX.m_iCols - 1, 0, 0);
			Matrix linearPred = matrixX.multiply(beta);
			Matrix omega = getEstimator().getParameterEstimates().getVariance().getSubMatrix(0, matrixX.m_iCols - 1, 0, matrixX.m_iCols - 1);
			for (int i = 0; i < nbObs; i++) {
				pred.setValueAt(i, 0, linearPred.getValueAt(i, 0));
				Matrix x_i = matrixX.getSubMatrix(i, i, 0, matrixX.m_iCols - 1);
				pred.setValueAt(i, 1, x_i.multiply(omega).multiply(x_i.transpose()).getValueAt(0, 0));
			}
			return pred;
		} else {
			return null;
		}
	}

	
	@Override
	public Matrix getPredicted() {
		if (getEstimator().isConvergenceAchieved()) {
			setParameters(getEstimator().getParameterEstimates().getMean());
//			getCompleteLogLikelihood().setBeta(getEstimator().getParameterEstimates().getMean());
			return getCompleteLogLikelihood().getPredictions();
		} else {
			return null;
		}
	}

	@Override
	public Matrix getResiduals() {
		if (getEstimator().isConvergenceAchieved()) {
			return y.subtract(getPredicted());
		} else {
			return null;
		}
	}

	@Override
	protected HierarchicalStatisticalDataStructure getDataStructureFromDataSet(DataSet dataSet) {
		return new GenericHierarchicalStatisticalDataStructure(dataSet);
	}

	@Override
	protected void setCompleteLLK() {completeLLK = new CompositeLogLikelihood(individualLLK, matrixX, y);}

}
