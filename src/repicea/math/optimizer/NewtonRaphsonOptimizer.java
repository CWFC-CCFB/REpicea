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
package repicea.math.optimizer;

import java.util.List;

import repicea.math.AbstractMathematicalFunction;
import repicea.math.Matrix;

/**
 * The NewtonRaphsonOptimizer class implements the Optimizer interface. It optimizes a log-likelihood function using the
 * Newton-Raphson algorithm. The vector of parameter changes is estimated as - inv(d2LLK/d2Beta) * dLLK/dBeta.
 * @author Mathieu Fortin - June 2011
 */
public class NewtonRaphsonOptimizer extends AbstractOptimizer {

	public final static String InnerIterationStarted = "InnerIterationStarted";
//	public final static String OuterIterationStarted = "OuterIterationStarted";
	
	protected int maxNumberOfIterations = 20;
	protected double gradientCriterion = 1E-3;
	private int iterationID;
	
	public NewtonRaphsonOptimizer() {}

	/**
	 * This method optimize the log-likelihood function using the Newton-Raphson optimisation step.
	 * @param model the model that provides the log-likelihood function
	 * @param originalBeta the vector that contains the parameters of the previous outer optimisation
	 * @param optimisationStep the optimisation step from the Newton-Raphson algorithm
	 * @param previousLogLikelihood the value of the log-likelihood function in the last outer optimisation
	 * @throws OptimisationException if the inner optimisation is not successful
	 */
	protected double runInnerOptimisation(AbstractMathematicalFunction<Integer, Double, Integer, Double> function, 
			List<Integer> indicesOfParametersToOptimize,
			Matrix originalBeta, 
			Matrix optimisationStep,
			double previousLogLikelihood,
			AbstractOptimizer.LineSearchMethod lineSearchMethod) throws OptimizationException {
		
		int numberSubIter = 0;
		int maxNumberOfSubiterations = (Integer) lineSearchMethod.getParameter();
		
		double scalingFactor = 1d;
		double currentLlkValue;
		
		do {
			fireOptimizerEvent(NewtonRaphsonOptimizer.InnerIterationStarted);
			Matrix newBeta = originalBeta.add(optimisationStep.scalarMultiply(scalingFactor - numberSubIter * .1));
			setParameters(function, indicesOfParametersToOptimize, newBeta);
			currentLlkValue = function.getValue();
			if (isVerboseEnabled()) {
				System.out.println("    Subiteration : " + numberSubIter + "; Log-likelihood : " + currentLlkValue);
			}
			numberSubIter++;
		} while ((Double.isNaN(currentLlkValue) || currentLlkValue < previousLogLikelihood) && numberSubIter < maxNumberOfSubiterations); // loop if the number of iterations is not over the maximum number and either the likelihood is still higher or non defined
		
		if (Double.isNaN(currentLlkValue) ||  currentLlkValue < previousLogLikelihood) {
			throw new OptimizationException("Failed to improve the log-likelihood function !");
		} else {
			return currentLlkValue;
		}
	}
	
	@Override
	public boolean optimize(AbstractMathematicalFunction<Integer, Double, Integer, Double> function, List<Integer> indicesOfParametersToOptimize) throws OptimizationException {

		if (function instanceof OptimizerListener) {
			addOptimizerListener((OptimizerListener) function);
		}
		
		fireOptimizerEvent(OptimizerListener.optimizationStarted);
		
		double llkValue0 = function.getValue();
		Matrix gradient = function.getGradient();
		Matrix hessian = function.getHessian();
				
		iterationID = 0;
		
		double gconv = calculateConvergence(gradient, hessian, llkValue0);
		
		convergenceAchieved = false;
		if (Math.abs(gconv) < convergenceCriterion) {
			convergenceAchieved = true;
		}
		Matrix currentBeta = null;
		try {
			while (!convergenceAchieved && iterationID <= maxNumberOfIterations) {
				iterationID++;
				Matrix optimisationStep = hessian.getInverseMatrix().multiply(gradient).scalarMultiply(-1d);
				
				Matrix originalBeta = extractParameters(function,indicesOfParametersToOptimize);

				llkValue0 = runInnerOptimisation(function, 
						indicesOfParametersToOptimize, 
						originalBeta, 
						optimisationStep, 
						llkValue0, 
						LineSearchMethod.TEN_EQUAL);		// if it does not throw an Exception, it means the inner optimisation was successful. 
				gradient = function.getGradient();
				hessian = function.getHessian();
				currentBeta = extractParameters(function, indicesOfParametersToOptimize);
				gconv = calculateConvergence(gradient, hessian, llkValue0);
				
				if (gconv < 0) {
					convergenceAchieved = !gradient.getAbsoluteValue().anyElementLargerThan(1E-5);
					System.out.println("GConv is negative!");
					throw new OptimizationException("Convergence could not be achieved after " + ((Integer) iterationID).toString() + " iterations!");
				} else if (gconv < convergenceCriterion) {
					convergenceAchieved = true;
				}

				System.out.println("Iteration : " + iterationID + "; Log-likelihood : " + llkValue0 + "; df : " + gconv + "; parms : " + currentBeta.toString());
			}  

			if (iterationID > maxNumberOfIterations && !convergenceAchieved) {
				throw new OptimizationException("Convergence could not be achieved after " + ((Integer) iterationID).toString() + " iterations!");
			}

			if (gradient.getAbsoluteValue().anyElementLargerThan(gradientCriterion)) {
				System.out.println("A least one element of the gradient vector is larger than 1E-3.");
			}

			optimalValue = llkValue0;

			fireOptimizerEvent(OptimizerListener.optimizationEnded);

			return convergenceAchieved;
		} catch (OptimizationException e) {
			throw e;
		} finally {
			betaVector = currentBeta;
			hessianMatrix = hessian;
			if (function instanceof OptimizerListener) {
				removeOptimizerListener((OptimizerListener) function);
			}
		}
	}

	
	private Matrix extractParameters(AbstractMathematicalFunction<Integer, Double, Integer, Double> function, List<Integer> indicesOfParametersToOptimize) {
		Matrix beta = new Matrix(indicesOfParametersToOptimize.size(),1);
		for (Integer index : indicesOfParametersToOptimize) {
			beta.m_afData[index][0] = function.getParameterValue(index);
		}
		return beta;
	}

	private void setParameters(AbstractMathematicalFunction<Integer, Double, Integer, Double> function, List<Integer> indicesOfParametersToOptimize, Matrix newBeta) {
		for (int i = 0; i < indicesOfParametersToOptimize.size(); i++) {
			int index = indicesOfParametersToOptimize.get(i);
			function.setParameterValue(index, newBeta.m_afData[i][0]);
		}
	}
	
	
	/**
	 * This method returns a double that is the convergence indicator based on the gradient, the hessian and the value of the log-likelihood function.
	 * @param gradient a Matrix instance
	 * @param hessian a Matrix instance
	 * @param llk a double that is the value of the log-likelihood function
	 * @return the indicator (a Double instance)
	 */
	protected double calculateConvergence(Matrix gradient, Matrix hessian, double llk) {
		return gradient.transpose().multiply(hessian.getInverseMatrix()).multiply(gradient).scalarMultiply(1d / llk).m_afData[0][0];
	}



	/**
	 * This method returns the number of iterations to reach convergence. It returns -1 if 
	 * convergence could not be reached.
	 * @return an integer
	 */
	public int getNumberOfIterations() {
		if (convergenceAchieved) {
			return iterationID;
		} else {
			return -1;
		}
	}

	
}
