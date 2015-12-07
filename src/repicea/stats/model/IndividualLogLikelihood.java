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
package repicea.stats.model;

import repicea.math.LogFunctionWrapper;


/**
 * The LogLikelihood interface provides the basic services for all LogLikelihood classes. A 
 * LogLLikelihood instance must be able to provide its original likelihood function and its 
 * derivatives.
 * @author Mathieu Fortin - June 2011
 */
@SuppressWarnings("serial")
public class IndividualLogLikelihood extends LogFunctionWrapper implements LikelihoodCompatible<Double> {
	
	public IndividualLogLikelihood(IndividualLikelihood originalFunction) {
		super(originalFunction);
	}
	
	@Override
	public IndividualLikelihood getOriginalFunction() {return (IndividualLikelihood) super.getOriginalFunction();}
	
	
	public double getPrediction() {return getOriginalFunction().getPrediction();}

	@Override
	public void setY(Double y) {getOriginalFunction().setY(y);}
	
}
