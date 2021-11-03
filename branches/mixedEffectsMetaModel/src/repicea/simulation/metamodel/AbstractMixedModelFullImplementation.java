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

import java.util.List;

import repicea.math.Matrix;
import repicea.stats.data.HierarchicalStatisticalDataStructure;
import repicea.stats.data.StatisticalDataException;
import repicea.stats.distributions.GaussianDistribution;

/**
 * The AbstractMixedModelImplementation class relies on a marginal log-likelihood instead of a pure 
 * log-likelihood.
 * @author Mathieu Fortin - October 2021
 */
abstract class AbstractMixedModelFullImplementation extends AbstractModelImplementation {

	int indexRandomEffectVariance;
	int indexFirstRandomEffect;
	
	GaussianDistribution randomEffectDistribution;
	
	AbstractMixedModelFullImplementation(String outputType, MetaModel metaModel) throws StatisticalDataException {
		super(outputType, metaModel);
	}

	@Override
	protected final AbstractDataBlockWrapper createWrapper(String k, List<Integer> indices, HierarchicalStatisticalDataStructure structure, Matrix varCov) {
		return new DataBlockWrapper(k, indices, structure, varCov);
	}

	@Override
	public final double getLogLikelihood(Matrix parameters) {
		setParameters(parameters);
		double logLikelihood = 0d;
		for (int i = 0; i < dataBlockWrappers.size(); i++) {
			AbstractDataBlockWrapper dbw = dataBlockWrappers.get(i);
			dbw.setParameterValue(0, parameters.getValueAt(indexFirstRandomEffect + i, 0));
			double logLikelihoodForThisBlock = dbw.getLogLikelihood();
			logLikelihood += logLikelihoodForThisBlock;
		}
		return logLikelihood;
	}


}
