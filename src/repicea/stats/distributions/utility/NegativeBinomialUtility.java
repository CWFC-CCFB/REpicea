/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2017 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.distributions.utility;

import java.security.InvalidParameterException;

import repicea.math.GammaFunction;

public class NegativeBinomialUtility {
	
	/**
	 * This method returns the mass probability from a negative binomial distribution for a particular integer. It follows the 
	 * SAS parameterization.
	 * @see<a href=http://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_genmod_sect030.htm> 
	 * SAS online documentation </a>
	 * @param y the count (must be equal to or greater than 0)
	 * @param mean the mean of the distribution
	 * @param dispersion the dispersion parameter
	 * @return a mass probability
	 */
	public static double getMassProbability(int y, double mean, double dispersion) {
		if (y < 0) {
			throw new InvalidParameterException("The binomial distribution is designed for integer equals to or greater than 0!");
		}
		double prob = 0.0;
		double dispersionTimesMean = dispersion * mean;
		double inverseDispersion = 1/dispersion;

		prob = Math.exp(GammaFunction.logGamma(y + inverseDispersion) - GammaFunction.logGamma(y + 1.0) - GammaFunction.logGamma(inverseDispersion)) 
				*  Math.pow(dispersionTimesMean,y) / (Math.pow(1+dispersionTimesMean,y + inverseDispersion));
		return prob;
	}
}
