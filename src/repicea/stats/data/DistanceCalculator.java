/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2021 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.data;

import java.util.List;

public abstract class DistanceCalculator {

	public abstract double calculateDistance(List<String> fields, DataSet dataSet, int indexA, int indexB); 
	
	public static class EuclidianDistanceCalculator extends DistanceCalculator {
		
		public double calculateDistance(List<String> fields, DataSet dataSet, int indexA, int indexB) {
			int nbDimensions = fields.size();
			double ssDifferences = 0;
			String fieldName;
			double diff;
			Object valueA;
			Object valueB;
			
			for (int i = 0; i < nbDimensions; i++) {
				fieldName = fields.get(i);
				int indexOfThisField = dataSet.getIndexOfThisField(fieldName);
				valueA = dataSet.getValueAt(indexA, indexOfThisField);
				valueB = dataSet.getValueAt(indexB, indexOfThisField);
				diff = ((Number) valueA).doubleValue() - ((Number) valueB).doubleValue();
				ssDifferences += diff * diff;
			}
			return Math.sqrt(ssDifferences);
		}
		
	}
	
}
