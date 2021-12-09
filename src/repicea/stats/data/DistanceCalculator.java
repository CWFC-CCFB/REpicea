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

/**
 * An abstract class for all distance calculator.
 * @author Mathieu Fortin - December 2021
 */
public abstract class DistanceCalculator {
	
	/**
	 * Calculate the distance between two observations. <br>
	 * <br>
	 * The distance is calculated using the fields represented by the fields argument.
	 * @param fields a List of strings
	 * @param dataSet the data set that contains the data
	 * @param indexA the index of the first observation
	 * @param indexB the index of the second observation
	 * @return a double
	 */
	public abstract double calculateDistance(List<String> fields, DataSet dataSet, int indexA, int indexB); 
	
	/**
	 * A distance calculator based on the Euclidian distance.
	 * @author Mathieu Fortin - December 2021
	 */
	public static class EuclidianDistanceCalculator extends DistanceCalculator {
		
		@Override
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

	/**
	 * A distance calculator based on geographical coordinates.
	 * @author Mathieu Fortin - December 2021
	 */
	public static class GeographicDistanceCalculator extends DistanceCalculator {
		
		private static final double EarthCircumferenceKmEquator = 40075.0167; 
		private static final double EarthCircumferenceKmPole = 40007.863; 
		private static final double OneDividedBy360 = 1d/360;

		private double getCircumferenceAtThisLat(double latitude) {
		  double latRad = latitude * OneDividedBy360 * 2 * Math.PI;
		  return Math.cos(latRad) * EarthCircumferenceKmEquator;
		} 

		private double getKmFromLatitudeDiff(double latDegDiff) {
			return latDegDiff * EarthCircumferenceKmPole * OneDividedBy360;
		}

		/**
		 * Calculate the distance in km between two geographic coordinates.
		 * @param fields a List of strings IMPORTANT the longitude field comes first and then the latitude field
		 * @param dataSet the DataSet instance that contains the data
		 * @param indexA the index of the first observation
		 * @param indexB the index of the second observation
		 * @return the distance (km)
		 */
		@Override
		public double calculateDistance(List<String> fields, DataSet dataSet, int indexA, int indexB) {
			String longitudeFieldName = fields.get(0);
			String latitudeFieldName = fields.get(1);
			
			double longitudeA = ((Number) dataSet.getValueAt(indexA, longitudeFieldName)).doubleValue();
			double longitudeB = ((Number) dataSet.getValueAt(indexB, longitudeFieldName)).doubleValue();
			double latitudeA = ((Number) dataSet.getValueAt(indexA, latitudeFieldName)).doubleValue();
			double latitudeB = ((Number) dataSet.getValueAt(indexB, latitudeFieldName)).doubleValue();
			
			double meanLatitude = (latitudeB + latitudeA) * .5;
			
			double diffX = (longitudeA - longitudeB) * getCircumferenceAtThisLat(meanLatitude) * OneDividedBy360;
			double diffY = getKmFromLatitudeDiff(latitudeA - latitudeB);
			return Math.sqrt(diffX * diffX + diffY * diffY) / 100;
		}
		
	}

}
