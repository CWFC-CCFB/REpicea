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
package repicea.stats.data;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * This class is the basic class for a hierarchical and spatial data structure.
 * @author Mathieu Fortin - October 2011
 */
public class GenericHierarchicalSpatialDataStructure extends GenericHierarchicalStatisticalDataStructure implements HierarchicalSpatialDataStructure {

	protected final List<Map<Integer, Map<Integer, Double>>> distanceMapList;
	protected List<List<String>> distanceFields;
	protected boolean distanceCalculated;
	protected final List<Double> distanceLimits;
	protected final Map<Integer, DistanceCalculator> distCalcMap;
//	protected Map<Integer, Map<Integer, Matrix>> angleMap;
//	protected boolean angleCalculated;
		
	/**
	 * General constructor. To be implemented in derived class.
	 * @param dataSet a DataSet instance
	 */
	public GenericHierarchicalSpatialDataStructure(DataSet dataSet) {
		super(dataSet);
		distanceMapList = new ArrayList<Map<Integer, Map<Integer, Double>>>();
//		angleMap = new HashMap<Integer, Map<Integer, Matrix>>();
		distanceFields = new ArrayList<List<String>>();
		distanceLimits = new ArrayList<Double>();
		distCalcMap = new HashMap<Integer, DistanceCalculator>();
}

	@Override
	public void setDistanceFields(List<List<String>> fields) {
		distanceLimits.clear();
		distanceFields = fields;
	}

	@Override
	public List<Map<Integer, Map<Integer, Double>>> getDistancesBetweenObservations() {
		if (distanceCalculated) {
			return distanceMapList;
 		} else {
 			if (isThereAnyHierarchicalStructure() && !distanceFields.isEmpty()) {
 				distanceMapList.clear();
 				int nbDistanceTypes = distanceFields.size();
 				for (int i = 0; i < nbDistanceTypes; i++) {
 					distanceMapList.add(new HashMap<Integer, Map<Integer,Double>>());
 				}
 				Map<String, DataBlock> hierarchicalStructure = getHierarchicalStructure();
 				for (String levelID : hierarchicalStructure.keySet()) {
 					List<Integer> observationIndex = hierarchicalStructure.get(levelID).getIndices();
 					if (observationIndex != null && !observationIndex.isEmpty()) {
 						int nbObs = observationIndex.size();
 						TreeMap<Integer, Double> tmpMap;
 						int indexA;
 						int indexB;
 						for (int i = 0; i < nbObs; i++) {
 							indexA = observationIndex.get(i);
 							outer:
 								for (int j = i + 1; j < nbObs; j++) {
 									indexB = observationIndex.get(j);
 									for (int dType = 0; dType < nbDistanceTypes; dType++) {
 										if (!distanceMapList.get(dType).containsKey(indexA)) {
 											distanceMapList.get(dType).put(indexA, new HashMap<Integer,Double>());
 										}
// 										tmpMap = new TreeMap<Integer, Double>();
 										double distance = getDistanceBetweenObservations(dType, indexA, indexB);
 										if (distance > getLimit(dType)) {
 											break outer;
 										}
// 										tmpMap.put(indexB, distance);
 										distanceMapList.get(dType).get(indexA).put(indexB, distance);
 									}
 								}
 						}
 					}
 				}
 				distanceCalculated = true;
 				return distanceMapList;
 			} else {
 				return null;
 			}
 		}
	}

	private double getLimit(int dType) {
		if (!distanceLimits.isEmpty()) {
			return distanceLimits.get(dType);
		} else {
			return Double.POSITIVE_INFINITY;
		}
	}
	
//	@Override
//	public Map<Integer, Map<Integer, Matrix>> getAngleBetweenObservations() {
//		if (angleCalculated) {
//			return angleMap;
// 		} else {
// 			if (isThereAnyHierarchicalStructure() && !distanceFields.isEmpty()) {
// 				angleMap.clear();
// 				Map<String, DataBlock> hierarchicalStructure = getHierarchicalStructure();
// 				for (String levelID : hierarchicalStructure.keySet()) {
// 					List<Integer> observationIndex = hierarchicalStructure.get(levelID).getIndices();
// 					if (observationIndex != null && !observationIndex.isEmpty()) {
// 						int nbObs = observationIndex.size();
// 						TreeMap<Integer, Matrix> tmpMap;
// 						int indexA;
// 						int indexB;
// 						for (int i = 0; i < nbObs; i++) {
// 							indexA = observationIndex.get(i);
// 							tmpMap = new TreeMap<Integer, Matrix>();
// 							angleMap.put(indexA, tmpMap);
// 							for (int j = 0; j < nbObs; j++) {
// 								indexB = observationIndex.get(j);
// 								Matrix angles = getAnglesBetweenObservations(indexA, indexB);
// 								tmpMap.put(indexB, angles);
// 							}
// 						}
// 					}
// 				}
// 				angleCalculated = true;
// 				return angleMap;
// 			} else {
// 				return null;
// 			}
// 		}
//	}

	
	
	
//	protected Matrix getAnglesBetweenObservations(int indexA, int indexB) {
//		int nbDimensions = distanceFields.size();
//		String fieldName;
//		double diff;
//		Object valueA;
//		Object valueB;
//
//		List<Double> listDiff = new ArrayList<Double>();
//		
//		for (int i = 0; i < nbDimensions; i++) {
//			fieldName = distanceFields.get(i);
//			int indexOfThisField = dataSet.getIndexOfThisField(fieldName);
//			valueA = dataSet.getValueAt(indexA, indexOfThisField);
//			valueB = dataSet.getValueAt(indexB, indexOfThisField);
//			diff = (Double) valueA - (Double) valueB;
//			listDiff.add(Math.abs(diff));
//		}
//
//		double angleRadians;
//		Matrix outputMatrix = new Matrix(nbDimensions,nbDimensions);
//		for (int i = 0; i < nbDimensions - 1; i++) {
//			for (int j = i + 1; j < nbDimensions; j++) {
//				angleRadians = Math.atan(listDiff.get(j) / listDiff.get(i));
//				outputMatrix.setValueAt(i, j, angleRadians);
//				outputMatrix.setValueAt(j, i, angleRadians);
//			}
//		}
//		
//		return outputMatrix;
//	}

	/**
	 * This method computes the Euclidian distance between two observations regardless of the dimensions. 
	 * @param indexA the index of the first observation
	 * @param indexB the index of the second observation
	 * @return a double the Euclidian distance
	 */
	protected double getDistanceBetweenObservations(int type, int indexA, int indexB) {
		if (type >= distanceFields.size()) {
			throw new InvalidParameterException("There are not enough distance types to compute the " + type + "th distance type!");
		}
		List<String> distanceTypeFields = distanceFields.get(type);
		return getDistanceCalculator(type).calculateDistance(distanceTypeFields, dataSet, indexA, indexB);
//		
//		
//		int nbDimensions = distanceType.size();
//		double ssDifferences = 0;
//		String fieldName;
//		double diff;
//		Object valueA;
//		Object valueB;
//		
//		for (int i = 0; i < nbDimensions; i++) {
//			fieldName = distanceType.get(i);
//			int indexOfThisField = dataSet.getIndexOfThisField(fieldName);
//			valueA = dataSet.getValueAt(indexA, indexOfThisField);
//			valueB = dataSet.getValueAt(indexB, indexOfThisField);
//			diff = ((Number) valueA).doubleValue() - ((Number) valueB).doubleValue();
//			ssDifferences += diff * diff;
//		}
//		
//		return Math.sqrt(ssDifferences);
	}

	@Override
	public void setDistanceLimits(List<Double> limits) {
		if (distanceFields.size() != limits.size()) {
			throw new InvalidParameterException("The number of limits in the limits argument is inconsistent: expected " + distanceFields.size() + " instead of " + limits.size());
		}
		this.distanceLimits.clear();
		this.distanceLimits.addAll(limits);
	}

	@Override
	public void setDistanceCalculators(DistanceCalculator... distanceCalculators) {	// could use a different calculator depending on the distance type
		if (distanceCalculators == null) {
			throw new InvalidParameterException("The distanceCalculators argument must be non null!");
		}
		if (distanceCalculators.length != distanceFields.size()) {
			throw new InvalidParameterException("The number of distance calculator instances is inconsistent! Please set the distance fields properly before setting the distance calculators!");
		}
		for (int i = 0; i < distanceCalculators.length; i++) {
			distCalcMap.put(i, distanceCalculators[i]);
		}
	}
	
	protected DistanceCalculator getDistanceCalculator(int distanceType) {
		if (distanceType >= distanceFields.size()) {
			throw new InvalidParameterException("The distanceType is inconsistent! Please set the distance fields properly before setting the distance calculators!");
		}
		if (!distCalcMap.containsKey(distanceType)) {
			distCalcMap.put(distanceType, new DistanceCalculator.EuclidianDistanceCalculator());
		}
		return distCalcMap.get(distanceType);
	}
	
}