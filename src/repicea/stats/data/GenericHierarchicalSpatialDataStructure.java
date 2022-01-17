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
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;

import repicea.stats.model.AbstractStatisticalModel;
import repicea.util.REpiceaLogManager;

/**
 * This class is the basic class for a hierarchical and spatial data structure.
 * @author Mathieu Fortin - October 2011
 */
public class GenericHierarchicalSpatialDataStructure extends GenericHierarchicalStatisticalDataStructure implements HierarchicalSpatialDataStructure {

	static class DistanceRecorderManager {
		final Map<Integer, DistanceRecorder> distanceRecorderMap;
		DistanceRecorderManager() {
			distanceRecorderMap = new HashMap<Integer, DistanceRecorder>();
		}
		
		void addIndices(List<Integer> indices) {
			DistanceRecorder dr = new DistanceRecorder(indices);
			for (Integer i : indices) {
				distanceRecorderMap.put(i, dr);
			}
		}
	}
	
	static class DistanceRecorder {
		final List<Integer> index;
		final int nbObs;
		final int offset;
		final double[] distances;
		
		DistanceRecorder(List<Integer> index) {
			this.index = new ArrayList<Integer>();
			this.index.addAll(index);
			Collections.sort(this.index);
			offset = this.index.get(0);
			nbObs = index.size();
			int arraySize = ((nbObs - 1) * nbObs) / 2;
			distances = new double[arraySize];
		}
		
		void setValueAt(int indexA, int indexB, double value) {
			int index = getIndex(indexA, indexB);
			distances[index] = value;
		}
		
		double getValueAt(int indexA, int indexB) {
			int index = getIndex(indexA, indexB);
			return distances[index];
		}
		
		private int getIndex(int indexA, int indexB) {
			int i = indexA - offset;
			int j = indexB - offset;
			if (j > i) {
				int rowOffset = i * (nbObs - 1) - (int) (((i - 1) * i) * .5);
				return rowOffset + (j - i - 1);
			} else {
				throw new InvalidParameterException("Index A must refer to an index smaller than that of index B!");
			}
		}
	}
	
	
	protected final List<DistanceRecorderManager> distanceRecorderManagers;
	protected List<List<String>> distanceFields;
	protected boolean distanceCalculated;
	protected final Map<Integer, DistanceCalculator> distCalcMap;
		
	/**
	 * General constructor. To be implemented in derived class.
	 * @param dataSet a DataSet instance
	 */
	public GenericHierarchicalSpatialDataStructure(DataSet dataSet) {
		super(dataSet);		// the data set is sorted according to the structure
		distanceRecorderManagers = new ArrayList<DistanceRecorderManager>();
		distanceFields = new ArrayList<List<String>>();
		distCalcMap = new HashMap<Integer, DistanceCalculator>();
}

	@Override
	public void setDistanceFields(List<List<String>> fields) {
		distanceFields = fields;
	}

	@Override
	public double getDistancesBetweenObservations(int type, int ii, int jj) {
		if (!distanceCalculated) {
			if (distanceFields.isEmpty()) {
				throw new InvalidParameterException("There is no available field to calculate the distances!");
			}
 			System.out.println("Calculating distances between the observations...");
 			if (isThereAnyHierarchicalStructure()) {	// then we should use the maps in distanceMapList 
 				distanceRecorderManagers.clear();
 				int nbDistanceTypes = distanceFields.size();
 				for (int i = 0; i < nbDistanceTypes; i++) {
 					distanceRecorderManagers.add(new DistanceRecorderManager());
 				}
 				Map<String, DataBlock> hierarchicalStructure = getHierarchicalStructure();
 				for (String levelID : hierarchicalStructure.keySet()) {
 					List<Integer> observationIndex = hierarchicalStructure.get(levelID).getIndices();
 					if (observationIndex != null && !observationIndex.isEmpty()) {
						for (int dType = 0; dType < nbDistanceTypes; dType++) {
							distanceRecorderManagers.get(dType).addIndices(observationIndex);
						}
						int nbObs = observationIndex.size();
 						int indexA, indexB;
 						for (int i = 0; i < nbObs; i++) {
 							indexA = observationIndex.get(i);
 	 						if (indexA%1000 == 0) {
 	 							REpiceaLogManager.logMessage(AbstractStatisticalModel.LoggerName,
 	 									Level.FINE, 
 	 									null,
 	 									"Processing observation " + indexA);
 	 						}
 							for (int j = i + 1; j < nbObs; j++) {
 								indexB = observationIndex.get(j);
 								for (int dType = 0; dType < nbDistanceTypes; dType++) {
 									double distance = calculateDistanceBetweenObservations(dType, indexA, indexB);
 									distanceRecorderManagers.get(dType).distanceRecorderMap.get(indexA).setValueAt(indexA, indexB, distance);
 								}
 							}
 						}
 					}
 				}
 			} else {
 				throw new InvalidParameterException("There is no hierarchical structure in the data!");
 			}
 			distanceCalculated = true;
			System.out.println("Distances have been calculated.");
		}				
 		
		return distanceRecorderManagers.get(type).distanceRecorderMap.get(ii).getValueAt(ii, jj);
	}

	
	/**
	 * This method computes the Euclidian distance between two observations regardless of the dimensions. 
	 * @param indexA the index of the first observation
	 * @param indexB the index of the second observation
	 * @return a double the Euclidian distance
	 */
	protected double calculateDistanceBetweenObservations(int type, int indexA, int indexB) {
		if (type >= distanceFields.size()) {
			throw new InvalidParameterException("There are not enough distance types to compute the " + type + "th distance type!");
		}
		List<String> distanceTypeFields = distanceFields.get(type);
		return getDistanceCalculator(type).calculateDistance(distanceTypeFields, dataSet, indexA, indexB);
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
