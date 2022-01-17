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

import java.util.Arrays;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import repicea.stats.data.GenericHierarchicalSpatialDataStructure.DistanceRecorder;

public class GenericHierarchicalSpatialDataStructureTest {

	@Test
	public void testDistanceRecording() {
		List<Integer> indices = Arrays.asList(new Integer[]{8, 9, 10, 11, 12, 13});
		DistanceRecorder dr = new DistanceRecorder(indices);
		Assert.assertEquals("Testing number of observations", 6, dr.nbObs);
		Assert.assertEquals("Testing array size", 6 * 5 / 2, dr.distances.length);
		double value = 1;
		for (int i = 0; i < indices.size(); i++) {
			int indexA = indices.get(i);
			for (int j = i + 1; j < indices.size(); j++) {
				int indexB = indices.get(j);
				dr.setValueAt(indexA, indexB, value++);
			}
		}
		
		double[] expectedArray = new double[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
		Assert.assertArrayEquals("Testing values in the array", expectedArray, dr.distances, 1E-12);
		
		double distance11to13 = dr.getValueAt(11, 13);
		Assert.assertEquals("Testing distance between 11 and 13", 14, distance11to13, 1E-12);
	}
}
