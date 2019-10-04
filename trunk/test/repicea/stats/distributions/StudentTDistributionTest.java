package repicea.stats.distributions;


import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.stats.estimates.ConfidenceInterval;
import repicea.stats.estimates.MonteCarloEstimate;

public class StudentTDistributionTest {

	
	@Test
	public void testRandomNumberGenerationWith3DegreesOfFreedom() {
		
		MonteCarloEstimate estimate = new MonteCarloEstimate();
		StudentTDistribution dist = new StudentTDistribution(3);
		for (int i = 0; i < 5000000; i++) {
			estimate.addRealization(dist.getRandomRealization());
		}
		double mean = estimate.getMean().m_afData[0][0];
		double variance = estimate.getVariance().m_afData[0][0];
		ConfidenceInterval ci = estimate.getConfidenceIntervalBounds(.95);
		double quantile025 = ci.getLowerLimit().m_afData[0][0];
		double quantile975 = ci.getUpperLimit().m_afData[0][0];
		double expectedMean = dist.getMean().m_afData[0][0];
		double expectedVariance = dist.getVariance().m_afData[0][0];
		Assert.assertEquals("Testing the mean", expectedMean, mean, 2E-3);
		Assert.assertEquals("Testing the variance", expectedVariance, variance, 5E-2);
		Assert.assertEquals("Testing quantile 0.025", quantile025, -3.182446, 5E-2);
		Assert.assertEquals("Testing quantile 0.975", quantile975, 3.182446, 5E-2);
	}

	@Test
	public void testRandomNumberGenerationWith10DegreesOfFreedom() {
		
		MonteCarloEstimate estimate = new MonteCarloEstimate();
		StudentTDistribution dist = new StudentTDistribution(10);
		for (int i = 0; i < 5000000; i++) {
			estimate.addRealization(dist.getRandomRealization());
		}
		double mean = estimate.getMean().m_afData[0][0];
		double variance = estimate.getVariance().m_afData[0][0];
		ConfidenceInterval ci = estimate.getConfidenceIntervalBounds(.95);
		double quantile025 = ci.getLowerLimit().m_afData[0][0];
		double quantile975 = ci.getUpperLimit().m_afData[0][0];
		double expectedMean = dist.getMean().m_afData[0][0];
		double expectedVariance = dist.getVariance().m_afData[0][0];
		Assert.assertEquals("Testing the mean", expectedMean, mean, 2E-3);
		Assert.assertEquals("Testing the variance", expectedVariance, variance, 5E-2);
		Assert.assertEquals("Testing quantile 0.025", quantile025, -2.228139, 5E-2);
		Assert.assertEquals("Testing quantile 0.975", quantile975, 2.228139, 5E-2);
	}

	@Test
	public void testRandomNumberGenerationWith20DegreesOfFreedom() {
		
		MonteCarloEstimate estimate = new MonteCarloEstimate();
		StudentTDistribution dist = new StudentTDistribution(20);
		for (int i = 0; i < 5000000; i++) {
			estimate.addRealization(dist.getRandomRealization());
		}
		double mean = estimate.getMean().m_afData[0][0];
		double variance = estimate.getVariance().m_afData[0][0];
		ConfidenceInterval ci = estimate.getConfidenceIntervalBounds(.95);
		double quantile025 = ci.getLowerLimit().m_afData[0][0];
		double quantile975 = ci.getUpperLimit().m_afData[0][0];
		double expectedMean = dist.getMean().m_afData[0][0];
		double expectedVariance = dist.getVariance().m_afData[0][0];
		Assert.assertEquals("Testing the mean", expectedMean, mean, 2E-3);
		Assert.assertEquals("Testing the variance", expectedVariance, variance, 5E-2);
		Assert.assertEquals("Testing quantile 0.025", quantile025, -2.085963, 5E-2);
		Assert.assertEquals("Testing quantile 0.975", quantile975, 2.085963, 5E-2);
	}
	

	@Test
	public void testRandomNumberGenerationWithNonCentered10DegreesOfFreedom() {
		
		MonteCarloEstimate estimate = new MonteCarloEstimate();
		Matrix mean = new Matrix(1,1);
		mean.m_afData[0][0] = 10;
		Matrix variance = new Matrix(1,1);
		variance.m_afData[0][0] = 20;
		StudentTDistribution dist = new StudentTDistribution(mean, variance, 10);
		for (int i = 0; i < 5000000; i++) {
			estimate.addRealization(dist.getRandomRealization());
		}
		double actualMean = estimate.getMean().m_afData[0][0];
		double actualVariance = estimate.getVariance().m_afData[0][0];
		ConfidenceInterval ci = estimate.getConfidenceIntervalBounds(.95);
		double quantile025 = ci.getLowerLimit().m_afData[0][0];
		double quantile975 = ci.getUpperLimit().m_afData[0][0];
		double expectedMean = dist.getMean().m_afData[0][0];
		double expectedVariance = dist.getVariance().m_afData[0][0];
		double expectedQuantile025 = -2.228139 * Math.sqrt(20) + 10;
		double expectedQuantile975 = 2.228139 * Math.sqrt(20) + 10;
		Assert.assertEquals("Testing the mean", expectedMean, actualMean, 2E-3);
		Assert.assertEquals("Testing the variance", expectedVariance, actualVariance, 5E-2);
		Assert.assertEquals("Testing quantile 0.025", quantile025, expectedQuantile025, 5E-2);
		Assert.assertEquals("Testing quantile 0.975", quantile975, expectedQuantile975, 5E-2);
	}


}