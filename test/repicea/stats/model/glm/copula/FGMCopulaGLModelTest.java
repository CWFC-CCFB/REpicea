package repicea.stats.model.glm.copula;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import org.junit.Ignore;
import org.junit.Test;

import repicea.math.ParameterBound;
import repicea.stats.data.DataSet;
import repicea.stats.data.DistanceCalculator.GeographicDistanceCalculator;
import repicea.stats.data.StatisticalDataException;
import repicea.stats.model.glm.GeneralizedLinearModel;
import repicea.stats.model.glm.LinkFunction.Type;
import repicea.stats.model.glm.copula.CopulaLibrary.DistanceLinkFunctionCopulaExpression;
import repicea.util.ObjectUtility;

public class FGMCopulaGLModelTest {

	/**
	 * This test implements a constant parameter copula model.
	 * @throws Exception
	 */
//	@Ignore
	@Test
    public void TestWithSimpleCopula() throws Exception {
		double expectedCopulaValue = 0.18072406716297676;
		double expectedLlk = -1077.162324348724;
		String filename = ObjectUtility.getPackagePath(FGMCopulaGLModelTest.class).concat("donneesR_min.csv");
		DataSet dataSet = new DataSet(filename, true);
		
		GeneralizedLinearModel glm = new GeneralizedLinearModel(dataSet, Type.Logit, "coupe ~ diffdhp + marchand:diffdhp + marchand:diffdhp2 +  essence");
		
		try {
			CopulaExpression copula = new CopulaLibrary.SimpleCopulaExpression(0.0, "IDENT");
			FGMCopulaGLModel copulaModel = new FGMCopulaGLModel(glm, copula);
			copulaModel.doEstimation();
//			copulaModel.getSummary();
			double actual = copula.getValue();
			assertEquals(expectedCopulaValue, actual, 1E-8);
			double actualLlk = copulaModel.getCompleteLogLikelihood().getValue();
			assertEquals(expectedLlk, actualLlk, 1E-8);
		} catch (StatisticalDataException e) {
			e.printStackTrace();
			throw e;
		}
		
	}

	/**
	 * This test implements a constant parameter copula model.
	 * @throws Exception
	 */
//	@Ignore
	@Test
    public void TestWithDistanceCopula() throws Exception {
		double expectedCopulaValue = -0.17036586673288104;
		double expectedLlk = -1062.1351502297755;
		String filename = ObjectUtility.getPackagePath(FGMCopulaGLModelTest.class).concat("donneesR_min.csv");
		DataSet dataSet = new DataSet(filename, true);
		
		GeneralizedLinearModel glm = new GeneralizedLinearModel(dataSet, Type.Logit, "coupe ~ diffdhp + marchand:diffdhp + marchand:diffdhp2 +  essence");
		
		try {
			CopulaExpression distanceCopula = new CopulaLibrary.DistanceLinkFunctionCopulaExpression(Type.Log, 
					"IDENT", 
					"X + Y", 
					false, // no intercept
					true, // strictly positive
					null,
					-.15);		
			((DistanceLinkFunctionCopulaExpression) distanceCopula).setBounds(0, new ParameterBound(null, 0d));
			FGMCopulaGLModel copulaModel = new FGMCopulaGLModel(glm, distanceCopula);
			copulaModel.setConvergenceCriterion(1E-8);
			copulaModel.gridSearch(copulaModel.getParameters().m_iRows - 1, -.25d, -.15d, .01);
			copulaModel.doEstimation();
//			copulaModel.getSummary();
			double actual = distanceCopula.getBeta().getValueAt(0, 0);
			assertEquals(expectedCopulaValue, actual, 1E-8);
			double actualLlk = copulaModel.getCompleteLogLikelihood().getValue();
			assertEquals(expectedLlk, actualLlk, 1E-8);
		} catch (StatisticalDataException e) {
			e.printStackTrace();
			throw e;
		}
	}
	
	
	// SIMPLE COPULA FOR NEWID_PE			
//	@Ignore
	@Test
    public void TestWithSimpleCopulaWithBeechData() throws Exception {
		double expectedCopulaValue = 0.5103720998957171;
		double expectedLlk = -2722.197862205543;
		String filename = ObjectUtility.getPackagePath(FGMCopulaGLModelTest.class).concat("copulaHEG.csv");
		DataSet dataSet = new DataSet(filename, true);
		GeneralizedLinearModel glm = new GeneralizedLinearModel(dataSet, Type.CLogLog, "occurred ~ lnDt + logDD + TotalPrcp + logPrcp + G_TOT + lnG_TOT + speciesThere + lnG_SpGr");
		glm.doEstimation();
		glm.getSummary();
		try {
			CopulaExpression simpleCopula = new CopulaLibrary.SimpleCopulaExpression(0, "newID_PE");
			FGMCopulaGLModel copulaModel = new FGMCopulaGLModel(glm, simpleCopula);
			copulaModel.setConvergenceCriterion(1E-8);
			copulaModel.gridSearch(copulaModel.getParameters().m_iRows - 1, 0d, .9d, .1);
			copulaModel.doEstimation();
//			copulaModel.getSummary();
			double actualCopulaValue = simpleCopula.getBeta().getValueAt(0, 0);
			assertEquals(expectedCopulaValue, actualCopulaValue, 1E-8);
			double actualLlk = copulaModel.getCompleteLogLikelihood().getValue();
			assertEquals(expectedLlk, actualLlk, 1E-8);
	
		} catch (Exception e) {
			e.printStackTrace();
			throw e;
		}
	}


	// DISTANCE COPULA (temporal) FOR NEWID_PE			
//	@Ignore
	@Test
    public void TestWithTemporalCorrelationOnBeechData() throws Exception {
		double expectedCopulaValue = -0.04040906569629807;
		double expectedLlk = -2717.8099755938074;
		String filename = ObjectUtility.getPackagePath(FGMCopulaGLModelTest.class).concat("copulaHEG.csv");
		DataSet dataSet = new DataSet(filename, true);
		GeneralizedLinearModel glm = new GeneralizedLinearModel(dataSet, Type.CLogLog, "occurred ~ lnDt + logDD + TotalPrcp + logPrcp + G_TOT + lnG_TOT + speciesThere + lnG_SpGr");
		glm.doEstimation();
		glm.getSummary();
		try {
			CopulaExpression distanceCopula2 = new CopulaLibrary.DistanceLinkFunctionCopulaExpression(Type.Log, 
					"newID_PE", 
					"year.y", 
					false, // no intercept
					true, // strictly positive
					null, 
					-.15);		
			((DistanceLinkFunctionCopulaExpression) distanceCopula2).setBounds(0, new ParameterBound(null, 0d));
			FGMCopulaGLModel copulaModel2 = new FGMCopulaGLModel(glm, distanceCopula2);
			copulaModel2.setConvergenceCriterion(1E-8);
			copulaModel2.gridSearch(copulaModel2.getParameters().m_iRows - 1, -.5d, -.1d, .05);
			copulaModel2.doEstimation();
//			copulaModel2.getSummary();
			double actual = distanceCopula2.getBeta().getValueAt(0, 0);
			double actualLlk = copulaModel2.getCompleteLogLikelihood().getValue();
			assertEquals(expectedCopulaValue, actual, 1E-8);
			assertEquals(expectedLlk, actualLlk, 1E-8);
		} catch (Exception e) {
			e.printStackTrace();
			throw e;
		}
	}

	
	@Ignore
	@Test
    public void TestWithSimpleCopulaHEG() throws Exception {
		String filename = ObjectUtility.getPackagePath(FGMCopulaGLModelTest.class).concat("copulaHEG.csv");
//		String filename = ObjectUtility.getPackagePath(FGMCopulaGLModelTest.class).concat("copulaHEGOneMesPerPlot.csv");
		DataSet dataSet = new DataSet(filename, true);
		GeneralizedLinearModel glm = new GeneralizedLinearModel(dataSet, Type.CLogLog, "occurred ~ lnDt + logDD + TotalPrcp + logPrcp + G_TOT + lnG_TOT + speciesThere + lnG_SpGr");
		glm.doEstimation();
		glm.getSummary();
		try {
// SIMPLE COPULA FOR NEWID_PE			
			CopulaExpression simpleCopula = new CopulaLibrary.SimpleCopulaExpression(0, "newID_PE");
			FGMCopulaGLModel copulaModel = new FGMCopulaGLModel(glm, simpleCopula);
			copulaModel.setConvergenceCriterion(1E-8);
			copulaModel.gridSearch(copulaModel.getParameters().m_iRows - 1, 0d, .9d, .1);
			copulaModel.doEstimation();
			copulaModel.getSummary();

// DISTANCE COPULA (temporal) FOR NEWID_PE			
			CopulaExpression distanceCopula2 = new CopulaLibrary.DistanceLinkFunctionCopulaExpression(Type.Log, 
					"newID_PE", 
					"year.y", 
					false, // no intercept
					true, // strictly positive
					null, 
					-.15);		
			((DistanceLinkFunctionCopulaExpression) distanceCopula2).setBounds(0, new ParameterBound(null, 0d));
			FGMCopulaGLModel copulaModel2 = new FGMCopulaGLModel(glm, distanceCopula2);
			copulaModel2.setConvergenceCriterion(1E-8);
			copulaModel2.gridSearch(copulaModel.getParameters().m_iRows - 1, -.5d, -.1d, .05);
			copulaModel2.doEstimation();
			copulaModel2.getSummary();

//// DISTANCE COPULA (spatial)
//			CopulaExpression distanceCopula = new CopulaLibrary.DistanceLinkFunctionCopulaExpression(Type.Log, 
//					"whole", 
//					"longitudeDeg + latitudeDeg", 
//					false,
//					null, 
//					Arrays.asList(new GeographicDistanceCalculator()), -.6);		// no intercept in the linear expression
//			((DistanceLinkFunctionCopulaExpression) distanceCopula).setBounds(0, new ParameterBound(null, 0d));
////			((DistanceLinkFunctionCopulaExpression) distanceCopula).setBounds(1, new ParameterBound(null, 0d));
//			FGMCopulaGLModel copulaModel = new FGMCopulaGLModel(glm, distanceCopula);
//			copulaModel.getEstimator().setVerboseEnabled(true);;
//			copulaModel.setConvergenceCriterion(1E-8);
////			copulaModel.gridSearch(copulaModel.getParameters().m_iRows - 1, -.4d, -.1d, .1);
//			copulaModel.doEstimation();
//			copulaModel.getSummary();

//// DISTANCE COPULA (temporal + spatial)
//			CopulaExpression distanceCopula = new CopulaLibrary.DistanceLinkFunctionCopulaExpression(Type.Log, "whole", "longitudeDeg + latitudeDeg, year.y", Arrays.asList(new Double[]{1d, 100d}), -.05, -.04);		// no intercept in the linear expression
//			((DistanceLinkFunctionCopulaExpression) distanceCopula).setBounds(0, new ParameterBound(null, 0d));
//			((DistanceLinkFunctionCopulaExpression) distanceCopula).setBounds(1, new ParameterBound(null, 0d));
//			FGMCopulaGLModel copulaModel = new FGMCopulaGLModel(glm, distanceCopula);
//			copulaModel.setConvergenceCriterion(1E-8);
////			copulaModel.gridSearch(copulaModel.getParameters().m_iRows - 1, -.5d, -.1d, .1);
//			copulaModel.doEstimation();
//			copulaModel.getSummary();

			
		} catch (Exception e) {
			e.printStackTrace();
			throw e;
		}
	}

	
	
	
	
}
