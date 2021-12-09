package repicea.stats.estimators;

import repicea.stats.data.StatisticalDataStructure;
import repicea.stats.estimates.Estimate;
import repicea.stats.model.StatisticalModel;

public interface Estimator {

	/**
	 * The EstimatorException class encompasses all the exception that can be thrown when the
	 * optimizer fails to reach convergence.
	 * @author Mathieu Fortin - November 2015
	 */
	public static class EstimatorException extends Exception {
		private static final long serialVersionUID = 20110614L;

		public EstimatorException(String message) {
			super(message);
		}
	}
	
	public boolean doEstimation(StatisticalModel<? extends StatisticalDataStructure> model)	throws EstimatorException;

	/**
	 * This method returns true if the estimator successfully estimated the parameters.
	 * @return a boolean
	 */
	public boolean isConvergenceAchieved();

	/**
	 * This method returns the parameter estimates.
	 * @return an Estimate instance
	 */
	public Estimate<?> getParameterEstimates();

	/**
	 * Enable/disable verbose
	 * @param bool
	 */
	public void setVerboseEnabled(boolean bool);
	
}
