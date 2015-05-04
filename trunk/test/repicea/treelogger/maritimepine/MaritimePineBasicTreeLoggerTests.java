package repicea.treelogger.maritimepine;

import java.util.ArrayList;
import java.util.Collection;

import org.junit.Assert;
import org.junit.Test;

import repicea.simulation.treelogger.WoodPiece;

public class MaritimePineBasicTreeLoggerTests {
	
	@SuppressWarnings({ "rawtypes", "unchecked" })
	@Test
	public void TestWithSimpleTree() {
		MaritimePineBasicTreeLogger treeLogger = new MaritimePineBasicTreeLogger();
		treeLogger.setTreeLoggerParameters(treeLogger.createDefaultTreeLoggerParameters());
		Collection trees = new ArrayList<MaritimePineTree>();
		MaritimePineTree tree = new MaritimePineTree();
		trees.add(tree);
		treeLogger.init(trees);
		treeLogger.run();
		double sum = 0;
		for (WoodPiece piece : treeLogger.getWoodPieces().get(tree)) {
			double volumeM3 = piece.getVolumeM3();
			sum += volumeM3;
		}
		Assert.assertEquals("Comparing bole volume", 1d, sum, 1E-8); 
	}
}
