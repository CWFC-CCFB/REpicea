package repicea.net.server;

import java.util.ArrayList;

@SuppressWarnings({ "rawtypes", "serial" })
public class FakeArrayList extends ArrayList {

	
	public double processThisDouble(Double d) {
		return d.doubleValue();
	}
	
}
