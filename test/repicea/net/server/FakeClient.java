package repicea.net.server;

import java.net.InetSocketAddress;
import java.net.SocketAddress;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class FakeClient extends BasicClient {

	protected FakeClient(SocketAddress socketAddress) throws BasicClientException {
		super(socketAddress, 10);
	}

	
	@SuppressWarnings("unused")
	public static void main(String[] args) {
		Random random = new Random();
//		InetSocketAddress socketAddress = new InetSocketAddress(InetAddress.getLoopbackAddress(), 18000);
		InetSocketAddress socketAddress = new InetSocketAddress("rouge-epicea.dyndns.org", 18000);
		try {
			FakeClient client = new FakeClient(socketAddress);

			List<Double[]> locations = new ArrayList<Double[]>();
			Double[] location;
			for (int i = 0; i < 50; i++) {
				location = new Double[3];
				location[0] = 46d + random.nextDouble();
				location[1] = -71d + random.nextDouble();
				location[2] = 300d + 100d * random.nextDouble();
				locations.add(location);
			}
			
			Object result = client.processRequest(locations);
			client.close();
			int u = 0;
			System.out.println("Successfully treated.");
		} catch (BasicClientException e) {
			e.printStackTrace();
		}
	}
}
;