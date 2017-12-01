/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2016 Mathieu Fortin for Rouge Epicea.
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
package repicea.net.server;

import java.io.Closeable;
import java.io.IOException;
import java.net.Socket;
import java.net.SocketAddress;
import java.net.UnknownHostException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.TimeoutException;

import repicea.net.SocketWrapper;
import repicea.net.server.AbstractServer.ServerReply;
import repicea.util.REpiceaTranslator;
import repicea.util.REpiceaTranslator.TextableEnum;

public class BasicClient implements Closeable {
	
	protected static enum ClientRequest {closeConnection}

	public static enum ExceptionType implements TextableEnum {
		ConnectionFailed("The client failed to connect to the server.", "La connexion au serveur n'a pas pu \u00EAtre \u00E9tablie."),
		ConnectionRejected("The server is busy and cannot process the requests.", "Le serveur est actuellement occup\u00E9 et ne peut r\u00E9pondre aux requ\u00EAtes."),
		ConnectionTooLong("The reply from the server has exceeded the allowed time.", "La r\u00E9ponse du serveur a ex\u00E9c\u00E9d\u00E9 le temps d'attente."),
		UnknownErrorWhileConnected("An exception occurred while processing the request.", "Une exception est survenue dans l'ex\u00E9cution de la requ\u00EAte.");
		
		ExceptionType(String englishText, String frenchText) {
			setText(englishText, frenchText);
		}
		
		@Override
		public void setText(String englishText, String frenchText) {
			REpiceaTranslator.setString(this, englishText, frenchText);
		}

		@Override
		public String toString() {
			return REpiceaTranslator.getString(this);
		}
	}
	
	
	@SuppressWarnings("serial")
	public static class BasicClientException extends Exception {
		private BasicClientException(ExceptionType exceptionType) {
			super(exceptionType.toString());
		}
	}
	
	
	
	private SocketWrapper socketWrapper;
	private boolean open;
	private final int timeout;
	
	
	/**
	 * Constructor.
	 * @param socketAddress the SocketAddress instance that corresponds to the server
	 * @param timeoutSeconds the number of seconds to wait for server reply before throwing a TimeoutException
	 * @throws BasicClientException 
	 * @throws UnknownHostException if the host is unknown
	 * @throws IOException if the connection failed
	 * @throws ClassNotFoundException if the reply from the server is incorrect
	 * @throws ClientException if the connection has been lost
	 */
	@SuppressWarnings("resource")
	protected BasicClient(SocketAddress socketAddress, int timeoutSeconds) throws BasicClientException {
		this.timeout = timeoutSeconds;
		Socket socket = new Socket();
		try {
			socket.connect(socketAddress, 5000);
		} catch (IOException e) {
			close();
			throw new BasicClientException(ExceptionType.ConnectionFailed);
		} 
		
		socketWrapper = new SocketWrapper(socket);

		ServerReply replyFromServer = (ServerReply) readObjectFromServer();
		if (replyFromServer == ServerReply.CallAccepted) {
			open = true;
		} else if (replyFromServer == ServerReply.IAmBusyCallBackLater) {
			open = false;
			close();
			throw new BasicClientException(ExceptionType.ConnectionRejected);
		}
		
	}

	
	protected Object readObjectFromServer() throws BasicClientException {
		try {
			return socketWrapper.readObject(timeout);
		} catch (Exception e) {
			e.printStackTrace();
			close();
			throw handleException(e);
		} 
	}
	
	
	
	private BasicClientException handleException(Exception e) {
		if (e instanceof InterruptedException || e instanceof ExecutionException || e instanceof TimeoutException) {
			return new BasicClientException(ExceptionType.ConnectionTooLong);
		} else {
			return new BasicClientException(ExceptionType.UnknownErrorWhileConnected);
		}
	}


	protected Object processRequest(Object obj) throws BasicClientException {
		if (open) {
			try {
				socketWrapper.writeObject(obj);
			} catch (IOException e) {
				close();
				e.printStackTrace();
			}
			return readObjectFromServer();
		}
		throw new BasicClientException(ExceptionType.ConnectionFailed);
	}
	
	/**
	 * This method sends a request to the server to close the connection and closes the output and input streams.
	 */
	@Override
	public void close() {
		try {
			if (open) {
				socketWrapper.writeObject(ClientRequest.closeConnection);
			}
		} catch (IOException e) {
		} finally {
			try {
				socketWrapper.close();
			} catch (IOException e) {}
		}
	}

	
}