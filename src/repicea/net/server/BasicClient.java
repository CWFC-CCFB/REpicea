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
import java.net.InetSocketAddress;
import java.net.Socket;
import java.net.SocketTimeoutException;
import java.security.InvalidParameterException;

import repicea.net.SocketWrapper;
import repicea.net.TCPSocketWrapper;
import repicea.net.UDPSocketWrapper;
import repicea.net.server.AbstractServer.ServerReply;
import repicea.net.server.ServerConfiguration.Protocol;
import repicea.util.REpiceaTranslator;
import repicea.util.REpiceaTranslator.TextableEnum;

public class BasicClient implements Closeable {
	
	public static enum ClientRequest {closeConnection}

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
	private final Protocol protocol;
	private boolean open;
	private final int timeout;
	private boolean bypassTimeout;

	
	/**
	 * Constructor.
	 * @param socketAddress the SocketAddress instance that corresponds to the server
	 * @param timeoutSeconds the number of seconds to wait for server reply before throwing a TimeoutException
	 * @throws BasicClientException 
	 */
	protected BasicClient(InetSocketAddress socketAddress, int timeoutSeconds) throws BasicClientException {
		this(socketAddress, timeoutSeconds, Protocol.TCP, true);
	}

	/**
	 * Constructor.
	 * @param socketAddress the SocketAddress instance that corresponds to the server
	 * @param timeoutSeconds the number of seconds to wait for server reply before throwing a TimeoutException
	 * @param isAJavaApplication a boolean that should be set to true by default (false mode is used only for test purpose)
	 * @throws BasicClientException 
	 */
	protected BasicClient(InetSocketAddress socketAddress, int timeoutSeconds, Protocol protocol, boolean isAJavaApplication) throws BasicClientException {
		this.protocol = protocol;
		if (timeoutSeconds < 0) {
			throw new InvalidParameterException("The timeout delay must be equal to or greater than 0!");
		}
		this.timeout = timeoutSeconds;
		this.bypassTimeout = false;
		try {
			if (this.protocol == Protocol.TCP) {
				Socket socket = new Socket();
				socket.connect(socketAddress, 5000);
				socketWrapper = new TCPSocketWrapper(socket, isAJavaApplication);
			} else {
				socketWrapper = new UDPSocketWrapper(socketAddress, isAJavaApplication);
			}
			
		} catch (Exception e) {
			close();
			throw new BasicClientException(ExceptionType.ConnectionFailed);
		} 
		

		if (isAJavaApplication) {
			ServerReply replyFromServer = (ServerReply) readObjectFromServer();
			if (replyFromServer == ServerReply.CallAccepted) {
				open = true;
			} else if (replyFromServer == ServerReply.IAmBusyCallBackLater) {
				open = false;
				close();
				throw new BasicClientException(ExceptionType.ConnectionRejected);
			}
		} else {
			String replyFromServer = (String) readObjectFromServer();
			if (replyFromServer.equals(ServerReply.CallAccepted.name())) {
				open = true;
			} else if (replyFromServer.equals(ServerReply.IAmBusyCallBackLater.name())) {
				open = false;
				close();
				throw new BasicClientException(ExceptionType.ConnectionRejected);
			}
		}
		
	}

	
	protected Object readObjectFromServer() throws BasicClientException {
		try {
			if (bypassTimeout) {
				return socketWrapper.readObject();
			} else {
				return socketWrapper.readObject(timeout);
			}
		} catch (Exception e) {
			e.printStackTrace();
			close();
			throw handleException(e);
		} 
	}
	
	
	
	private BasicClientException handleException(Exception e) {
		if (e instanceof SocketTimeoutException) {
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
			if (socketWrapper != null & open) {
				socketWrapper.writeObject(ClientRequest.closeConnection);
			}
		} catch (IOException e) {
		} finally {
			try {
				if (socketWrapper != null) {
					socketWrapper.close();
				}
			} catch (IOException e) {}
		}
	}

	protected void setBypassTimeout(boolean bypass) {
		this.bypassTimeout = bypass;
	}
	
}