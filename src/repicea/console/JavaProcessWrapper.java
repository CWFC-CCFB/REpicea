/*
 * This file is part of the repicea-console library.
 *
 * Copyright (C) 2012 Mathieu Fortin for Rouge-Epicea
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed with the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * Please see the license at http://www.gnu.org/copyleft/lesser.html.
 */
package repicea.console;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.List;

import repicea.app.AbstractGenericTask;
import repicea.multiprocess.JavaProcess;
import repicea.multiprocess.JavaProcess.JVM_OPTION;

public class JavaProcessWrapper extends AbstractGenericTask implements PropertyChangeListener {

	private static final long serialVersionUID = 20120218L;

	private JavaProcess internalProcess;
	private boolean atLeastOneMessageReceived;

	
	public JavaProcessWrapper(String taskName, List<String> commands, File workingDirectory) {
		setName(taskName);
		internalProcess = new JavaProcess(commands, workingDirectory);
		internalProcess.redirectOutputStream(false);
		internalProcess.addPropertyChangeListener(this);
	}
	
	
	public JavaProcess getInternalProcess() {return internalProcess;}
	
	@Override
	public void cancel() {
		internalProcess.cancel(true);
		super.cancel();
	}
	
	@Override
	public void propertyChange(PropertyChangeEvent arg0) {
		if (arg0.getSource().equals(internalProcess)) {
			String propertyName = arg0.getPropertyName();
			if (propertyName.equals("MessageReceived")) {
				if (!atLeastOneMessageReceived) {
					atLeastOneMessageReceived = true;
				}
				System.out.println((String) arg0.getNewValue());
			}
		}
	}


	@Override
	protected void doThisJob() throws Exception {
		System.out.println("Launching " + getName() + "...");
		for (JVM_OPTION option : internalProcess.getJVMSettings().keySet()) {
			System.out.println(option.name() + " = " + internalProcess.getJVMSettings().get(option).toString());
		}
		internalProcess.execute();
		int output = -1;
		try {
			output = internalProcess.get();
		} catch (Exception e) {
			if (hasBeenCancelled()) {
				output = 0;
			} else {
				throw e;
			}
		} 
		if (output == 0 && atLeastOneMessageReceived) {
			return;
		} else {
			throw new Exception("Process exited with value " + output);
		}
	}

	
	
}
