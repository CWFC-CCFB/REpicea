/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2014 Mathieu Fortin for Rouge Epicea.
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
package repicea.lang;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import repicea.util.REpiceaTranslator;
import repicea.util.REpiceaTranslator.Language;

/**
 * The REpiceaSystem offers some additional features to the System class.
 * 
 * @author Mathieu Fortin - November 2014
 */
public class REpiceaSystem {

	private static final double ConversionFromBytestoMegabytes = 1d / (1024 * 1024);
	
	private static String jreVersion;
	
	private static String revision;

	
	static {
		String completeJREVersion = System.getProperty("java.version");
		try {
			jreVersion = completeJREVersion.substring(0, completeJREVersion.indexOf("_"));
			revision = completeJREVersion.substring(completeJREVersion.indexOf("_") + 1);
		} catch (Exception e) {
			jreVersion = completeJREVersion;
			revision = "unknown";
		}
	}

	
	/**
	 * This method returns the temporary input/output directory. It is preferable to the System.getProperty("java.io.tmpdir") method
	 * because it makes sure the path ends with a file separator.
	 * @return a String
	 */
	public static String getJavaIOTmpDir() {
		String directory = System.getProperty("java.io.tmpdir");
		String osName = System.getProperty("os.name");
		if (osName.equals("Linux")) {
			directory = directory.concat(File.separator);
		}
		return directory;
	}
	
	/**
	 * This method scans the arguments given to the main method in order to 
	 * set the language. Does not do anything if the arguments do not contain
	 * "-l" followed by the language code, typically "fr" for French
	 * @param args the arguments given to the main method
	 * @param defaultLanguage a default language can be null
	 */
	public static void setLanguageFromMain(String[] args, Language defaultLanguage) {
		Language language = null;
		List<String> argumentList = Arrays.asList(args);
		if (argumentList.contains("-l")) {
			int indexLanguage = argumentList.indexOf("-l") + 1;
			if (argumentList.size() > indexLanguage) {
				String languageCode = argumentList.get(indexLanguage);
				language = REpiceaTranslator.Language.getLanguage(languageCode);
			}
		}
		if (language == null && defaultLanguage != null) {
			language = defaultLanguage;
		}
		if (language != null) {
			REpiceaTranslator.setCurrentLanguage(language);
		}
	}
	

	/**
	 * This method scans the arguments given to the main method in order to 
	 * set the language. Does not do anything if the arguments do not contain
	 * "-l" followed by the language code, typically "fr" for French
	 * @param args the arguments given to the main method
	 */
	public static void setLanguageFromMain(String[] args) {
		setLanguageFromMain(args, null);
	}

//	/**
//	 * This method returns the three digits of the JVM version
//	 * @return a 1x3 array of integers
//	 */
//	public static int[] getJVMVersion() {
//		String jvmVersion = ObjectUtility.getJVMVersion();
//		return parseJVMVersion(jvmVersion);
//	}

	private static int parseJVMVersion(String jvmVersion) {
		String[] splittedDigits = jvmVersion.split("\\.");
		int firstInteger = Integer.parseInt(splittedDigits[0]);
		if (firstInteger != 1) {
			return firstInteger;	// prior to Java 11 the string is something like 1.8.0. After that it is something like 11.0.6
		} else {
			return Integer.parseInt(splittedDigits[1]);
		}
	}
	
	/**
	 * This method returns true if the current version of the JVM is more recent than the parameter targetJVM
	 * @param targetJVM a String
	 * @param upToThirdDigit a boolean true to test up to the third digit
	 * @return a boolean
	 */
	public static boolean isCurrentJVMLaterThanThisVersion(String targetJVM, boolean upToThirdDigit) {
		int currentVersion = parseJVMVersion(getJVMVersion());
		int targetVersion = parseJVMVersion(targetJVM);
		if (currentVersion > targetVersion) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * This method returns true if the current version of the JVM is more recent than the parameter targetJVM. The 
	 * test is carried out on the first and second digit only. For instance, versions 1.7.12 and 1.7.17 would be equally
	 * recent.
	 * @param targetJVM a String
	 * @return a boolean
	 */
	public static boolean isCurrentJVMLaterThanThisVersion(String targetJVM) {
		return isCurrentJVMLaterThanThisVersion(targetJVM, false);
	}


	/**
	 * This static methods processes the arguments given to the main function in a very classical manner. 
	 * The language is set if the args parameter contains -l followed by either en or fr. Then the arguments
	 * are returned as a List of String instances for further specific processing.
	 * @param args the input array of String instances
	 * @return a List of String instances
	 */
	public static List<String> setClassicalOptions(String[] args) {
		if (args.length == 0) {
			System.out.println("No parameters received.");
			return new ArrayList<String>();
		} else {
			String inputString = "";
			for (String str : args) {
				inputString = inputString + str + "; ";
			}
			System.out.println("Parameters received:" + inputString);
			REpiceaSystem.setLanguageFromMain(args, Language.English);
			System.out.println("Language set to: " + REpiceaTranslator.getCurrentLanguage().name());
			return Arrays.asList(args);
		}
	}
	
	/**
	 * This method scan the arguments for a particular option and returns the value of that option.
	 * @param option the option
	 * @param argumentList a List of String instances
	 * @return a String or null if the option was not found
	 */
	public static String retrieveArgument(String option, List<String> argumentList) {
		if (argumentList.contains(option) && (argumentList.indexOf(option) + 1 < argumentList.size())) {
			return argumentList.get(argumentList.indexOf(option) + 1);
		} else {
			return null;
		}
	}
	
	
	/**
	 * This method returns the version of the virtual machine.
	 * @return a String
	 */
	public static String getJVMVersion() {return jreVersion;}
	
	/**
	 * This method returns the revision of the virtual machine.
	 * @return a String
	 */
	public static String getJVMRevision() {return revision;}

	/**
	 * Returns true if the OS is Windows or false otherwise.
	 * @return a boolean
	 */
	public static boolean isRunningOnWindows() {
		if (System.getProperty("os.name").startsWith("Windows")) {
			return true;
		} else {
			return false;
		}
	}
	
	/**
	 * Returns the architecture of the JVM, i.e. either 32-Bit, 64-Bit or unknown.
	 * @return a String
	 */
	public static String getJavaArchitecture() {
		return System.getProperty("sun.arch.data.model");
	}
	


	/**
	 * Call the garbage collector and compute the current memory load in Megabytes. <br>
	 * <br>
	 * This method actually relies on the difference between to total memory and
	 * the free memory as provided by the Runtime class.
	 * @return the memory load in Megabytes
	 */
	public static double getCurrentMemoryLoadMb() {
		System.gc();
		long usedMem = (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
		return usedMem * ConversionFromBytestoMegabytes;
	}
}
