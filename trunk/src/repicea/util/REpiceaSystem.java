/*
 * This file is part of the repicea-util library.
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
package repicea.util;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import repicea.util.REpiceaTranslator.Language;

/**
 * The REpiceaSystem offers some additional features to the System class.
 * @author Mathieu Fortin - November 2014
 */
public class REpiceaSystem {

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

	/**
	 * This method returns the three digits of the JVM version
	 * @return a 1x3 array of integers
	 */
	public static int[] getJVMVersion() {
		String jvmVersion = ObjectUtility.getJVMVersion();
		return parseJVMVersion(jvmVersion);
	}

	private static int[] parseJVMVersion(String jvmVersion) {
		int[] version = new int[3];
		String[] splittedDigits = jvmVersion.split("\\.");
		version[0] = Integer.parseInt(splittedDigits[0]);
		version[1] = Integer.parseInt(splittedDigits[1]);
		version[2] = 0;
		if (splittedDigits.length > 2) {
			version[2] = Integer.parseInt(splittedDigits[2]);
		}
		return version;
	}
	
	/**
	 * This method returns true if the current version of the JVM is more recent than the parameter targetJVM
	 * @param targetJVM a String
	 * @param upToThirdDigit a boolean true to test up to the third digit
	 * @return a boolean
	 */
	public static boolean isCurrentJVMGreaterThanThisVersion(String targetJVM, boolean upToThirdDigit) {
		int[] currentVersion = getJVMVersion();
		int[] targetVersion = parseJVMVersion(targetJVM);
		if (currentVersion[0] > targetVersion[0]) {
			return true;
		} else if (currentVersion[1] > targetVersion[1]) {
			return true;
		} else if (upToThirdDigit && currentVersion[2] > targetVersion[2]) {
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
	public static boolean isCurrentJVMGreaterThanThisVersion(String targetJVM) {
		return isCurrentJVMGreaterThanThisVersion(targetJVM, false);
	}

//	public static void main(String[] args) {
//		String[] argTest = "repicea-console.jar -l cd".split(" ");
//			REpiceaSystem.setLanguageFromMain(argTest);
//	}
	
}
