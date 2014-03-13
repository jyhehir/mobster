package org.umcn.me.simulation;

import org.apache.log4j.BasicConfigurator;


/**
 * Short test class for MobileInserterProperties
 * 
 *TODO make a JUNIT test class out of this
 */
public class PropTest {
	public static void main(String[] args) {
		BasicConfigurator.configure();
		MobileInserterProperties props = new MobileInserterProperties();
		System.out.println(props.getL1InsertionNrs());
		System.out.println(props.getAluInsertionNrs());
		System.out.println(props.getSVAInsertionNrs());
		props.setAluInsertionNrs(10).setSVAInsertionNrs(7).setL1InsertionNrs(20);
		System.out.println(props.getL1InsertionNrs());
		System.out.println(props.getAluInsertionNrs());
		System.out.println(props.getSVAInsertionNrs());
		props.setAluInsertionNrs(-5).setSVAInsertionNrs(17);
		System.out.println(props.getL1InsertionNrs());
		System.out.println(props.getAluInsertionNrs());
		System.out.println(props.getSVAInsertionNrs());
		System.out.println(props.toString());

	}
}
