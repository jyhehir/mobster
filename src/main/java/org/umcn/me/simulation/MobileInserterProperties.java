package org.umcn.me.simulation;

import java.lang.reflect.Field;

import org.apache.log4j.Logger;

/**
 * Class containing settings for MobileInserter class:
 * 
 * How many l1, sva, alu elements should be inserted and how many of them should
 * contain tsds, inversions or truncations
 * 
 * This class makes use of a fluent interface
 * (http://en.wikipedia.org/wiki/Fluent_interface)
 * 
 * So when initializing MobileInserterProperties you do not have to pass any
 * parameters. All private attributes are set to their default values and you
 * can subsequently set the value of any number of attributes using the
 * following syntax:
 * 
 * //initialize with default values MobileInserterProperties props = new
 * MobileInserterProperties();
 * 
 * //set 3 parameters by method chaining
 * props.setAluInsertionNrs(10).setSVAInsertionNrs(7).setL1InsertionNrs(20);
 * 
 * TODO
 * not for all private variables a setter or getter is provided, should still
 * be implemented
 * 
 * @author Djie
 * 
 */
public class MobileInserterProperties {
	public static Logger logger = Logger.getLogger("MobileInserterProperties");
	
	//IMPORTANT:
	//WHEN ADDING EXTRA MOBILE ELEMENT FAMILIES LIKE HERV-K ELEMENTS
	//DO NOT FORGET TO UPDATE THE GETTOTAL METHODS (LIKE GETTOTALINSERTIONS)

	private int l1_insertions = 700;
	private int sva_insertions = 700;
	private int alu_insertions = 700;
	private int herv_insertions = 0;

	private int l1_tsds = 0;
	private int sva_tsds = 0;
	private int alu_tsds = 0;
	private int herv_tsds = 0;
	private int min_tsd = 3;
	private int max_tsd = 5;

	private int l1_inversions = 0;
	private int sva_inversions = 0;
	private int alu_inversions = 0;
	private double max_inversion_perc = 0.5;

	private int l1_left_truncations = 0;
	private int sva_left_truncations = 0;
	private int alu_left_truncations = 0;
	private double max_left_truncation = 0.2;

	private int l1_right_truncations = 0;
	private int sva_right_truncations = 0;
	private int alu_right_truncations = 0;
	private double max_right_truncation = 0.2;
	
	private int exclusion_window = 100;
	private int inclusion_window = 35;

	public MobileInserterProperties() {

	}

	public int getAluInsertionNrs() {
		return this.alu_insertions;
	}

	public int getL1InsertionNrs() {
		return this.l1_insertions;
	}
	
	public int getHervInsertionNrs() {
		return this.herv_insertions;
	}

	public int getSVAInsertionNrs() {
		return this.sva_insertions;
	}

	public int getAluTsdNrs() {
		return this.alu_tsds;
	}
	
	public int getHervTsdNrs() {
		return this.herv_tsds;
	}

	public int getL1TsdNrs() {
		return this.l1_tsds;
	}

	public int getSVATsdNrs() {
		return this.sva_tsds;
	}

	public int getAluInversionNrs() {
		return this.alu_inversions;
	}

	public int getL1InversionNrs() {
		return this.l1_inversions;
	}

	public int getSVAInversionNrs() {
		return this.sva_inversions;
	}

	public int getExclusionWindow(){
		return this.exclusion_window;
	}
	
	public double getMaxInversion() {
		return this.max_inversion_perc;
	}

	public int getMinTsdSize() {
		return this.min_tsd;
	}

	public int getMaxTsdSize() {
		return this.max_tsd;
	}

	public int getTotalInsertions() {
		return this.alu_insertions + this.l1_insertions + this.sva_insertions + this.herv_insertions;
	}
	
	public int getTotalTsds(){
		return this.alu_tsds + this.l1_tsds + this.sva_tsds + this.herv_tsds;
	}
	
	public int getInclusionWindow(){
		return this.inclusion_window;
	}
	

	public MobileInserterProperties setAluInsertionNrs(int nr) {
		if (nr >= 0) {
			this.alu_insertions = nr;
			return this;
		}
		logger.warn("Illegal value for Alu insertion nr: no setting will be performed");
		return this;
	}

	public MobileInserterProperties setL1InsertionNrs(int nr) {
		if (nr >= 0) {
			this.l1_insertions = nr;
			return this;
		}
		logger.warn("Illegal value for L1 insertion nr: no setting will be performed");
		return this;
	}

	public MobileInserterProperties setSVAInsertionNrs(int nr) {
		if (nr >= 0) {
			this.sva_insertions = nr;
			return this;
		}
		logger.warn("Illegal value for SVA insertion nr: no setting will be performed");
		return this;
	}

	public MobileInserterProperties setAluTsdNrs(int nr) {
		if (nr >= 0 && nr <= this.alu_insertions) {
			this.alu_tsds = nr;
			return this;
		}
		logger.warn("Illegal value for Alu TSD nr: no setting will be performed");
		return this;
	}

	public MobileInserterProperties setL1TsdNrs(int nr) {
		if (nr >= 0 && nr <= this.l1_insertions) {
			this.l1_tsds = nr;
			return this;
		}
		logger.warn("Illegal value for L1 TSD nr: no setting will be performed");
		return this;
	}

	public MobileInserterProperties setSVATsdNrs(int nr) {
		if (nr >= 0 && nr <= this.sva_insertions) {
			this.sva_tsds = nr;
			return this;
		}
		logger.warn("Illegal value for SVA TSD nr: no setting will be performed");
		return this;
	}

	public MobileInserterProperties setAluInversionNrs(int nr) {
		if (nr >= 0 && nr <= this.alu_insertions) {
			this.alu_inversions = nr;
			return this;
		}
		logger.warn("Illegal value for Alu Inversion nr: no setting will be performed");
		return this;
	}

	public MobileInserterProperties setL1InversionNrs(int nr) {
		if (nr >= 0 && nr <= this.l1_insertions) {
			this.l1_inversions = nr;
			return this;
		}
		logger.warn("Illegal value for L1 Inversion nr: no setting will be performed");
		return this;
	}

	public MobileInserterProperties setSVAInversionNrs(int nr) {
		if (nr >= 0 && nr <= this.sva_insertions) {
			this.sva_inversions = nr;
			return this;
		}
		logger.warn("Illegal value for SVA Inversion nr: no setting will be performed");
		return this;
	}
	
	public MobileInserterProperties setInclusionWindow(int nr){
		this.inclusion_window = nr;
		return this;
	}
	
	public MobileInserterProperties setExclusionWindow(int nr){
		this.exclusion_window = nr;
		return this;
	}

	@Override
	public String toString() {
		String s = "Settings set in MobileInserterProperties:\n";
		String fieldName;
		Field[] fields = MobileInserterProperties.class.getDeclaredFields();
		for (Field i : fields) {
			fieldName = i.getName();
			try {
				if (fieldName != "logger") {
					s += i.getName();
					s += ": ";
					s += i.get(this);
					s += "\n";
				}
			} catch (IllegalArgumentException e) {
				logger.warn("error in toString " + e.getMessage());
			} catch (IllegalAccessException e) {
				logger.warn("error in to String " + e.getMessage());
			}

		}
		return s;
	}
}
