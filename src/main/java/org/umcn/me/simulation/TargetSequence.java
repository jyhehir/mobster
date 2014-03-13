package org.umcn.me.simulation;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Arrays;
import java.util.Collections;

import org.apache.log4j.Logger;
import org.umcn.gen.sequence.InvalidSequenceException;
import org.umcn.gen.sequence.Sequence;

/**
 * This Class is a blueprint for TargetSequence objects. A TargetSequence object
 * is an object in which Sequences, like MobileSequences can be inserted. A
 * TargetSequence has optional exclusion_regions in which insertions may not
 * occur. Additionally a TargetSequence has an exclusion_window of 100 as
 * default. This means no insertion may occur in exclusion_regions +- 100 nt.
 * Also as of yet the insertSingle and insertMultiple method may only be called
 * once successfully for a given TargetSequence object. In this way
 * exclusion_regions do not have to get updated after each insertion, because
 * the sequence length and positions change after insertion. The updating of
 * exclusion_regions may be a future enhancement of this class.
 * 
 * @author Djie
 * 
 */
public class TargetSequence<S extends Sequence> extends Sequence {
	public static Logger logger = Logger.getLogger("TargetSequence");
	
	protected String name = "";
	protected int original_sequence_length;
	protected List<int[]> exclusion_regions = new ArrayList<int[]>();
	protected List<Integer> insertion_positions = new ArrayList<Integer>();
	protected int exclusion_window = 100;

	// after single or multiple insertion this should be set to false
	// because exclusion regions do not get updated after insertion events
	// leading to wrong values for some exclusion regions.
	protected boolean insertable = true;

	/**
	 * Make new TargetSequence. Important notes: default exclusion window for
	 * target region is 100. Use setExclusionWindow() after initialization to
	 * change this exclusion window.
	 * 
	 * @param seq
	 *            String containing nucleotide sequence.
	 * @throws InvalidSequenceException
	 */
	public TargetSequence(String seq) throws InvalidSequenceException {
		super(seq);
		this.original_sequence_length = this.sequence.length();
	}

	/**
	 * Add one extra exclusion region without overwriting previous exclusion regions.
	 * @param region integer array of size 2. First int specifies start (inclusive)
	 * of exclusion region. Second int specified end (exclusive) of exclusion region.
	 * @return wheter adding an exclusion region was succesfull
	 */
	public boolean addExclusionRegion(int[] region){
		int[] exclRegion = new int[2];
		
		if (region.length != 2) {
			logger.warn("exclusion regions not succesfully set: int[]" +
					" array exceeds size of 2 per region.");
			return false;
		}
		if (region[1] <= region[0]) {
			logger.warn("exclusion regions not succesfully set:" +
					" end of region should be higher int than beginning of region");
			return false;
		}
		if (region[0] < 0 || region[1] > this.original_sequence_length) {
			logger.warn("exclusion regions not succesfully set:" +
					" illegal values for coordinates: " + 
					"lower than 0 or larger than sequence length");
			return false;
		}
		exclRegion[0] = region[0];
		exclRegion[1] = region[1];
		this.exclusion_regions.add(exclRegion);
		return true;
	}
	
	public List<int[]> getExclusionRegions() {
		return this.exclusion_regions;
	}

	public int getExclusionWindow() {
		return this.exclusion_window;
	}

	public List<Integer> getInsertionPositions() {
		return this.insertion_positions;
	}

	public String getName() {
		return this.name;
	}

	public int getOriginalSequenceLength() {
		return this.original_sequence_length;
	}

	/**
	 * Get Target Site Duplication (TSD) Sequence at insertion position. TSD
	 * will always be left from insertion position. If TSD size is larger than
	 * sequence bounds than a TSD is returned until the boundary of the
	 * sequence. TSD can not be obtained from exclusion region + exclusion
	 * window.
	 * 
	 * @param insertionPos
	 *            position where insertion will take place (0-based)
	 * @param tsdSize
	 *            size of TSD to be returned.
	 * @return Sequence of TSD. Null when invalid insertionPos (i.e. larger than
	 *         sequence length or in exclusion region) or invalid size (i.e.
	 *         negative).
	 */
	public Sequence getTSD(int insertionPos, int tsdSize) {
		int tsdEnd = insertionPos;
		int tsdBeg = tsdEnd - tsdSize;
		Sequence tsd;

		if (tsdBeg < 0) {
			tsdBeg = 0;
		}
		try {
			tsd = new Sequence(this.sequence.substring(tsdBeg, tsdEnd));
			return tsd;
		} catch (Exception e) {
			// catches both StringIndexOutOfBounds as InvalidSequenceException
			System.out.println(e.getMessage() + " from method getTSD");
			return null;
		}
	}

	/**
	 * Insert a Sequence (s) into TargetSequence at position (0-based).
	 * insertSingle and insertMultiple can only be called once for a
	 * TargetSequence object.
	 * 
	 * @param position
	 *            Position into which Sequence (s) will be inserted. The first
	 *            nt of s will be at index specified by position in the original
	 *            Sequence.
	 * @param s
	 *            Sequence to insert into TargetSequence object.
	 * @return Boolean whether insertion was successful. Will for example return
	 *         false when insertion is in exclusion region + exclusion window.
	 */
	public boolean insertSingle(int position, Sequence s) {
		boolean succes = true;
		if (!isValidInsertion(position)) {
			succes = false;
		} else if (!this.insertable) {
			succes = false;
		} else {
			this.insert(position, s);
			this.insertable = false;
			this.insertion_positions.add(position);
		}
		return succes;
	}

	/**
	 * Insert multiple Sequences in TargetSequence object at specified
	 * positions. insertSingle and insertMultiple can only be called once for a
	 * TargetSequence object. This method inserts Sequences in reverse order
	 * (from highest insertion site to lowest). In this way it is guaranteed
	 * that all Sequences insert into the original positions of the
	 * TargetSequence.
	 * 
	 * @param insertMap
	 *            A map containing insertion positions (0-based) as keys and
	 *            Sequences (or subclasses) as values, to insert into the
	 *            TargetSequence at the specified position captured by the key.
	 * @param checkValidity
	 *            Boolean indicating whether the method should check if all
	 *            insertion sites are valid (i.e. not in exclusion areas +
	 *            exclusion windows). When false it is recommended you call
	 *            isValidInsertion() yourself.
	 * @return Boolean whether insertion was successful.
	 */
	public boolean insertMultiple(Map<Integer, S> insertMap,
			boolean checkValidity) {
		boolean succes = true;
		Integer[] insertPositions = insertMap.keySet().toArray(new Integer[0]);
		Arrays.sort(insertPositions, Collections.reverseOrder());
		StringBuilder sb = new StringBuilder(this.sequence);

		if (!this.insertable) {
			succes = false;
			logger.warn("No insertions performed: sequence is not insertable anymore");
		}

		if (checkValidity && succes) {
			for (int position : insertPositions) {
				if (!isValidInsertion(position)) {
					logger.warn("Not a valid insertion coordinate: insertion will not be performed");
					succes = false;
					break;
				}
			}
		}
		if (succes) {
			for (int position : insertPositions) {
				sb.insert(position, insertMap.get(position));
				//this.insert(position, insertMap.get(position));
				this.insertion_positions.add(position);
			}
			this.insertable = false;
		}
		
		this.sequence = sb.toString();
		logger.info("Inserted: " + insertPositions.length + " in " + this.getName());
		return succes;
	}

	public boolean isInsertable() {
		return this.insertable;
	}

	/**
	 * Function which checks if the supplied position to insert is valid: - Not
	 * smaller than 0 or larger than sequence length. - Not in any exclusion
	 * region + exclusion window.
	 * 
	 * @param pos
	 *            0-based insertion position
	 * @return whether insertion position is valid
	 * 
	 */
	public boolean isValidInsertion(int pos) {
		boolean valid = true;
		if (pos > this.original_sequence_length || pos < 0) {
			valid = false;
		} else {
			for (int[] exclRegion : this.exclusion_regions) {
				int exclBeg = exclRegion[0] - this.exclusion_window;
				int exclEnd = exclRegion[1] + this.exclusion_window;
				if (pos >= exclBeg && pos < exclEnd) {
					valid = false;
					break;
				}
			}
		}
		return valid;
	}

	/**
	 * Set Exclusion Regions: regions in which no insertions may take place.
	 * Previous set exclusion regions will be overridden.
	 * 
	 * @param exclRegions
	 *            List containing integer arrays of size 2, each marking the
	 *            beginning and end region (0-based, beginning is inclusive, end
	 *            is non-inclusive) of an exclusion region.
	 * @return Boolean whether setting of exclusion regions was successful.
	 *         Returns false when supplied List for exclRegions is not valid.
	 */
	public boolean setExclusionRegions(List<int[]> exclRegions) {
		for (int[] region : exclRegions) {
			if (region.length != 2) {
				logger.warn("exclusion regions not succesfully set: int[]" +
						" array exceeds size of 2 per region.");
				return false;
			}
			if (region[1] <= region[0]) {
				logger.warn("exclusion regions not succesfully set:" +
						" end of region should be higher int than beginning of region");
				return false;
			}
			if (region[0] < 0 || region[1] > this.original_sequence_length) {
				logger.warn("exclusion regions not succesfully set:" +
						" illegal values for coordinates: " + 
						"lower than 0 or larger than sequence length");
				return false;
			}
		}
		this.exclusion_regions = new ArrayList<int[]>(exclRegions);
		logger.info("Nr of exclusion regions set: " + exclRegions.size());
		return true;
	}

	/**
	 * Set exclusion window. Insertions are only valid when they do NOT occur in
	 * exclusion regions +/- the exclusion window. When a TargetSequence is
	 * initialized a default exclusion window of 100 is set.
	 * 
	 * @param exclWindow
	 *            integer specifying the exclusion window.
	 * @return boolean whether setting the exclusion window was succesful.
	 */
	public boolean setExclusionWindow(int exclWindow) {
		if (exclWindow < 0 || exclWindow >= this.original_sequence_length) {
			return false;
		} else {
			this.exclusion_window = exclWindow;
			return true;
		}
	}

	public void setName(String name) {
		this.name = name;
	}
	
	public void printExclusionRegions(){
		StringBuilder sb = new StringBuilder();
		for (int[] region : this.exclusion_regions){
			sb.append("{");
			sb.append(region[0]);
			sb.append(", ");
			sb.append(region[1]);
			sb.append("}");
			sb.append(" ");
			
		}
		System.out.println(sb);
	}
}
