package org.umcn.me.simulation;

import java.util.ArrayList;
import java.util.List;

import org.umcn.gen.sequence.InvalidSequenceException;
import org.umcn.gen.sequence.Sequence;

/**
 * MobileSequence Class is a blueprint for Mobile Element Sequences. As Mobile
 * Elements often get 5' truncated or have inversions, the MobileSequence Class
 * provides these types of functionalities.
 * 
 * @author Djie
 * 
 */
public class MobileSequence extends Sequence {
	//NOTE: WHEN ADDING/REMOVING ATTRIBUTES
	//DO NOT FORGET TO ALSO ADD/REMOVE THEM FROM COPY CONSTRUCTOR
	protected String family_name = "";
	protected String name = "";
	protected String description = "";

	protected int left_truncation_size = 0;
	protected int right_truncation_size = 0;
	protected boolean inserts_with_tsd = false;
	protected int tsd_size = 0;

	protected List<int[]> inversion_coordinates = new ArrayList<int[]>();

	public MobileSequence(String seq) throws InvalidSequenceException {
		super(seq);
	}
	
	/**
	 * Copy constructor. For deep copying a MobileSequence object
	 * @param s MobileSequence to be deep copied
	 * @throws InvalidSequenceException
	 */
	public MobileSequence(MobileSequence s) throws InvalidSequenceException {
		//Make sure all attributes get copied when changing attributes in this class!
		super(s.getSequence());
		this.family_name = s.getFamilyName();
		this.name = s.getName();
		this.description = s.getDescription();
		this.left_truncation_size = s.getLeftTruncationSize();
		this.right_truncation_size = s.getRightTruncationSize();
		this.inserts_with_tsd = s.insertsWithTsd();
		this.inversion_coordinates = s.getInversionCoordinates();
		this.tsd_size = s.getTsdSize();
	}

	public String getDescription() {
		return this.description;
	}

	public String getFamilyName() {
		return this.family_name;
	}

	public List<int[]> getInversionCoordinates() {
		return this.inversion_coordinates;
	}

	public int getLeftTruncationSize() {
		return this.left_truncation_size;
	}

	public String getName() {
		return this.name;
	}

	public int getRightTruncationSize() {
		return this.right_truncation_size;
	}
	
	public int getTsdSize(){
		return this.tsd_size;
	}
	
	public boolean insertsWithTsd(){
		return this.inserts_with_tsd;
	}

	/**
	 * Inverse (reverse) part of the sequence.<br>
	 * The part of the sequence to be reversed is indicated<br>
	 * by the from parameter (which includes the integer provided by from)<br>
	 * and end parameter (which does not include the integer provided by end).<br>
	 * Positions in from and end are 0-based. I.e. first nucleotide is numbered
	 * 0.<br>
	 * Inversion will not be performed when inversion coordinates overlap with
	 * previous inversion events.
	 * 
	 * @param from
	 *            sequence position from (including integer)
	 * @param to
	 *            sequence position to (up to, but not including integer)
	 * @return boolean indicating whether inversion was successful.
	 */
	public boolean inverseSubSequence(int from, int to) {
		int[] newCoordinates = { from, to };
		List<int[]> coordinates = this.getInversionCoordinates();
		if (!isValidInversion(from, to)) {
			return false;
		}
		String subSeq = this.sequence.substring(from, to);
		StringBuffer revSeq = new StringBuffer(subSeq);

		try {
			this.substitute(from, new Sequence(revSeq.reverse().toString()));
			coordinates.add(newCoordinates);
			this.setInversionCoordinates(coordinates);
			return true;
		} catch (InvalidSequenceException e) {
			e.printStackTrace();
		}
		return false;

	}

	/**
	 * Checks whether insertion coordinates are valid: Whether they are not
	 * negative, larger than sequence size, whether first coordinate is smaller
	 * than second coordinate and whether coordinates do not overlap with
	 * previous inversion events
	 * 
	 * @param from
	 *            0-based integer: first base of subsequence to inverse
	 * @param to
	 *            0-based integer: up to but not including last base of
	 *            subsequence to inverse
	 * @return whether Inversion event is valid
	 */
	private boolean isValidInversion(int from, int to) {
		boolean valid = true;
		if (from >= to - 1) {
			System.out.println("Failed from >= to - 1");
			valid = false;
		} else if (from < 0 || from > this.sequence.length()) {
			valid = false;
			System.out.println("from is invalid (<0 or > seq.length)");
		} else if (to < 0 || to > this.sequence.length()) {
			valid = false;
			System.out.println("to is invalid (<0 or > seq.length)");
		} else {
			// Check for each existing inversion coordinate whether
			// there is an overlap with the new coordinates (from and to)
			for (int[] coordinate : this.inversion_coordinates) {
				int begin = coordinate[0];
				int end = coordinate[1];

				if (to <= begin) {
					continue;
				} else if (from >= end) {
					continue;
				} else {
					valid = false;
					System.out
							.println("Failed because of overlap in inversion event");
					break;
				}
			}
		}
		return valid;
	}

	public void setDescription(String desc) {
		this.description = desc;
	}

	public void setFamilyName(String fname) {
		this.family_name = fname;
	}
	
	public void setInsertWithTsd(boolean tsd){
		this.inserts_with_tsd = tsd;
	}

	/**
	 * @param coordinates
	 *            note this should be a List containing both old and new
	 *            inversion coordinates! Otherwise old inversion coordinates get
	 *            deleted. Each inversion event has to be represented by one
	 *            int[] Array of size 2. First int should be lower than second
	 *            int.
	 */
	private void setInversionCoordinates(List<int[]> coordinates) {
		this.inversion_coordinates = coordinates;
	}

	public void setName(String name) {
		this.name = name;
	}
	
	public boolean setTsdSize(int size){
		if(size > 0){
			this.tsd_size = size;
			return true;
		}
		return false;
	}

	/**
	 * TODO: do not truncate TSD from right site of ME
	 * Truncates the sequence from the left or right site. <b>Example:</b><br>
	 * MobileSequence.getSequence() = "ACGT"<br>
	 * MobileSequence.truncate(true, 2).getSequence() = "GT"<br>
	 * 
	 * When original sequence = "ACGT" again:<br>
	 * MobileSequence.truncate(false, 2).getSequence() = "AC"
	 * 
	 * @param fromLeft
	 *            Truncate sequence from left side when true
	 * @param size
	 *            Size of sequence that will be truncated.
	 * @return boolean whether truncation was successful. False will be returned
	 *         for instance when truncation size is bigger than sequence size.
	 */
	public boolean truncate(boolean fromLeft, int size) {
		boolean willTruncate = false;
		int originalLength = this.sequence.length() - this.tsd_size;
		if (size < originalLength && size > 0) {
			if (fromLeft) {
				this.sequence = this.sequence.substring(size);
				this.left_truncation_size += size;
				willTruncate = true;

			} else {
				this.sequence = this.sequence.substring(0,
						this.sequence.length() - size);
				this.right_truncation_size += size;
				willTruncate = true;
			}
		}
		return willTruncate;
	}
}
