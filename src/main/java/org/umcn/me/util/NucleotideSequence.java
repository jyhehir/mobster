package org.umcn.me.util;

public class NucleotideSequence{
	
	public static final char AMBIGUOUS_BASE_CHARACTER = 'N';
	
	protected String sequence = "";
	
	public NucleotideSequence(String seq) throws InvalidNucleotideSequenceException{
		if (! isValidSequence(seq)) {
			throw new InvalidNucleotideSequenceException("Sequence not valid: " + seq);
		}
		
		this.sequence = seq;
	}

	public String getSequence(){
		return this.sequence;
	}

	public int compareTo(NucleotideSequence o) {
		return this.sequence.compareTo(o.sequence);
	}
	
	public int getLength() {
		return sequence.length();
	}
	
	public int getNumberOf(char bp) {
		if (this.sequence == null) {
			return -1;
		}
		
		int c = 0;
		for(char s : sequence.toCharArray()) {
			if (s == bp) {
				c++;
			}
		}
		return c;
	}
	
	
	public boolean hasSequenceAtPosition(String seq, int position) {
		return sequence.substring(position, seq.length()).equals(seq);
	}
	
	public void complement() {
		this.sequence = getComplementSequence().getSequence();
	}
	
	// reverts this sequence
	public void reverse() {
		StringBuilder tmp_seq = new StringBuilder(this.sequence);
		this.sequence = tmp_seq.reverse().toString();
	}
	
	public void reverseComplement() {
		this.reverse();
		this.complement();
	}
	
	/**
	 * 
	 * @return Complementary sequence:
	 * a or A --> T
	 * c or C --> G
	 * g or G --> C
	 * t or T --> A
	 * N --> N
	 * 
	 * Any other unknown character will remain the same!
	 * 
	 */
	public NucleotideSequence getComplementSequence() {
		char[] complement = new char[sequence.length()];
		
		int i =0;
		char new_base = 'z';
		for (char ch : sequence.toCharArray()) {			
			if (ch == AMBIGUOUS_BASE_CHARACTER ) {
				new_base = AMBIGUOUS_BASE_CHARACTER;
			}else if (ch == 'a' || ch == 'A') {
				new_base = 'T';
			}else if (ch == 'c' || ch == 'C') {
				new_base = 'G';
			}else if (ch == 'g' || ch == 'G') {
				new_base = 'C';
			}else if (ch == 't' || ch == 'T') {
				new_base = 'A';
			}else {
				new_base = ch;
				System.err.println("[S] Could not get complement base of " + String.valueOf(ch) + " in sequence" + sequence);
			} 
			complement[i++] = new_base;
		}
		try {
			return new NucleotideSequence(String.copyValueOf(complement));
		} catch (InvalidNucleotideSequenceException e) {
			//should not occur
			return null;
		}
	}
	
	
	public boolean equals(Object o) {
		NucleotideSequence s = (NucleotideSequence) o;
		return s.getSequence().equalsIgnoreCase(this.getSequence());
	}
	
	
	
	public String toString () {
		return this.sequence;
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((sequence == null) ? 0 : sequence.hashCode());
		return result;
	}
	


	public static boolean isValidSequence(String seq) {
		seq = seq.toUpperCase();
		for (char c : seq.toCharArray()) {
			if (c != 'C' && c != 'T' && c != 'G' && c != 'A' && c != AMBIGUOUS_BASE_CHARACTER) {
				return false;
			}
		}
		return true;
	}
	

	/**
	 * 
	 * @param min_length minimum length of poly A or poly T stretch at start or end of sequence
	 * @param maxperc_non_a_t maximum bases in poly A or poly T stretch which is respectively
	 * not a or not T
	 * @return polyA or polyT when sequence applies to poly A or poly T thresholds.
	 */
	public String getPolyAOrTMapping(int min_length, int max_non_a_t, boolean start){
		String polymer_stretch = "";
		NucleotideSequence seq;
		if (this.sequence.length() >= min_length){
			try {
				if (start){
					seq = new NucleotideSequence(sequence.substring(0, min_length));
				}else{
					seq = new NucleotideSequence(sequence.substring(sequence.length() - min_length, sequence.length()));
				}

				if (seq.getNumberOf('A') >= (min_length - max_non_a_t)){
					polymer_stretch = "polyA";
				}else if (seq.getNumberOf('T') >= (min_length - max_non_a_t)){
					polymer_stretch = "polyT";
				}
			} catch (InvalidNucleotideSequenceException e) {
				//Should not occur
				return polymer_stretch;
			}
		}		
		return polymer_stretch;	
	}

	
	
	// this is a 1 based 
	public NucleotideSequence getSubSequence(int start, int end) {
		if (start <1 || end >= this.sequence.length()) {
			return null;
		}
		try {
			return new NucleotideSequence(this.sequence.substring(start, end));
		} catch (InvalidNucleotideSequenceException e) {
			return null;
		}
	}
	


}