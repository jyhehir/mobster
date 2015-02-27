package org.umcn.me.pairedend;

import java.io.Serializable;


import net.sf.samtools.SAMRecord;

public class SAMRecordSerializable implements Serializable {

	public final SAMRecord record;
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -8142733422799759236L;
	
	public SAMRecordSerializable(SAMRecord record){
		this.record = record;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((record == null) ? 0 : record.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		SAMRecordSerializable other = (SAMRecordSerializable) obj;
		if (record == null) {
			if (other.record != null)
				return false;
		} else if (!record.equals(other.record))
			return false;
		return true;
	}


}
