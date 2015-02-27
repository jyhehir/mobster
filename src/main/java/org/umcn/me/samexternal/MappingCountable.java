package org.umcn.me.samexternal;

/**
 * @author Djie
 */

public interface MappingCountable {
	
	public boolean isMappedUniquely();
	
	public boolean isMappedMultiple();
	
	public int getNumberOfMappings();
}
