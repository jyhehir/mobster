package org.umcn.me.util;

import java.util.Iterator;

/**
 * Iterator which adds a method to close a resource
 * 
 * @author Djie Tjwan Thung
 *
 * @param <E>
 */
public interface CloseableIterator<E> extends Iterator<E> {

	void close();
	
}
