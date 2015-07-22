package org.umcn.me.util;

import java.io.File;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import org.umcn.me.samexternal.SAMSilentReader;
import org.umcn.me.util.ReadName;

public class ReadNameRetriever implements Iterable<ReadName> {

	public final File inFile;
	public final ReadNameOption option;
	
	
	public ReadNameRetriever(File inFile, ReadNameOption option){
		this.inFile = inFile;
		this.option = option;
	}
	

	/**
	 * This iterator will invoke close() when hasNext() returns false,
	 * if you need to dispose the iterator before make sure you call close your self
	 */
	@Override
	public CloseableIterator<ReadName> iterator() {
		final SAMSilentReader reader = new SAMSilentReader(this.inFile);
		final SAMRecordIterator recordIter = reader.iterator();
		final ReadNameOption readOption = this.option;
		
		CloseableIterator<ReadName> it = new CloseableIterator<ReadName>() {

			@Override
			public boolean hasNext() {
				if (! recordIter.hasNext()){
					this.close();
				}
				return recordIter.hasNext();
			}

			@Override
			public ReadName next() {
				SAMRecord rec = recordIter.next();
				return new ReadName(rec, readOption);
			}

			@Override
			public void remove() {
				throw new UnsupportedOperationException();
			}

			@Override
			public void close() {
				recordIter.close();
				reader.close();
			}
			
		};
		return it;
		
	}
	
}
