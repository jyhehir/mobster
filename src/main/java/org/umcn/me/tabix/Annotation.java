package org.umcn.me.tabix;

import java.text.ParseException;

public abstract class Annotation {

	abstract Annotation parseFromLine(String line) throws ParseException;

}
