package org.umcn.me.pairedend;

import java.io.File;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.sf.samtools.*;

import org.umcn.gen.sam.SAMSilentReader;
import org.umcn.me.sam.InvalidCategoryException;
import org.umcn.me.util.CollectionUtil;
import org.umcn.me.util.MobileDefinitions;

public class Testing {
	public static void main(String[] args) throws InvalidCategoryException {
		String[] categorySplit;
		String mobileCategory;
		String category = "Sadhu10-1.ME_CAT.SINE";
		
		categorySplit = category.split(MobileDefinitions.OTHER_MOBILE_CATEGORY_ATTRIBUTE, -1);
		mobileCategory = categorySplit[categorySplit.length - 1].trim();
		if ( mobileCategory.length() == 0 ) {
			throw new InvalidCategoryException (category + "can not be categorized\n." +  
		 " Tried to split on " + MobileDefinitions.OTHER_MOBILE_CATEGORY_ATTRIBUTE);
		}
		System.out.println(mobileCategory);
		
	}

}
