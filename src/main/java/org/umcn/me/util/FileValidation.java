package org.umcn.me.util;

import java.io.File;

public class FileValidation {
    public static boolean fileExists(String filePath){
        return new File(filePath).isFile();
    }

    public static boolean fileReadable(String filePath){
        return new File(filePath).canRead();
    }

    public static boolean fileValid(String filePath){
        return fileExists(filePath) && fileReadable(filePath);
    }
}
