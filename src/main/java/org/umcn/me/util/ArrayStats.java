package org.umcn.me.util;

import org.omg.CORBA.DynAnyPackage.InvalidValue;

public class ArrayStats {

    //Function to add a number to an array
    public static void add(int[] a, int start, int end, int num){
        for(int i = start; i <= end; i++){
            a[i] += num;
        }
    }
    public static void add(double[] a, int start, int end, double num){
        for(int i = start; i <= end; i++){
            a[i] += num;
        }
    }

    //Function to get the maximum in an array
    public static int max(int[] a) throws InvalidValue {

        if(a.length==0){
            throw new InvalidValue("Array is empty");
        }
        int max = a[0];

        for(int i=1; i < a.length; i++){
            if((a[i] > max)){
                max = a[i];
            }
        }

        return max;
    }
    public static double max(double[] a) throws InvalidValue {

        if(a.length==0){
            throw new InvalidValue("Array is empty");
        }
        double max = a[0];

        for(int i=1; i < a.length; i++){
            if(!Double.isNaN(a[i]) && (a[i] > max || Double.isNaN(max))){
                max = a[i];
            }
        }

        return max;
    }

    //Function to take the mean of the values in an array
    public static int mean(int[] a) throws InvalidValue {
        if(a.length==0){
            throw new InvalidValue("Array is empty");
        }

        int total = 0;
        for(int i=0; i < a.length; i++){
            total += a[i];
        }
        return (int) ((double)total / (double)a.length);
    }
    public static int mean(double[] a) throws InvalidValue {
        if(a.length==0){
            throw new InvalidValue("Array is empty");
        }

        int total = 0;
        for(int i=0; i < a.length; i++){
            total += a[i];
        }
        return total / a.length;
    }

    public static boolean isNumeric(String str) {
        try {
            Double.parseDouble(str);
            return true;
        } catch(NumberFormatException e){
            return false;
        }
    }
}
