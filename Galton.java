// c:/ddc/Java/Galton.java
// (C) Dennis de Champeaux/ OntoOO
// Tue Jan 18 19:57:12 2022

import java.util.*;
// Show regression TO and FROM the mean assuming m=f (with assortative mating)
public class Galton { 
    static double c  = 1/3.0; // c = 1-h2; h2 is parent-child heritability; 0 < h2 < 1
    // static double c  = 0.5; // c = 1-h2; h2 is parent-child heritability; 0 < h2 < 1
    static final int lng = 201; // interval size
    // static final int lng = 20001; // alternative interval
    static final int mid = (lng-1)/2; // median index
    static final int sigma = mid/20;  // for the bell-curve
    static double[] pop = new double[lng]; // x-axis for initial normal bell-curve
    static int size = 100000000, repeat = 300000; // # elements & # parent-child replacements
    public static void main (String [] args) {
	Random randomno = new Random(); // get random # generator
	for (int i = 0; i < lng; i++)   // initialize the bell-curve normal distribution 
	    pop[i] = size * f(i - mid); // scale up with the size # of the elements
	// for (int i = 0; i < lng; i++) // for tracing
	//     if ( ( i % 10 ) == 0 ) System.out.println(i + " " + pop[i]);
	double sum = sum(pop); // integrate the curve
	int maxi = maxI(pop);  // get max index
	double out = out(pop); // get median
	System.out.println("sum= " + sum);  // show the sum, should be size
	System.out.println("maxI " + maxi); // show the max/median, should be mid
	System.out.println("out  " + out);  // show the expected value: is 100
	int cnt = 0; // iteration counter
	int childXparent = 0; // qualifying parents counter 
	double childXindex = 0; // sum of parent indices
	int parentXchild = 0; // qualifying children counter 
	double parentXindex = 0;// sum of children indices
	// replace parent by child repeatedly 
	while ( cnt < repeat ) {
	    double pd = (double)(mid + sigma*randomno.nextGaussian()); // select random parent
	    int p = (int) Math.round(pd); // parent
	    // System.out.print("Parent " + p);
	    if (p < 0 || lng <= p ) continue; // skip out of 0-200 range 
	    // Get a child index using the regression equation PLUS a random component that
	    // captures the random mix of the two contributing DNAs.
	    // The 0.8660282 factor is a function of the c param to preserve the sigma, x=0.5
	    // double childd = c*mid + (1-c)*p + 0.8660282*sigma*randomno.nextGaussian();
	    // The 0.7454137 factor is a function of the c param to preserve the sigma, x=1/3
	    double childd = c*mid + (1-c)*p + 0.7454137*sigma*randomno.nextGaussian();
	    int child = (int) Math.round(childd); 
	    if (child < 0 || lng <= child ) continue; // skip out of range
	    cnt++; // increment counter
	    // System.out.println(" Child: " + child);
	    pop[p] = pop[p]-1; // remove parent
	    pop[child] = pop[child]+1; // add child
	    if ( 110 == p ) { // keep track of these children
		// if ( 90 == p ) { // keep track of these children
		childXparent++; childXindex = childXindex+childd;
	    }
	    if (110 == child ) { // keep track of these parents
		// if ( 90 == child ) { // keep track of these parents
		parentXchild++; parentXindex = parentXindex + pd;
	    }
	}
	// Show the averages for the tracked children and parents
	System.out.println("child index: " + (int)(Math.round(childXindex/childXparent)));
	System.out.println("parent index: " + (int)(Math.round(parentXindex/parentXchild)));

	// Check the distribution after the parent-child replacements:
	sum = sum(pop); 
	maxi = maxI(pop);
	out = out(pop);
	System.out.println("sum= " + sum);  // the result, should be size
	System.out.println("maxI " + maxi); // show the max index, should be mid
	System.out.println("out  " + out);  // show the expected value: is 100

    } // end main

    static double f(double x) { // generator for normal distribution values 
	return (1/ (Math.sqrt(2.0*Math.PI))) * Math.exp(- (x * x/2.0)); 
    }
    static double sum(double[] pop) { // integrate the function over total interval
	return sum(pop, 0, pop.length);
    }
    static double sum(double[] pop, int p, int q) { // integrate the function over [p q)
	double sum = 0;
	for (int i = p; i < q; i++) sum += pop[i];
	return sum;
    }
    static int maxI(double[] pop) { // return the index of the max value
	int lng = pop.length;
	double max = -1;
	int maxi = -1;
	for (int i = 0; i < lng; i++) 
	    if ( max < pop[i] ) { max = pop[i]; maxi = i; }
	return maxi;
    }
    static double out(double[] pop) {
	double sum = sum(pop);
	double x = 0;
	int lng = pop.length;
	for (int i = 0; i < lng; i++) 
	    x += i*pop[i];
	return x/sum;
    }
} // end Galton
