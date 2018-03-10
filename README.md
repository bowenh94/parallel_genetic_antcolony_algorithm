# parallel_genetic_antcolony_algorithm

*parallel hybrid of genetic antcolony algorithm to solve TSP problem. 

*code files: ga-ac-seq.c ; ga-ac-parallel.c

*ga-ac-seq.c:
	**This code aims to sequentially solve the program.

	**The running process of this code can roughly divide into two parts, the first part is using GA algorithm to calculate a set of solutions and initializing pheromone matrix with these solutions, the next part is using AC algorithm to get the best solution.

	**Running this code:
	  Generally, you can run compile.sh and run_seq.sh to test our program. However, you can also test our program by changing some variables: 
	  1) you can change variable ’N’(the number of input cities).
	  2) you can change input file name.
	  3) you can change variables ‘MAX_GEN’(the number of population);
			      	      ‘NUM_ANT’(the number of ants);
			      	      ‘RC_ITER’(iteration times);
	     these variables effect running time.

*ga-ac-parallel.c:
	**This code aims to achieve the program in parallel.
	
	**The running process of this code is similar with the former one, but it runs those two parts in parallel.

	**Running this code:
	  Generally, you can run compile.sh and run_par.sh to test our program. However, you can also test our program by changing some variables:
	  1) you can change variable ’N’(the number of input cities).
	  2) you can change input file name.
	  3) you can change variables ‘MAX_GEN’(the number of population);
			      	      ‘NUM_ANT’(the number of ants);
			      	      ‘RC_ITER’(iteration times);
	     these variables effect running time.
