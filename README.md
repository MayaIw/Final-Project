This is the final project of "Software Project" course, in Tel Aviv University. </br>
In this project, we implemented a clustering algorithm that is based on symmetric Non-negative Matrix Factorization (symNMF), and compared it to Kmeans. </br>
Further explanation on symNMF: Da Kuang, Chris Ding, and Haesun Park. Symmetric nonnegative matrix factorization for graph clustering. </br>
In Proceedings of the 2012 SIAM International Conference on Data Mining (SDM), Proceedings, pages 106–117. Society for Industrial and Applied Mathematics, April 2012. 

The project consists of the following files:
1. symnmf.py: Python interface of our code.
2. symnmf.h: C header file.
3. symnmf.c: C interface of our code.
4. symnmfmodule.c: Python C API wrapper.
5. analysis.py: Analyze the algorithm- compare it to Kmeans.
6. setup.py: The setup file.
7. Makefile: make script to build the C interface.

Python Program (symnmf.py) </br>
Reading user CMD arguments: </br>
(a) k (int, < N): Number of required clusters. </br>
(b) goal: Can get the following values: </br> 
  i. symnmf: Performs full symNMF and outputs the decomposition matrix H. </br>
  ii. sym: Calculates and outputs the similarity matrix. </br>
  iii. ddg: Calculates and outputs the Diagonal Degree Matrix. </br>
  iv. norm: Calculates and outputs the normalized similarity matrix. </br>
(c) file_name (.txt): The path to the Input file, it will contain N data points for all above goals, the file extension is .txt. </br>
Build the extension by running the following command: python3 setup.py build_ext --inplace

C Program (symnmf.c) </br>
Reading user CMD arguments: </br>
(a) goal: Can get the following values: </br>
  i. sym: Calculates and outputs the similarity matrix. </br>
  ii. ddg: Calculates and outputs the Diagonal Degree Matrix. </br>
  iii. norm: Calculates and outputs the normalized similarity matrix. </br>
(b) file_name (.txt): The path to the Input file, it will contain N data points for all above goals, the file extension is .txt. </br>
Compile the C program by running the following command: make 

analysis.py </br>
Compares SymNMF to Kmeans. Applies both methods to a given dataset and reports the silhouette_score from the sklearn.metrics. </br>
Reading user CMD arguments: </br>
(a) k (int, < N): Number of required clusters. </br>
(b) file_name (.txt): The path to the Input file, it will contain N data points, the file extension is .txt.
