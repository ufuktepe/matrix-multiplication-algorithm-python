# *Strassen's Divide and Conquer Matrix Multiplication Algorithm*

Implemented the conventional and Strassen's matrix multiplication algorithms for 𝑛 × 𝑛 matrices and determined the optimal cross-over point both analytically and experimentally. For 𝑛 × 𝑛 matrices, the cross-over point between the two algorithms is the value of 𝑛 for which we stop using Strassen's algorithm and switch to conventional matrix multiplication.

## **Data Layout Optimization**

Splitting matrices takes up a significant portion of the actual runtime in Strassen’s algorithm. To speed up this process, instead of using a standard row-major ordering, Morton ordering is used to represent 𝑛 × 𝑛 matrices. Morton ordering takes a 2D array stored in row-major order and arranges the matrix in 𝑚 × 𝑚 block arrays where 𝑚 is the size of the “base case”.

## **Padding**

Strassen’s algorithm recursively divides 𝑛 × 𝑛 matrices into four 𝑛/2 × 𝑛/2 matrices. That means at each recursive call 𝑛 must be divisible by 2. The obvious method to resolve this issue would be to pad the original matrix to the next power of 2 with zero rows and zero columns. However, this would be an expensive approach because we would have to almost double the size of the original matrix if its original size is just over a power of 2 (i.e 𝑛=1025).
Instead, we use the following approach: we first find the size of the “base case” for our Morton-ordered array. Then we keep doubling that size until we reach or exceed the size of our original matrix to find the minimum required size that we need to perform Strassen’s algorithm. That is:
1. Let 𝑝=𝑛. Repeatedly divide 𝑝 in half, each time taking the ceiling, until 𝑝 is less than or equal to the cross-over point. This will be equal to the size of the “base case” for our Morton-ordered array.
2. Then repeatedly double 𝑝 until 𝑝≥𝑛.
3. Pad the original matrix with zero rows and zero columns until we obtain a matrix of size 𝑝 ×𝑝.
In other words, we are padding the matrix just enough so that Strassen’s algorithm can reach the “base case”. One of the main advantages of this approach is, since the padding is done upfront, we do not have to worry about odd-sized matrices while performing Strassen’s algorithm.

## **Multiprocessing**

To reduce the overall runtime, Python’s multiprocessing module is utilized while performing the experimental analysis to find the optimum cross-over point.
