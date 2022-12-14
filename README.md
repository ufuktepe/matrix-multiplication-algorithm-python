# *Strassen's Divide and Conquer Matrix Multiplication Algorithm*

Implemented the conventional and Strassen's matrix multiplication algorithms for π Γ π matrices and determined the optimal cross-over point both analytically and experimentally. For π Γ π matrices, the cross-over point between the two algorithms is the value of π for which we stop using Strassen's algorithm and switch to conventional matrix multiplication.

## **Data Layout Optimization**

Splitting matrices takes up a significant portion of the actual runtime in Strassenβs algorithm. To speed up this process, instead of using a standard row-major ordering, Morton ordering is used to represent π Γ π matrices. Morton ordering takes a 2D array stored in row-major order and arranges the matrix in π Γ π block arrays where π is the size of the βbase caseβ.

## **Padding**

Strassenβs algorithm recursively divides π Γ π matrices into four π/2 Γ π/2 matrices. That means at each recursive call π must be divisible by 2. The obvious method to resolve this issue would be to pad the original matrix to the next power of 2 with zero rows and zero columns. However, this would be an expensive approach because we would have to almost double the size of the original matrix if its original size is just over a power of 2 (i.e π=1025).
Instead, we use the following approach: we first find the size of the βbase caseβ for our Morton-ordered array. Then we keep doubling that size until we reach or exceed the size of our original matrix to find the minimum required size that we need to perform Strassenβs algorithm. That is:
1. Let π=π. Repeatedly divide π in half, each time taking the ceiling, until π is less than or equal to the cross-over point. This will be equal to the size of the βbase caseβ for our Morton-ordered array.
2. Then repeatedly double π until πβ₯π.
3. Pad the original matrix with zero rows and zero columns until we obtain a matrix of size π Γπ.
In other words, we are padding the matrix just enough so that Strassenβs algorithm can reach the βbase caseβ. One of the main advantages of this approach is, since the padding is done upfront, we do not have to worry about odd-sized matrices while performing Strassenβs algorithm.

## **Multiprocessing**

To reduce the overall runtime, Pythonβs multiprocessing module is utilized while performing the experimental analysis to find the optimum cross-over point.
