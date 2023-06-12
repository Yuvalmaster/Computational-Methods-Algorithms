# Numerical Computational Methods Algorithms Implementation and Performance Evaluation 🏫👩‍🏫📚🔢📐
This repository showcases a project that focuses on the implementation and performance evaluation of various numerical computational methods algorithms. The project aims to analyze and compare the efficiency, accuracy, and convergence properties of different computational methods commonly used in scientific and engineering applications.

## Project Overview
The main objective of this project is to provide a comprehensive analysis of numerical algorithms used in scientific computing. By implementing and evaluating a diverse set of algorithms, the project aims to assess their performance in solving simple mathematical problems. The evaluation includes comparing their efficiency, accuracy, and convergence properties to identify the most suitable algorithm for different scenarios.

### Included Algorithms
The project includes the following numerical algorithms:

1. AdamBashforth Algorithm: A numerical integration method for solving ordinary differential equations (ODEs) using an explicit multi-step approach that approximates the solution based on previous solution values and their derivatives.

2. Bisection Algorithm: A root-finding method used to find the roots of a continuous function by iteratively bisecting an interval and narrowing down the range until the desired accuracy is achieved.

3. Fixed Point Algorithm: An iterative method used to find the fixed points of a given function by repeatedly applying the function to an initial guess until convergence to a fixed point is achieved.

4. Gauss Seidel Algorithm: An iterative method for solving systems of linear equations that updates the solution iteratively by considering the most recently computed values, providing a more efficient approach compared to the Gauss-Jordan elimination method.

5. Henstenes Steifel Algorithm: An iterative method for solving systems of linear equations that combines ideas from the Gauss Seidel and conjugate gradient methods to achieve faster convergence and improved stability.

6. Jacobi Algorithm: An iterative method for solving systems of linear equations that updates the solution iteratively by using the previous iteration's values, making it a simple and widely used method for solving large linear systems.

7. Newton Raphson Algorithm: A root-finding method used to find the roots of a differentiable function by iteratively refining an initial guess using the function's value and derivative, converging to the root with quadratic convergence rate under certain conditions.

8. Predictor Corrector Algorithm: A numerical integration method for solving ordinary differential equations that combines the predictions from a lower-order method (predictor) and the corrections from a higher-order method (corrector) to improve accuracy and stability.

9. QR Factorization Algorithm: A numerical method for decomposing a matrix into the product of an orthogonal matrix and an upper triangular matrix. It is widely used in solving linear least squares problems and eigenvalue computations.

10. Runge Kutta Algorithm: A numerical method for solving ordinary differential equations using an explicit multi-stage approach that estimates the solution at each time step by considering weighted averages of function evaluations at multiple points within the time interval.

11. SOR (Successive Over-Relaxation) Algorithm: An iterative method for solving systems of linear equations that improves upon the Gauss Seidel method by introducing an additional relaxation parameter to accelerate convergence.

12. Secant Algorithm: A root-finding method used to find the roots of a continuous function by employing a secant line between two points instead of bisecting intervals, similar to the Bisection algorithm.

13. Taylor Algorithm: A numerical method for approximating functions by expanding them into a Taylor series, providing a polynomial representation that enables accurate approximations within a given range of inputs.

## Implementation Details
The entire project is implemented in MATLAB. The MATLAB code for each algorithm is organized within the repository's folder structure, with the necessary functions and scripts contained in the Functions folder.

The implementation of each algorithm follows established mathematical principles and algorithms specific to their respective domains.


Feel free to explore the repository, review the implemented algorithms, and analyze the performance evaluation results. The project provides valuable insights into the efficiency, accuracy, and convergence properties of various numerical methods commonly employed in scientific and engineering applications.
