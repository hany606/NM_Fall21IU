#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>

// As we did not cover in the course how to get the eigenvalues and the eigenvectors of the matrix numerically (Topic 7: Matrix spectral problem)
//      QR algorithm
// And I tried to calculate the eignensvalues and eigenvectors symbolically (Analytically) using matlab and sympy 
//      but it was getting weird and huge results and sometimes not converging for a solution even after 30 minutes
// Thus, I chose to use eigen3 library to calculate the eigenvalues and eigenvectors
#include <eigen3/Eigen/Eigenvalues>

// Sources:
// https://www.arxiv-vanity.com/papers/1110.3397/
// https://www.codeproject.com/Articles/268589/odeint-v2-Solving-ordinary-differential-equations
// https://github.com/Thierry-Dumont/Odes-NG
// https://www.youtube.com/watch?v=KS_6mxdzQws
// https://au.mathworks.com/help/matlab/math/solve-stiff-odes.html
// https://www.wias-berlin.de/people/john/LEHRE/NUMERIK_II_20_21/ode_2.pdf
// https://www.wikiwand.com/en/Oregonator

#define eps 1e-6
#define N 6

using namespace std;

void display_mat(double *A){
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++)
            printf("%.12f\t", A[i * N + j]);
        printf("\n");
    }
    cout<<"------------------\n";
}

void display_vec(double *A){
    for (int i=0; i<N; i++){
            printf("%.12f\t", A[i]);
    }
    cout<<"\n------------------\n";

}


void set_zero(double *A){
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            A[i * N + j] = 0.0;
}
void set_eye(double *A){
    for (int i=0; i<N; i++)
        A[i * N + i] = 1.0;
}


// Function to get cofactor 
void getCofactor(double *A, double *temp, int p, int q, int n)
{
    int i = 0, j = 0;
 
    // Looping for each element of the matrix
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            // Copying into temporary matrix only those element
            // which are not in given row and column
            if (row != p && col != q)
            {
                temp[i * N + j++] = A[row * N + col];
 
                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

// Recursive function for finding determinant of matrix.
int determinant(double *A, int n)
{
    int D = 0; // Initialize result
 
    // Base case : if matrix contains single element
    if (n == 1)
        return A[0];
 
    double temp[N*N]; // To store cofactors
 
    int sign = 1; // To store sign multiplier
 
    // Iterate for each element of first row
    for (int f = 0; f < n; f++)
    {
        // Getting Cofactor of A[0][f]
        getCofactor(A, temp, 0, f, n);
        D += sign * A[0 * N + f] * determinant(temp, n - 1);
 
        // terms are to be added with alternate sign
        sign = -sign;
    }
 
    return D;
}
 
// Function to get adjoint
void adjoint(double *A,double *adj)
{
    if (N == 1)
    {
        adj[0] = 1;
        return;
    }
 
    // temp is used to store cofactors 
    int sign = 1;
    double temp[N*N];
 
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            // Get cofactor
            getCofactor(A, temp, i, j, N);
 
            // sign of adj positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;
 
            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j * N + i] = (sign)*(determinant(temp, N-1));
        }
    }
}

// Function to calculate and store inverse, returns false if
// matrix is singular
int inverse(double *inverse, double *A)
{
    // Find determinant of A[][]
    int det = determinant(A, N); 
    if (det == 0)
    {
        // printf("Singular matrix, can't find its inverse");

        return 0;
    }
 
    // Find adjoint
    double adj[N*N];
    adjoint(A, adj);
 
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            inverse[i * N + j]= adj[i * N + j]/((double) det);
 
    return 1;
}



void initialize_u_vec(double *u, double *u_0){
    for(int i = 0; i < N; i++){
        u[i] = u_0[i];
    }
}

// Trying to make it as much as possible diagonal

// a*k3 - k2*y - 2*k4*x,  a*k1,  0,  0,  0,  0,  
// 0,  -a*k1 - k2*x,  b*k5,  0,  0,  0,  
// 0,  0,  -b*k5,  2*k3*x,  0,  0,  
// k4*x,  0,  0,  -k1*y - k3*x,  0,  0,  
// 0,  0,  0,  0,  -k5*z,  0,  
// k4*x,  2*k2*x,  0,  k1*y,  0,  0,

void construct_A_mat(double *A, double *u, double *k){
    A[0] = -k[1]*u[1]+k[2]*u[3]-2*k[3]*u[0];
    A[1] = k[0]*u[3];
    A[2] = 0;
    A[3] = 0;
    A[4] = 0;
    A[5] = 0;


    A[6] = 0;
    A[7] = -k[1]*u[0]-k[0]*u[3];
    A[8] = k[4]*u[4];
    A[9] = 0;
    A[10] = 0;
    A[11] = 0;


    A[12] = 0;
    A[13] = 0;
    A[14] = -k[4]*u[4];
    A[15] = 2*k[2]*u[0];
    A[16] = 0;
    A[17] = 0;


    A[18] = k[3]*u[0];
    // A[18] = k[3]*u[0]-k[2]*u[3];
    A[19] = 0;
    A[20] = 0;
    // A[21] = -k[0]*u[1];
    A[21] = -k[0]*u[1]-k[2]*u[0];
    A[22] = 0;
    A[23] = 0;


    A[24] = 0;
    A[25] = 0;
    A[26] = 0;
    A[27] = 0;
    A[28] = -k[4]*u[2];
    A[29] =  0;


    A[30] = k[3]*u[0];
    A[31] = 2*k[1]*u[0];
    // A[31] = 2*k[1]*u[0]+k[0]*u[3];
    // A[32] = 0;
    A[33] = k[0]*u[1];
    A[33] = 0;
    A[34] = 0;
    A[35] = 0;


}

void get_eigenvalues(double *L, double *A){
    // ?
}
void get_eigenvectors(double *S, double *A){
    // ?
}

// Using Eigen3 lib:
void construct_eigenvalues(double *L, Eigen::EigenSolver<Eigen::Matrix<double, N,N> > s){
    for (int i=0; i<N; i++)
        L[i * N + i] = s.eigenvalues()(i).real();
}

void construct_eigenvectors(double *S, Eigen::EigenSolver<Eigen::Matrix<double, N,N> > s){
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            S[i * N + j] = s.eigenvectors()(i,j).real();
}

void construct_mat(double *A, Eigen::Matrix<double,N,N,Eigen::RowMajor> mat){
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            A[i * N + j] = mat(i*N+j);

}

void mat_vec_mult(double *out, double *A, double *v){
    for (int i=0; i<N; i++){
        out[i] = 0.0;
        for (int j=0; j<N; j++){
            out[i] += A[i*N+j] * v[j];
        }
    }
}
void mat_mat_mult(double *matrix_1, double *matrix_2, double *matrix_product) {

    for(int k = 0; k < N; k++){
        matrix_product[k] = 0;
            for (int j = 0; j<N; j++) {
                matrix_product[k] += matrix_1[k*N+j]* matrix_2[j];
            }
    }

}



void update(double *u, double *L, double *S, double *S_inv, double tau){
    // u^{n+1} = S@e^{\tau\Lambda}@S^{-1}@u^n
    // Formed to:
    // v^n = S^{-1}u^n
    double *v_vec= (double *) malloc(N*sizeof(double)); // To store 4xn @ nx4
    mat_vec_mult(v_vec, S_inv, u);

    Eigen::Map<Eigen::Matrix<double,N,N,Eigen::RowMajor> > tmp(L);
    // tau = tmp.determinant();
    // cout<<"Tau: "<< tau<<endl;
    // v^{n+1}_i = e^{\tau\lambda_i)}v^n_i
    for (int i=0; i<N; i++){
        // cout<<v_vec[i]<<"\t";
        v_vec[i] *= exp(tau*L[i*N+i]);
        // cout<<v_vec[i]<<"\n";

    }

    //-------------------------------------------
    // Debug:
    // double *test= (double *) malloc(N*N*sizeof(double));
    // set_zero(S_inv);
    // S @ S^{-1} != zero_mat ??? Because it is an approximation for the inverse
    // mat_mat_mult(S, S_inv, test);
    // cout<<"^^^^^^^^^^^^^^^\n";
    // display_mat(test);
    // cout<<"^^^^^^^^^^^^^^^\n";
    //-------------------------------------------

    // And again go back from v to u
    // u^{n+1} = Sv^{n+1}
    mat_vec_mult(u, S, v_vec);

    // Put v in u and return it
    // double *eye= (double *) malloc(N*N*sizeof(double)); // To store 4xn @ nx4
    // set_eye(eye);
    // mat_vec_mult(u, eye, v_vec);

    free(v_vec);
}


int main() {
    string line;
    ifstream in_file;
    int T;
    double initial_concentration[6];
    double rate_const[5];

    stringstream ss;
    stringstream ss2;   // Bad way to fix the problem 


    // istream ret_flag;
    in_file.open("input.txt");
    getline(in_file, line);
    T = stoi(line);
    
    getline(in_file, line);
    ss.str(line);
    // ss<< line;

    for(int i = 0; i < 6; i++)
        ss >> initial_concentration[i];

    getline(in_file, line);
    // ss.str(line);
    // ss<< line;
    ss2.str(line);

    for(int i = 0; i < 5; i++)
        ss2 >> rate_const[i];

    cout<< "T: "<< T <<endl;
    cout<<"Initial Concentrations: ";
    for(int i = 0; i < 6; i++){
        cout<<initial_concentration[i]<<" ";
    }
    cout<<"\nRate constants: ";

    for(int i = 0; i < 5; i++){
        cout<<rate_const[i]<<" ";
    }
    cout<<"\n";
    in_file.close();

    double *A_mat= (double *) malloc(N*N*sizeof(double)); // To store 4xn @ nx4
    double *u_vec= (double *) malloc(N*sizeof(double)); // To store 4xn @ nx4
    double *eigenvalues= (double *) malloc(N*N*sizeof(double)); // To store 4xn @ nx4
    double *eigenvectors= (double *) malloc(N*N*sizeof(double)); // To store 4xn @ nx4
    double *inv_eigenvectors= (double *) malloc(N*N*sizeof(double)); // To store 4xn @ nx4

    set_zero(A_mat);
    set_zero(eigenvalues);
    set_zero(eigenvectors);
    set_zero(inv_eigenvectors);

    initialize_u_vec(u_vec, initial_concentration);
    cout<<"------------------\n";
    for(int i = 0; i < T; i++){
        // if(i > 0){    
        //     // And again go back from v to u
        //     // u^{n+1} = Sv^{n+1}
        //     mat_vec_mult(u, S, v_vec);

        // }
        // Construct A matrix
        construct_A_mat(A_mat, u_vec, rate_const);
        // Eigen::Matrix<double, N, N> A_eig;
        // A_eig << A_mat[0], A_mat[1], A_mat[2], A_mat[3], A_mat[4], A_mat[5], A_mat[6], A_mat[7], A_mat[8], A_mat[9], A_mat[10], A_mat[11], A_mat[12], A_mat[13], A_mat[14], A_mat[15], A_mat[16], A_mat[17], A_mat[18], A_mat[19], A_mat[20], A_mat[21], A_mat[22], A_mat[23], A_mat[24], A_mat[25], A_mat[26], A_mat[27], A_mat[28], A_mat[29], A_mat[30], A_mat[31], A_mat[32], A_mat[33], A_mat[34], A_mat[35];
        Eigen::Map<Eigen::Matrix<double,N,N,Eigen::RowMajor> > A_eig(A_mat);
        Eigen::EigenSolver<Eigen::Matrix<double, N,N> > s(A_eig); // the instance s(A) includes the eigensystem
        // std::cout << A_eig << std::endl;
        // std::cout << "eigenvalues:" << std::endl;
        // std::cout << s.eigenvalues() << std::endl;
        // std::cout << "eigenvectors=" << std::endl;
        // std::cout << s.eigenvectors() << std::endl;

        // display_mat(A_mat);
        // Get eigenvalues
        // Construct eignensvalues matrix \Lambda
        // get_eigenvalues(eigenvalues, A_mat);   
        construct_eigenvalues(eigenvalues, s);
        // cout<<"Eigenvalues matrix: \n";
        // display_mat(eigenvalues);

        // Check stability
        // Find eigenvectors
        // Construct eigenvectors matrix S
        // get_eigenvectors(eigenvectors, A_mat);
        construct_eigenvectors(eigenvectors, s);
        // cout<<"Eigenvectors matrix: \n";

        // display_mat(eigenvectors);

        Eigen::Map<Eigen::Matrix<double,N,N,Eigen::RowMajor> > inv_mat(eigenvectors);
        inv_mat = inv_mat.inverse().eval();
        // std::cout << "inv:" << std::endl;
        // std::cout << inv_mat << std::endl;

        construct_mat(inv_eigenvectors, inv_mat);
        // display_mat(inv_eigenvectors);

        // Get the inverse matrix of the previous matrix
        // cout<< inverse(inv_eigenvectors, eigenvectors)<< endl;
        // cout<<"Inverse Eigenvectors matrix: \n";

        // display_mat(inv_eigenvectors);
        // Compute tau (step)
        double tau = 0.00001;//1/abs(s.eigenvalues()(0).real()); //0.00001;
        // Update step
        // Get v_vec inside u_vec
        update(u_vec, eigenvalues, eigenvectors, inv_eigenvectors, tau);
    
        display_vec(u_vec);
        // break;

    }

    // system of equations
    // u vector
    // f vector
    // Construct A matrix

    // Get eigenvalues
    
    // Construct eignensvalues matrix \Lambda
    // Check stability
    // Find eigenvectors
    // Construct eigenvectors matrix S
    // Get the inverse matrix of the previous matrix
    // u^{n+1} = S@e^{\tau\Lambda}@S^{-1}@u^n

    ofstream out_file;
    out_file.open("output.txt");

    for(int i = 0; i < N-1; i++){
        out_file << std::setprecision(10) << std::scientific <<  u_vec[i]<<" ";
    }

    out_file << std::setprecision(10) << std::scientific << u_vec[N-1];

    out_file.close();
    
    free(A_mat);
    free(u_vec);
    free(eigenvalues);
    free(eigenvectors);
    free(inv_eigenvectors);

}