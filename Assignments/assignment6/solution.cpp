#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


// TODO:
// Check: https://www.arxiv-vanity.com/papers/1110.3397/
// Check: https://www.codeproject.com/Articles/268589/odeint-v2-Solving-ordinary-differential-equations
// Check: https://github.com/Thierry-Dumont/Odes-NG
// https://www.youtube.com/watch?v=KS_6mxdzQws
// https://au.mathworks.com/help/matlab/math/solve-stiff-odes.html
// https://www.wias-berlin.de/people/john/LEHRE/NUMERIK_II_20_21/ode_2.pdf
// https://www.wikiwand.com/en/Oregonator
// https://www.wikiwand.com/en/Oregonator

#define eps 1e-6

using namespace std;

int main() {
    string line;
    ifstream in_file;
    int T;
    float initial_concentration[6];
    float rate_const[5];

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



}