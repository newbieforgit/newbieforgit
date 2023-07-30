#include<iostream>
#include<ctime>
#include<cstdlib>
#include<fstream>
#include<cmath>
#include<vector>
using namespace std; 

int main () { 
    ofstream output;
    output.open("TvsM.txt");
    const int N = 20;
    const double J = -1.0;
    const int M = 1000000;
    double M_avg = 0.0;
    const int sample = 10000;
    vector<double> T_list ;
    for (int i=0; i < 21; i++) { 
        T_list.push_back ((5.*i)/21);
        //cout << T_list[i] << endl;
    }
    srand(time(NULL));
    for (int i =0; i < T_list.size(); i++) { 
        //cout << T_list[i];
        int S[N][N];
        for (int k =0; k < N; k++){
            for (int j =0; j< N; j++) {
                double term = (rand()/(double)RAND_MAX)-0.5;
                if (term < 0) {
                    S[k][j] = -1;
                }
                else {
                    S[k][j] = 1;
                }
            //cout << S[i][j] << endl;
            }
        } //set magnetic moment  
        for (int k=0; k< M; k++ ){
            int ifl = rand() % (N ) ; 
            int jfl = rand() % (N ); //set flip index

            double Eflip = (S[ifl][jfl-1] + S[ifl][(jfl+1)%N] + S[(ifl+1)%N][jfl-1] + S[(ifl+1)%N][jfl] + S[ifl-1][jfl] + S[ifl-1][(jfl+1)%N]) * (-2) * J * S[ifl][jfl];
            double p = rand()/(double)RAND_MAX;
            if (p < exp(-Eflip /(T_list[i]))) {
                S[ifl][jfl] *= (-1);
            }
            if (i%sample ==0) {
                int sum =0;
                for (int l =0; l<N; l++) {
                    for (int j=0; j<N; j++) {
                        sum+= S[l][j];
                    }
                }
                double result = sqrt( (sum*sum)/(N*N)*(sum*sum)/(N*N)) ;
                M_avg += result;
                 }
        }
        M_avg /= (M / sample);
        if (output.is_open()){
            output << T_list[i] << "\t" << M_avg << endl;
        }
        
        }
        return 0;

 }