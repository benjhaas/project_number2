#include <iostream>
#include <cmath>

using namespace std;


int main()
{

   //still to do: put h2, counter for iterations

    //dimension d
    //number of steps n
    //steplengh h

    double tolerance = 10^-8;
    int dimension, n;

    cout << "Type in dimension: ";
    cin >> dimension;
    n=dimension;


    double pmax=5; // later cin????
    double h=pmax/n;

    double off; // sum over all offdiagonal elements
    double maxoff; // maximum of all offdiagonal elements
    //double maxoffde; // MaximumOFFDiagonalElement
    int l, k; // position of largest element
    double tau, t, c, s; //tan(tetha), cos(theta), sin(theta)







    //define matrix
   double** matrix = new double*[n];
   for(int i=0; i<n; i++){
   matrix[i] = new double[n];
   }

   for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){
         if(i==j){matrix[i][j]=2./(h*h)+i*i*h*h;}
         if(i!=j){matrix[i][j]=0;}
         if(i==j+1||j==i+1){matrix[i][j]= -1./(h*h);}
      cout << matrix[i][j] << "    ";
       }
   cout <<endl << endl;
   } //end of matrix definition


   // ///////////////////////////
   // ///////
   // /////    diagonalize matrix
   // ////     with jacobi
   // ///


   maxoff = 0; // calculate first value for maximum of offdiagonal elements
   for(int i=0; i<n; i++){
       for(int j=0; j<n; j++){
           if(i!=j){
               if(matrix[i][j] > maxoff){
                   maxoff = matrix[i][j];
                   l = i;
                   k = j;
               }
           }
       }
   } //end of calculating first maxoff
   cout << "first maxoff= " << maxoff << endl;

// diagonalize till maximum is smaller than tolerance
   while (maxoff > tolerance){
       //first part
       if(matrix[k][l] != 0){
          tau = (matrix[l][l] - matrix[k][k]) / (2*matrix[k][l]);
          if(tau > 0){
             t = 1.0/(tau + sqrt(1.0 + tau*tau)); //  where from???!!!!
          }
          else{
             t = -1.0/(-tau + sqrt(1.0 + tau*tau)); // where from?????!!!
          }
          c = 1/(sqrt(1+t*t));
          s = c*t;
       }
       else{
          c = 1.0;
          s = 0.0;
       }
       // calculating new k,l matrix elements
       double helpkk = matrix[k][k];
       double helpll = matrix[l][l];
       matrix[k][k] = helpkk*c*c - 2*c*s*matrix[k][l]+s*s*helpll;
       matrix[l][l] = helpll*c*c - 2*c*s*matrix[k][l]+s*s*helpkk;
       // calculate new non k,l matrix elememts
       for(int i=0; i<n; i++){
           if(i!=k && i!= l){
               double helpik = matrix[i][k];
               double helpil =  matrix[i][l];
               matrix[i][k] = helpik*c - helpil*s;
               matrix[k][i] = matrix[i][k];
               matrix[i][l] = helpil*c + helpik*s;
               matrix[l][i] = matrix[i][l];
           }
       }
       //end first part

       maxoff = 0; // calculate the next value for maximum of offdiagonal elements
       for(int i=0; i<n; i++){
           for(int j=0; j<n; j++){
               if(i!=j){
                   if(matrix[i][j] > maxoff){
                       maxoff = matrix[i][j];
                       l = i;
                       k = j;
                   }
               }
           }
       } //end of calculating next maxoff
   }

// diagonalize till sum over offd. elements is smaller than tolerance
   off = 0; // calculate first value for off
   for(int i=0; i<n; i++){
       for(int j=0; j<n; j++){
           if(j!=i){
               off = sqrt(off + matrix[i][j] * matrix[i][j]);
           }
       }
   } // end of calculating first off
   while(off > tolerance){

       //second part
       if(matrix[k][l] != 0){
          tau = (matrix[l][l] - matrix[k][k]) / (2*matrix[k][l]);
          if(tau > 0){
             t = 1.0/(tau + sqrt(1.0 + tau*tau)); //  where from???!!!!
          }
          else{
             t = -1.0/(-tau + sqrt(1.0 + tau*tau)); // where from?????!!!
          }
          c = 1/(sqrt(1+t*t));
          s = c*t;
       }
       else{
          c = 1.0;
          s = 0.0;
       }
       // calculating new k,l matrix elements
       double helpkk = matrix[k][k];
       double helpll = matrix[l][l];
       matrix[k][k] = helpkk*c*c - 2*c*s*matrix[k][l]+s*s*helpll;
       matrix[l][l] = helpll*c*c - 2*c*s*matrix[k][l]+s*s*helpkk;
       // calculate new non k,l matrix elememts
       for(int i=0; i<n; i++){
           if(i!=k && i!= l){
               double helpik = matrix[i][k];
               double helpil =  matrix[i][l];
               matrix[i][k] = helpik*c - helpil*s;
               matrix[k][i] = matrix[i][k];
               matrix[i][l] = helpil*c + helpik*s;
               matrix[l][i] = matrix[i][l];
           }
       }
       //end second part


       off = 0; // calculate next value for off
       for(int i=0; i<n; i++){
           for(int j=0; j<n; j++){
               if(j!=i){
                   off = sqrt(off + matrix[i][j] * matrix[i][j]);
               }
           }
       } // end of calculating next off

   }//end of matrix diagonalization
   for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){
      cout << matrix[i][j] << "    ";
       }
   cout <<endl << endl;
   }
return 0;

}






