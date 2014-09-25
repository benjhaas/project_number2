#include <iostream>
#include <cmath>

using namespace std;


int main()
{

   //still to do: put h2, counter for iterations

    //dimension d
    //number of steps n
    //steplengh h

    double tolerance = pow(10, -8);
    cout << "The tolerance is " << tolerance << endl;
    int dimension, n;

    cout << "Type in dimension: ";
    cin >> dimension;
    n=dimension;


    double pmax = 5.0; // later cin????
    double h=pmax/double(n);

//    double off; // sum over all offdiagonal elements
    double maxoff; // maximum of all offdiagonal elements
    //double maxoffde; // MaximumOFFDiagonalElement
    int l, k; // position of largest element
    double tau, t, c, s; //tan(tetha), cos(theta), sin(theta)
    double rik, ril; // temps for eigenvectors






    //define main matrix
   double** matrix = new double*[n];
   for(int i=0; i<n; i++){
   matrix[i] = new double[n];
   }

   for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){
         if(i==j){matrix[i][j]=2./(h*h)+i*i*h*h;}
         if(i!=j){matrix[i][j]=0;}
         if(i==j+1||j==i+1){matrix[i][j]= -1./(h*h);}
     // cout << matrix[i][j] << "    ";
       }
  // cout <<endl << endl;
   } //end of matrix definition

//cout << endl << endl << endl;

   // define unity matrix
   double** R = new double*[n];
   for(int i=0; i<n; i++){
       R[i] = new double[n];
   }
   for(int i=0; i<n; i++){
       for(int j=0; j<n; j++){
           if(i==j){
               R[i][j] = 1;
           } else {
               R[i][j] = 0;
           } //cout << R[i][j] << "    ";
       } //cout << endl << endl;
   } // end of unity matrix definition


   // ///////////////////////////
   // ///////
   // /////    diagonalize matrix
   // ////     with jacobi
   // ///


   maxoff = -1./(h*h) ; // set initial values
   k = 0;
   l = 1;

   //cout << "first maxoff= " << maxoff << endl;

// diagonalize till maximum is smaller than tolerance
   while (fabs(maxoff) > tolerance){
       //first part
      // if(matrix[k][l] != 0){
          tau = (matrix[l][l] - matrix[k][k]) / (2*matrix[k][l]);        
          if(tau > 0){
            t = 1.0 / (tau + sqrt(1.0 + tau*tau));
          }
          else{
            t = -1.0 /(-tau + sqrt(1.0 + tau*tau));
          }
          c = 1/(sqrt(1+t*t));
          s = c*t;

    //  }else{
      //    c = 1.0;
        //  s = 0.0;
       //}
       // calculating new k,l matrix elements
       double helpkk = matrix[k][k];
       double helpll = matrix[l][l];
       matrix[k][k] = helpkk*c*c - 2*c*s*matrix[k][l]+s*s*helpll;
       matrix[l][l] = helpkk*s*s + 2*c*s*matrix[k][l]+c*c*helpll;
       matrix[l][k]= 0.0;
       matrix[k][l]= 0.0;
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
           rik = R[i][k];
           ril = R[i][l];
           R[i][k] = c*rik - s*ril;
           R[i][l] = c*ril + s*rik;
       }
       //end first part

       maxoff = 0.0; // calculate the next value for maximum of offdiagonal elements
       for(int i=0; i<n; i++){
           for(int j=i+1; j<n; j++){
                   if(fabs(matrix[i][j]) > fabs(maxoff)){
                       maxoff = matrix[i][j];
                       k = i;
                       l = j;
                   }
           }
          // cout << maxoff << endl;
       } //end of calculating next maxoff
   }




   //print stuff

/*   cout << endl << endl;
   //print diagonalized matrix
   cout << "this is the diagonalized matrix:" << endl;
   for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){
      cout << matrix[i][j] << "    ";
       }
   cout <<endl << endl;
   }
*/
  double eigenvalues[n];
   //print eigenvalues
   for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){
          if(i==j){
        //  cout << matrix[i][j] << endl;
          eigenvalues[i] = matrix[i][j];
          //cout << sortedelemnts[i] << endl;
          }
       }
   }

    cout << endl << endl << endl << endl << endl << "this are the first 20 eigenvalues" << endl;

   double lowest, temp;
   double lowestv, tempv;
   for(int i=0; i<n ; i++){
           lowest = temp = eigenvalues[i];
           int k = i;
           for(int j=i+1; j<n; j++){
               if( fabs(eigenvalues[j]) < fabs(lowest) ){
               lowest = eigenvalues[j];
               k = j;
               }
           }
           eigenvalues[i] = lowest;
           eigenvalues[k] = temp;



           if(i < 20){
           cout << eigenvalues[i] << endl;
           }

           for(int m=0; m<n ; m++){
                      temp = R[m][i];
                      R[m][i] = R[m][k];
                      R[m][k] = temp;
                  }
   }
   cout << endl << endl;
   for(int i=0; i<n; i++){//just print first 20 eigenvectors (2)
   cout << i << " "<< R[i][0] << endl;
   }

   cout << endl << endl << endl;

/*   //print basis
   cout << "this is the basis" << endl;
   for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){
      cout << R[i][j] << "    ";
       }
   cout <<endl << endl;
   }
*/





 /*  //crap:
   double help = 0;
   double position;
   for(int i=0; i<n; i++){
       for(int j=i; j<n; j++){
           if(sortedelemnts[j] >= sortedelemnts[i-1])
           help = sortedelemnts[j];
           position = j;
       }
   } //end of crap
*/



   //end of print out

/*
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
   */





return 0;

}
