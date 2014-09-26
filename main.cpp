#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;


void write(int n, char *file, double** R){
ofstream resout;
resout.open(file);
for (int i=0; i<n; i++){
    resout << i;
    for(int j=0; j<10; j++){
         resout  << "   " << R[i][j] << "   ";
    }
    resout << endl;
}
resout.close();
}


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




    // /////////////////////////////
    // /////
    // ////      print stuff
    // ///


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

//    cout << endl << endl << endl << endl << endl << "this are the first 20 eigenvalues" << endl;

   double lowest, temp;
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
          // cout << eigenvalues[i] << endl;
           }

           for(int m=0; m<n ; m++){
                      temp = R[m][i];
                      R[m][i] = R[m][k];
                      R[m][k] = temp;
                  }
   }
   cout << endl << endl;
   for(int i=0; i<n; i++){//just print first 20 eigenvectors (2)
//   cout << i << " "<< R[i][0] << endl;
   }

   cout << endl << endl << endl;

   //write the file
   write(n, "try.txt", R);



/*   //print basis
   cout << "this is the basis" << endl;
   for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){
      cout << R[i][j] << "    ";
       }
   cout <<endl << endl;
   }
*/
//end of printing stuff



   // ///////////////////////////////////////////////////
   // //////
   // /////      testing errors for lowest 3 eigenvalues
   // ////
   // ///
   // //




   double analytic1 = 3;
   double analytic2 = 7;
   double analytic3 = 11;
   double runningpmax[10];

   for(int i=3; i<13; i++){
       runningpmax[i] = i*1.0;
       cout << runningpmax[i];
   }





return 0;

}
