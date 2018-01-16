#include "TimeScheme.h"
#include <string>
#include <iostream>
#include <cmath>
#include <chrono>
using namespace std;
using namespace Eigen;

int main()
{
  double t0(0.), tfinal(10.0), dt(0.001); // temps initial, final, pas de temps
  int nb_iterations = int(ceil(tfinal/dt)); // Définition du nombre d'itérations
  dt = tfinal / nb_iterations; // Recalcul de dt
  string results; // nom du fichier résultat
  int userChoiceSys;
  int userChoiceScheme;// Choix de l'utilisateur
  VectorXd sol0, exactSol; // Condition initiale et Solution exacte
  cout << "------------------------------------" << endl;
  cout << "Choississez le système : " << endl;
  cout << "1) X' = X"<< endl;
  cout << "2) x' = -y et y' = x" << endl;
  cout << "3) x' = x^2 t" << endl;
  cout << "4) Lotka-Volterra" << endl;
  cout << "5) Pendule non amorti" << endl;
  cout << "6) Pendule amorti" << endl;
  cin >> userChoiceSys;
  OdeSystem* sys(0);
  switch(userChoiceSys)
    {
    case 1:
      sys = new FirstExampleOdeSystem();
      sol0.resize(4);
      sol0(0) = 2.; sol0(1) = 3.1; sol0(2) = -5.1; sol0(3) = 0.1;
      exactSol = sol0*exp(tfinal);
      results = "first_ex";
      break;
    case 2:
      sys = new SecondExampleOdeSystem();
      sol0.resize(2); exactSol.resize(2); sol0(0) = 1; sol0(1) = -1;
      exactSol(0) = sol0(0)*cos(tfinal)-sol0(1)*sin(tfinal);
      exactSol(1) = sol0(1)*cos(tfinal)+sol0(0)*sin(tfinal);
      results = "second_ex";
      break;
    case 3:
      sys = new ThirdExampleOdeSystem();
      sol0.resize(1); exactSol.resize(1); sol0(0) = -2;
      exactSol(0) = 2.0*sol0(0)/(2.0-tfinal*tfinal*sol0(0));
      results = "third_ex";
      break;
    case 4:
      sys = new LotkaVolterraOdeSystem(1,2,3,4); // lorsque qu'on initialise à x0 = d/c et y0 = a/b le systeme n'évolue pas au cour du temps. (Voir photos pour la suite)
      // Dans le 2eme cas -> ça tourne.
      sol0.resize(2); sol0(0) = 1.; sol0(1) = 1./4.; // Conditions utilisé pour le 2eme cas
      results = "LotkaVolterra";
      break;
    case 5:
      sys = new PendulumOdeSystem(0.1,1.);
      sol0.resize(2);  sol0(0) = 3.14159265358979323/5.;  sol0(1) = 0.;
      results = "Pendulum_na";
      break;
    case 6:
      sys = new PendulumOdeSystem(0.1,1.,0.004);
      sol0.resize(2);  sol0(0) = 3.14159265358979323/5.;  sol0(1) = 0.;
      results = "Pendulum_a";
      break;
    default:
      cout << "Ce choix n'est pas possible ! Veuillez recommencer !" << endl;
      exit(0);
    }

  //-----------------------------------------------------------------------------

  cout << "Choississez le schema : " << endl;
  cout << "1) Euler explicite"<< endl;
  cout << "2) Runge-Kutta d'ordre 3" << endl;
  cout << "3) Runge-Kutta d'ordre 4" << endl;
  cout << "4) Adam Bashforth" << endl;
  cin >> userChoiceScheme;
  TimeScheme* time_scheme(0);
  switch(userChoiceScheme)
    {
    case 1 :
      time_scheme = new EulerScheme (*sys);
      results = results + "_Euler.txt";
      break;
    case 2 :
      time_scheme = new RungeKuttaScheme3 (*sys);
      results = results + "_RK3.txt";
      break;
    case 3 :
      time_scheme = new RungeKuttaScheme4 (*sys);
      results = results + "_RK4.txt";
      break;
    case 4 :
      time_scheme = new AdBashforthScheme3 (*sys);
      results = results + "_AB3.txt";
      break;
    default:
      cout << "Ce choix n'est pas possible ! Veuillez recommencer !" << endl;
      exit(0);
    }



  auto start = chrono::high_resolution_clock::now(); // début du chrono

  time_scheme->Initialize(t0, dt, sol0, results); // Initialisation
  time_scheme->SaveSolution(); // Sauvegarde condition initiale
  for (int n = 0; n < nb_iterations; n++)
    { // Boucle en temps
      time_scheme->Advance(n);
      time_scheme->SaveSolution();
    }
auto finish = chrono::high_resolution_clock::now(); // Fin du chrono
// Différence entre les deux (affichage en microsecondes demandé)
double t = chrono::duration_cast<chrono::microseconds>(finish-start).count();
cout << "Cela a pris "<< t << " microsecondes" << endl; // Affichage du résultat


  if ((userChoiceSys == 1) || (userChoiceSys == 2) || (userChoiceSys == 3))
    {
      VectorXd approxSol = time_scheme->GetIterateSolution(); // Itere au temps final
      double error = ((approxSol-exactSol).array().abs()).sum();
      cout << "Erreur = " << error<< " pour dt = " << dt << endl;

      time_scheme->Initialize(t0, dt/2., sol0, results);
      for (int n = 0; n < nb_iterations*2; n++)
	time_scheme->Advance(n);

      approxSol = time_scheme->GetIterateSolution(); // Itere au temps final
      double error2 = ((approxSol-exactSol).array().abs()).sum();
      cout << "Erreur = " << error2<< " pour dt = " << dt/2. << endl;
      cout << "Ordre de la méthode = " << log2(error/error2) << endl;
    }

  // on supprime le pointeur cree
  delete sys;
  delete time_scheme;

  return 0;
}
