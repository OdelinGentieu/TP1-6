#ifndef FILE_TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur prenant en argument le system d'EDO
TimeScheme::TimeScheme(OdeSystem& sys) : _sys(sys)
{
}

// Destructeur
TimeScheme::~TimeScheme()
{
}

EulerScheme::EulerScheme(OdeSystem& sys) : TimeScheme(sys)
{
}

RungeKuttaScheme3::RungeKuttaScheme3(OdeSystem& sys) : TimeScheme(sys)
{}

RungeKuttaScheme4::RungeKuttaScheme4(OdeSystem& sys) : TimeScheme(sys)
{}

AdBashforthScheme3::AdBashforthScheme3(OdeSystem& sys) : TimeScheme(sys)
{}


// Initialisation de vos différentes variables
void TimeScheme::Initialize(double t0, double dt, VectorXd & sol0, string results)
{
  _dt = dt;
  _t = t0 ;
  _sol = sol0;
  _f.resize(_sol.rows());
  // on ouvre un fichier que si le nom est non-nul
  if (results.size() > 0)
    _sys.InitializeFileName(results);
}

// Schéma en temps par défaut : ici Euler Explicite
// Avancer d'un pas de temps
void EulerScheme::Advance(int n)
{
  _sys.BuildF(_t, _sol, _f);
  _sol += _dt*_f;
  _t += _dt;
}


void RungeKuttaScheme3::Advance(int n)
{
  _sys.BuildF(_t, _sol, _f);
  Eigen::VectorXd k1 = _f;
  _sys.BuildF (_t+(_dt/2),_sol + _dt*k1/2, _f);
  Eigen::VectorXd k2 = _f;
  _sys.BuildF (_t+_dt, _sol-_dt*k1+2*_dt*k2, _f);
  Eigen::VectorXd k3 = _f;
  _sol += (_dt/6.)*(k1+4*k2+k3);
  _t += _dt;
}



void RungeKuttaScheme4::Advance(int n)
{
  _sys.BuildF(_t, _sol, _f);
  Eigen::VectorXd k1 = _f;
  _sys.BuildF(_t+(_dt/2),_sol + _dt*k1/2, _f);
  Eigen::VectorXd k2 = _f;
  _sys.BuildF(_t+(_dt/2.), _sol+_dt*k2/2., _f);
  Eigen::VectorXd k3 = _f;
  _sys.BuildF(_t + _dt , _sol + _dt*k3, _f);
  Eigen::VectorXd k4 = _f;    
  _sol += (_dt/6.)*(k1 + 2*k2 + 2*k3 + k4);
  _t += _dt;
}



void AdBashforthScheme3::Advance(int n)
{
  if (n == 0)
    {
      _sys.BuildF(_t, _sol, _f);
      Eigen::VectorXd k1 = _f;
      _f_2 = _f;
      _sys.BuildF (_t+(_dt/2),_sol + _dt*k1/2, _f);
      Eigen::VectorXd k2 = _f;
      _sys.BuildF (_t+_dt, _sol-_dt*k1+2*_dt*k2, _f);
      Eigen::VectorXd k3 = _f;
      _sol += (_dt/6.)*(k1+4*k2+k3);
      _t += _dt;
    }
  else if (n==1)
    {
      _sys.BuildF(_t, _sol, _f);
      Eigen::VectorXd k1 = _f;
      _f_1 = _f;
      _sys.BuildF (_t+(_dt/2),_sol + _dt*k1/2, _f);
      Eigen::VectorXd k2 = _f;
      _sys.BuildF (_t+_dt, _sol-_dt*k1+2*_dt*k2, _f);
      Eigen::VectorXd k3 = _f;
      _sol += (_dt/6.)*(k1+4*k2+k3);
      _t += _dt;

    }
  else 
    {
      _sys.BuildF(_t, _sol, _f);
      _sol += (_dt/12.)*(23*_f - 16*_f_1 + 5*_f_2);
      _t += _dt;
      _f_2 = _f_1;
      _f_1 = _f;     
    }



}





// Enregistre la solution : fait appel à OdeSystem car la solution
// que l'on souhaite sauvegarder peut être différente de _sol SaveSolution
// le système
void TimeScheme::SaveSolution()
{
  _sys.SaveSolution(_t, _sol);
}

// Renvoie _sol (pratique pour calculer l'ordre de votre méthode)
const VectorXd & TimeScheme::GetIterateSolution() const
{
  return _sol;
}







#define FILE_TIME_SCHEME_CPP
#endif
