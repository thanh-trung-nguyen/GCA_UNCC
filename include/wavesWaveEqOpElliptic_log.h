// Copyright (C) 2000 Larisa Beilina
//
// This file is part of WavES project.
//
// WavES is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with WavES. If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2000-01-01 by Larisa Beilina
// Last changed: 2012-05-29 by Larisa Beilina



#ifndef _WAVEEQOPELLIPTIC_log_
#define _WAVEEQOPELLIPTIC_log_

#include "wavesSDOperator.h"

/* 
   Here, the operators needed for application 2 is collected.

   */

// In the interior, y_i := 2dt laplace x_i - y_i

/**
  WavesWaveEqInterior is an operator for the scalar wave equation.

  The scalar wave equation is discretized with centred differences. The {\tt SDField y} represnting the unknown on time level {\tt n-1} is overwritten to represent values on time level {\tt n+1}, using {\tt SDField x} representing the unknown on time level {\tt n}, according to 
\begin{verbatim}
  y_i := 2dt laplace x_i - y_i
\end{verbatim}
*/

class WavesWaveEqInteriorElliptic : public WavesSDOperator 
{
public:
  /// SDindexes covering the interior is created. The time step dt shall be supplied.
  WavesWaveEqInteriorElliptic(WavesSDGeometry *sdg); 

  /// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
  WavesWaveEqInteriorElliptic(WavesSDIndexes  *sdi);

  /// Destructor frees memory if needed.
  ~WavesWaveEqInteriorElliptic() {
    if(freeIndexes)
      delete indexes;
  }


  bool applyWaveElliptic(real *x, real *y, 
			 const real &omega ) {
    
      omega_ = omega;
      
// cout<<"dt in  applyPlaneWaveOpexpl"<< step<<endl; 
      return doApply(x, y);
  }
  



 protected:
  ///
  bool doApply(const double *x, double *y) const ;
  

  /// Utility
  void doInitialize();

 private:
  double omega_;
  //
  bool freeIndexes;
};

//===============================================================
class WavesWaveEqInteriorLaplace : public WavesSDOperator 
{
public:
  /// SDindexes covering the interior is created. The time step dt shall be supplied.
  WavesWaveEqInteriorLaplace(WavesSDGeometry *sdg); 

  /// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
  WavesWaveEqInteriorLaplace(WavesSDIndexes  *sdi);

  /// Destructor frees memory if needed.
  ~WavesWaveEqInteriorLaplace() {
    if(freeIndexes)
      delete indexes;
  }


  bool applyWaveElliptic(real *x, real *y ) {
    
    
      
// cout<<"dt in  applyPlaneWaveOpexpl"<< step<<endl; 
      return doApply(x, y);
  }
  



 protected:
  ///
  bool doApply(const double *x, double *y) const ;
  

  /// Utility
  void doInitialize();

 private:
 
  //
  bool freeIndexes;
};

//=========================================

class WavesWaveEqInteriorElliptic_laplace : public WavesSDOperator 
{
public:
  /// SDindexes covering the interior is created. The time step dt shall be supplied.
  WavesWaveEqInteriorElliptic_laplace(WavesSDGeometry *sdg); 

  /// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
  WavesWaveEqInteriorElliptic_laplace(WavesSDIndexes  *sdi);

  /// Destructor frees memory if needed.
  ~WavesWaveEqInteriorElliptic_laplace() {
    if(freeIndexes)
      delete indexes;
  }


  bool applyWaveElliptic(real *x, real *y, 
			 const real &omega ) {
    
      omega_ = omega;
      
// cout<<"dt in  applyPlaneWaveOpexpl"<< step<<endl; 
      return doApply(x, y);
  }
  



 protected:
  ///
  bool doApply(const double *x, double *y) const ;
  

  /// Utility
  void doInitialize();

 private:
  double omega_;
  //
  bool freeIndexes;
};



class WavesWaveEqBoundaryRight : public WavesSDOperator 
{
public:
 
  /// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
  WavesWaveEqBoundaryRight(WavesSDIndexes  *sdi);

  /// Destructor frees memory if needed.
  ~WavesWaveEqBoundaryRight() {
    if(freeIndexes)
      delete indexes;
  }


  bool applyWaveElliptic(real *x, real *y, 
			 const real &omega ) {
    
      omega_ = omega;
      
// cout<<"dt in  applyPlaneWaveOpexpl"<< step<<endl; 
      return doApply(x, y);
  }
  



 protected:
  ///
  bool doApply(const double *x, double *y) const ;
  

  /// Utility
  void doInitialize();

 private:
  double omega_;
  //
  bool freeIndexes;
};


//=======================================================================

class WavesWaveEqBoundaryBot : public WavesSDOperator 
{
public:
 
  /// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
  WavesWaveEqBoundaryBot(WavesSDIndexes  *sdi);

  /// Destructor frees memory if needed.
  ~WavesWaveEqBoundaryBot() {
    if(freeIndexes)
      delete indexes;
  }


  bool applyWaveElliptic(real *x, real *y, 
			 const real &omega ) {
    
      omega_ = omega;
      
// cout<<"dt in  applyPlaneWaveOpexpl"<< step<<endl; 
      return doApply(x, y);
  }
  



 protected:
  ///
  bool doApply(const double *x, double *y) const ;
  

  /// Utility
  void doInitialize();

 private:
  double omega_;
  //
  bool freeIndexes;
};



class WavesWaveEqBoundaryLeft : public WavesSDOperator 
{
public:
 
  /// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
  WavesWaveEqBoundaryLeft(WavesSDIndexes  *sdi);

  /// Destructor frees memory if needed.
  ~WavesWaveEqBoundaryLeft() {
    if(freeIndexes)
      delete indexes;
  }


  bool applyWaveElliptic(real *x, real *y, 
			 const real &omega ) {
    
      omega_ = omega;
      
// cout<<"dt in  applyPlaneWaveOpexpl"<< step<<endl; 
      return doApply(x, y);
  }
  



 protected:
  ///
  bool doApply(const double *x, double *y) const ;
  

  /// Utility
  void doInitialize();

 private:
  double omega_;
  //
  bool freeIndexes;
};


class WavesWaveEqBoundaryTop : public WavesSDOperator 
{
public:
 
  /// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
  WavesWaveEqBoundaryTop(WavesSDIndexes  *sdi);

  /// Destructor frees memory if needed.
  ~WavesWaveEqBoundaryTop() {
    if(freeIndexes)
      delete indexes;
  }


  bool applyWaveElliptic(real *x, real *y, 
			 const real &omega ) {
    
      omega_ = omega;
      
// cout<<"dt in  applyPlaneWaveOpexpl"<< step<<endl; 
      return doApply(x, y);
  }
  



 protected:
  ///
  bool doApply(const double *x, double *y) const ;
  

  /// Utility
  void doInitialize();

 private:
  double omega_;
  //
  bool freeIndexes;
};




#endif
  
