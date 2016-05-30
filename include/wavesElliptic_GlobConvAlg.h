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


#ifndef _WAVESELLIPTICGLOBCONVALG_
#define _WAVESELLIPTICGLOBCONVALG_

//#include "wavesSDOperator.h"


#include "wavesSDDefs.h"
#include "waveskraftwerk.hh"
#include "wavesfindBoundaryNodes.h"
#include "wavesutil2.h"
#include "wavesFuncDefs.h"
#include "wavesBCOperators.h"
#include "wavesOptional.h"
#include "wavesOutputs.h"

class WavesElliptic_GlobConvAlgI3 : public WavesSDOperator 
{
public:
  /// SDindexes covering the interior is created. The time step dt shall be supplied.
  WavesElliptic_GlobConvAlgI3(WavesSDGeometry *sdg); 

  /// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
  WavesElliptic_GlobConvAlgI3(WavesSDIndexes  *sdi);

  /// Destructor frees memory if needed.
  ~WavesElliptic_GlobConvAlgI3() {
    if(freeIndexes)
      delete indexes;
  }


  bool applyWaveElliptic(real *x, real *y,  real* grad_q_i,   real* grad_q_square ) {
   
      grad_q_i_ =  grad_q_i;
      grad_q_square_ =  grad_q_square;
// cout<<"dt in  applyPlaneWaveOpexpl"<< step<<endl; 
      return doApply(x, y);
  }
  

int delta(const int& s, const int& n) const
    {
      if (s==n)
	return 1;
      else
	return 0; 
      
    }

 protected:
  ///
  bool doApply(const double *x, double *y) const ;
  

  /// Utility
  void doInitialize();

 private:
  real* grad_q_i_;
  real* grad_q_square_;
  //
  bool freeIndexes;
};

// ===========================================================================

class WavesElliptic_GlobConvAlgI2 : public WavesSDOperator 
{
public:
  /// SDindexes covering the interior is created. The time step dt shall be supplied.
  WavesElliptic_GlobConvAlgI2(WavesSDGeometry *sdg); 

  /// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
  WavesElliptic_GlobConvAlgI2(WavesSDIndexes  *sdi);

  /// Destructor frees memory if needed.
  ~WavesElliptic_GlobConvAlgI2() {
    if(freeIndexes)
      delete indexes;
  }


  bool applyWaveElliptic(real *x, real *y,  real* grad_q_i,  real* grad_q_i_CWF ) {
   
      grad_q_i_ =  grad_q_i;
      grad_q_i_CWF_ =  grad_q_i_CWF;

      
// cout<<"dt in  applyPlaneWaveOpexpl"<< step<<endl; 
      return doApply(x, y);
  }
  

  int delta(int& s, int n) const
    {
      if (s==n)
	return 1;
      else
	return 0; 
      
    }

 protected:
  ///
  bool doApply(const double *x, double *y) const ;
  

  /// Utility
  void doInitialize();

 private:
  real* grad_q_i_;
  real* grad_q_i_CWF_;
  //
  bool freeIndexes;
};

//=================================================================================0

class WavesElliptic_GlobConvAlgI1 : public WavesSDOperator 
{
public:
  /// SDindexes covering the interior is created. The time step dt shall be supplied.
  WavesElliptic_GlobConvAlgI1(WavesSDGeometry *sdg); 

  /// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
  WavesElliptic_GlobConvAlgI1(WavesSDIndexes  *sdi);

  /// Destructor frees memory if needed.
  ~WavesElliptic_GlobConvAlgI1() {
    if(freeIndexes)
      delete indexes;
  }


  bool applyWaveElliptic(real *x, real *y,  real* grad_q_i ) {
   
      grad_q_i_ =  grad_q_i;
      
// cout<<"dt in  applyPlaneWaveOpexpl"<< step<<endl; 
      return doApply(x, y);
  }
  

  int delta(int& s, int n) const
    {
      if (s==n)
	return 1;
      else
	return 0; 
      
    }

 protected:
  ///
  bool doApply(const double *x, double *y) const ;
  

  /// Utility
  void doInitialize();

 private:
  real* grad_q_i_;
  //
  bool freeIndexes;
};

//==========================================================================
class WavesElliptic_GlobConvAlg_CWF : public WavesSDOperator 
{
public:
  /// SDindexes covering the interior is created. The time step dt shall be supplied.
  WavesElliptic_GlobConvAlg_CWF(WavesSDGeometry *sdg); 
  WavesElliptic_GlobConvAlg_CWF(WavesSDGeometry *sdg, int code); 

  /// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
  WavesElliptic_GlobConvAlg_CWF(WavesSDIndexes  *sdi);

  /// Destructor frees memory if needed.
  ~WavesElliptic_GlobConvAlg_CWF() {
    if(freeIndexes)
      delete indexes;
  }


  int delta(int& s, int n)  const
    {
      if (s==n)
	return 1;
      else
	return 0; 
      
    }

  bool applyWaveElliptic(real *x, real *y) {
    

    //      func_q_i_ = func_q_i;

          
// cout<<"dt in  applyPlaneWaveOpexpl"<< step<<endl; 
     return doApply(x, y);
    
  }

  void doInitGuess(const double *x,  double *y);
 
 protected:
  ///
  bool doApply(const double *x, double *y) const ;
  

  /// Utility
  void doInitialize();

 private:
  real *func_q_i_;
  //
  bool freeIndexes;
};


//==========================================================================
class WavesElliptic_GlobConvAlg : public WavesSDOperator 
{
public:
  /// SDindexes covering the interior is created. The time step dt shall be supplied.
  WavesElliptic_GlobConvAlg(WavesSDGeometry *sdg); 
  WavesElliptic_GlobConvAlg(WavesSDGeometry *sdg, int code); 

  /// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
  WavesElliptic_GlobConvAlg(WavesSDIndexes  *sdi);

  /// Destructor frees memory if needed.
  ~WavesElliptic_GlobConvAlg() {
    if(freeIndexes)
      delete indexes;
  }


  int delta(int& s, int n) const
    {
      if (s==n)
	return 1;
      else
	return 0; 
      
    }

  bool applyWaveElliptic(real *x, real *y) {
    

    //      func_q_i_ = func_q_i;

          
// cout<<"dt in  applyPlaneWaveOpexpl"<< step<<endl; 
     return doApply(x, y);
    
  }

  void doInitGuess(const double *x,  double *y);
 
 protected:
  /// 
  bool doApply(const double *x, double *y) const ;
  

  /// Utility
  void doInitialize();

 private:
  real *func_q_i_;
  //
  bool freeIndexes;
};

//********************************************************************************

class WavesElliptic_gradient : public WavesSDOperator 
{
public:
  /// SDindexes covering the interior is created. The time step dt shall be supplied.
  WavesElliptic_gradient(WavesSDGeometry *sdg); 
  WavesElliptic_gradient(WavesSDGeometry *sdg, int code); 

  /// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
  WavesElliptic_gradient(WavesSDIndexes  *sdi);

  /// Destructor frees memory if needed.
  ~WavesElliptic_gradient() {
    if(freeIndexes)
      delete indexes;
  }


  int delta(int& s, int n) const
    {
      if (s==n)
	return 1;
      else
	return 0; 
      
    }

  bool applyWaveElliptic(real *x, real *y) {
    

    //      func_q_i_ = func_q_i;

          
// cout<<"dt in  applyPlaneWaveOpexpl"<< step<<endl; 
     return doApply(x, y);
    
  }

  void compute_gradient(const double *x,  double *y);

  void compute_gradient(const double *x,
			double *y_1,
			double *y_2,
			double *y_3);

  void compute_gradient(MV_Vector<double> x, MV_Vector<double>& y_1, MV_Vector<double>& y_2, MV_Vector<double>& y_3); //as the above function, just different data type, added by Thanh

  void compute_square_gradient(const double *x,  double *y);
 
 protected:
  ///
  bool doApply(const double *x, double *y) const ;
  

  /// Utility
  void doInitialize();

 private:
  // real *func_q_i_;
  //
  bool freeIndexes;
};


//===============================

#endif
  
