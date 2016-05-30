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
// Last changed: 2012-05-15 by Larisa Beilina
      

#include <math.h>

#include "include/wavesWaveEqOp123D.h"

//----------------------------------------------------------------------
//
// class WavesWaveEq123D in the case when velosity speed is function
//
//----------------------------------------------------------------------


WavesWaveEq123D::WavesWaveEq123D(WavesSDGeometry *sdg, double dt_, double  *coef_)
  : WavesSDOperator(new WavesSDInterior(*sdg)), dt(dt_), coef(coef_), freeIndexes(true)
{
  cout<<" present sdindexes in  WaveWaveEq123D"<<endl;

  //indexes->presentWavesSDIndexes(); 

  doInitialize();
}

WavesWaveEq123D::WavesWaveEq123D(WavesSDIndexes *sdi, double dt_, double *coef)
  : WavesSDOperator(sdi), dt(dt_), freeIndexes(false)
{
  doInitialize();
}

void WavesWaveEq123D::doInitialize() 
{
//  cout<<"inside doInitialize"<<endl;
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 

  // 1) create Laplace operator
  const int nsd = geometry->getNoSpaceDim();

   assert( nsd>=1 && nsd <= 3);

  nrOfWeights = (nsd==2) ? 5 : 7;

  if (nsd==1)
    nrOfWeights = 3;

  weights = new real[nrOfWeights];
  offsets = new int[nrOfWeights];

  real dx2 = geometry->getDx()*geometry->getDx();
  real dy2 = geometry->getDy()*geometry->getDy();
  real dz2 = geometry->getDz()*geometry->getDz();

  if(nsd == 3) {
    weights[0] = -2.0/dx2 - 2.0/dy2 -2.0/dz2 ; 
    weights[1] = weights[2] = 1.0/dx2;
    weights[3] = weights[4] = 1.0/dy2;
    weights[5] = weights[6] = 1.0/dz2;
    
    offsets[0] = 0;
    offsets[1] = geometry->node( -1, 0, 0 );
    offsets[2] = geometry->node(  1, 0, 0 );
    offsets[3] = geometry->node(  0,-1, 0 );
    offsets[4] = geometry->node(  0, 1, 0 );
    offsets[5] = geometry->node(  0, 0,-1 );
    offsets[6] = geometry->node(  0, 0, 1 );
  }
  else if (nsd==2){
    weights[0] = -2.0/dx2 - 2.0/dy2; 
    weights[1] = weights[2] = 1.0/dx2;
    weights[3] = weights[4] = 1.0/dy2;

    offsets[0] = 0;
    offsets[1] = geometry->node( -1, 0, 0 );
    offsets[2] = geometry->node(  1, 0, 0 );
    offsets[3] = geometry->node(  0,-1, 0 );
    offsets[4] = geometry->node(  0, 1, 0 );
  }  
  else if (nsd==1)
    {
      weights[0] = -2.0/dx2; 
      weights[1] = weights[2] = 1.0/dx2;
      
      offsets[0] = 0;
      offsets[1] = geometry->node( -1, 0, 0 );
      offsets[2] = geometry->node(  1, 0, 0 );
    
    }  

  // 2) y := dt*dt Laplace + 2
  for ( int i = 0; i<nrOfWeights; i++ )
    weights[i] *= dt*dt;
  
   
  //  weights[0] += 2;
}


bool WavesWaveEq123D::doApply(const double *x, double *y) const
{
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 
  cout<<"nsd = "<< geometry->getNoSpaceDim()<<endl;
  indexes->presentWavesSDIndexes();
  const int nsd = geometry->getNoSpaceDim();
  assert( nsd>=1 && nsd <= 3);

  // We apply y <- (dt^2 Laplace + 2) x - y

  // For faster loops:
  real weights0 = weights[0];
  real weights1 = weights[1];
  real weights2 = weights[2];

  real weights3 = weights[3]; 
  real weights4 = weights[4];


  real weights5 = (nsd == 3)? weights[5] : 0.;
  real weights6 = (nsd == 3)? weights[6] : 0.;

  int offsets1 = offsets[1];
  int offsets2 = offsets[2];
  int offsets3 = offsets[3];
  int offsets4 = offsets[4];
  int offsets5 = (nsd == 3)? offsets[5] : 0;
  int offsets6 = (nsd == 3)? offsets[6] : 0;

  if (nsd==1)
    {
      weights3 = 0.0;
      weights4 = 0.0;
      weights5 = 0.0;
      weights6 = 0.0;
      offsets3 = 0;
      offsets4 = 0;
      offsets5 = 0;
      offsets6 = 0;
    }

  for(int loop = 0; loop< indexes->nOfLoopIndex(); loop++ )
    {
      int lstart = indexes->loopStart( loop );
      int lstop  = indexes->loopStop( loop );

      //  cout << "### WaveEqInterior " << loop << " from " << lstart << " to "  << lstop << endl;
 

      if(nsd == 3)
	for(int n = lstart; n<= lstop; n++ )
	  {
	    //    cout<<"uOuter("<<n<<") = "<<y[n]<<endl;
	  y[n] = -y[n] +
	    weights0 * x[ n ] +
	    weights1 * x[ n+offsets1 ] +
	    weights2 * x[ n+offsets2 ] +
	    weights3 * x[ n+offsets3 ] +
	    weights4 * x[ n+offsets4 ] +
	    weights5 * x[ n+offsets5 ] +
	    weights6 * x[ n+offsets6 ];
	  //  cout<<"uOuter("<<n<<") = "<<y[n]<<endl;
	  //cout<<"x("<<n<<") = "<<x[n]<<endl;
	  //cout<<"x("<<n+offsets1<<") = "<<x[n+offsets1]<<endl;
	  //cout<<"x("<<n+offsets2<<") = "<<x[n+offsets2]<<endl;
	  //cout<<"x("<<n+offsets3<<") = "<<x[n+offsets3]<<endl;
	  // cout<<"x("<<n+offsets4<<") = "<<x[n+offsets4]<<endl;
	  //cout<<"x("<<n+offsets5<<") = "<<x[n+offsets5]<<endl;
	  //cout<<"x("<<n+offsets6<<") = "<<x[n+offsets6]<<endl;
	  //  cout<<"offsets1 = "<<offsets1<<"offsets2 = "<<offsets2<<"offsets3         // = "<<offsets3<<"offsets4 = "<<offsets4<<" offsets5 ="<<offsets5<<"          //offsets6 = "<<offsets6<<endl;                                             
	  }	  
      else if (nsd ==2)
	for(int n = lstart; n<= lstop; n++ )
	  {
	  y[n] =  -y[n] +
	    weights0 * x[ n ] +
	    weights1 * x[ n+offsets1 ] +
	    weights2 * x[ n+offsets2 ] +
	    weights3 * x[ n+offsets3 ] +
	    weights4 * x[ n+offsets4 ];	
      cout<<" u_fdm("<<n<<")="<<y[n]<<endl;	  
	  }    
      else if (nsd == 1)
	for(int n = lstart; n<= lstop; n++ )
	  {
	    y[n] =  -y[n] + (weights0 * x[ n ] +  weights1 * x[ n+offsets1 ] +
	       weights2 * x[ n+offsets2 ])/coef[n] + 2*x[n];

	    cout<<" coefficient ("<<n<<")="<<coef[n]<<endl;	  
	  }    
      
}
 
  return true;
}

