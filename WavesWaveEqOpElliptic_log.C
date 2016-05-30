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

#include <math.h>
#include "include/wavesWaveEqOpElliptic_log.h"

//----------------------------------------------------------------------
//
// class WaveEqInteriorElliptic
//
//----------------------------------------------------------------------


WavesWaveEqInteriorElliptic::WavesWaveEqInteriorElliptic(WavesSDGeometry *sdg)
    : WavesSDOperator( new WavesSDInterior(*sdg) ), freeIndexes(true)
{
  doInitialize();


}

WavesWaveEqInteriorElliptic::WavesWaveEqInteriorElliptic(WavesSDIndexes *sdi)
  : WavesSDOperator( sdi ), freeIndexes(false)
{
  doInitialize();
}

void WavesWaveEqInteriorElliptic::doInitialize() 
{
 // cout<<"inside doInitialize"<<endl;
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 

  // 1) create Laplace operator
  const int nsd = geometry->getNoSpaceDim();

   assert( nsd>=2 && nsd <= 3);

  nrOfWeights = (nsd==2) ? 7 : 10;
  
  weights = new real[nrOfWeights];
  offsets = new int[nrOfWeights];

  real dx2 = geometry->getDx()*geometry->getDx();
  real dy2 = geometry->getDy()*geometry->getDy();
  real dz2 = geometry->getDz()*geometry->getDz();

  real dx = geometry->getDx();
  real dy = geometry->getDy();
  real dz = geometry->getDz();

  if(nsd == 3) {
    weights[0] = -2.0/dx2 - 2.0/dy2 -2.0/dz2 ; 
    weights[1] = weights[2] = 1.0/dx2;
    weights[3] = weights[4] = 1.0/dy2;
    weights[5] = weights[6] = 1.0/dz2;
    weights[7] = 1.0/dx;
    weights[8] = 1.0/dy;
    weights[9] = 1.0/dz;

    offsets[0] = 0;
    offsets[1] = geometry->node( -1, 0, 0 );
    offsets[2] = geometry->node(  1, 0, 0 );
    offsets[3] = geometry->node(  0,-1, 0 );
    offsets[4] = geometry->node(  0, 1, 0 );
    offsets[5] = geometry->node(  0, 0,-1 );
    offsets[6] = geometry->node(  0, 0, 1 );
  }
  else {
    weights[0] = -2.0/dx2 - 2.0/dy2; 
    weights[1] = weights[2] = 1.0/dx2;
    weights[3] = weights[4] = 1.0/dy2;
    weights[5] = 1.0/dx;
    weights[6] = 1.0/dy;


    offsets[0] = 0;
    offsets[1] = geometry->node( -1, 0, 0 );
    offsets[2] = geometry->node(  1, 0, 0 );
    offsets[3] = geometry->node(  0,-1, 0 );
    offsets[4] = geometry->node(  0, 1, 0 );

    cout<<" dx = "<<dx<<"  and dy = "<<dy<<endl;

  }  

}


bool WavesWaveEqInteriorElliptic::doApply(const double *x,  double *y) const
{
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 
  //cout<<"nsd = "<< geometry->getNoSpaceDim()<<endl;
  // indexes->presentWavesSDIndexes();
  const int nsd = geometry->getNoSpaceDim();
  assert( nsd>=2 && nsd <= 3);

  // We apply y <- (dt^2 Laplace + 2) x - y

  // For faster loops:
  real weights0 = weights[0];
  real weights1 = weights[1];
  real weights2 = weights[2];
  real weights3 = weights[3];
  real weights4 = weights[4];
  real weights5 = weights[5];
  real weights6 = weights[6];
  real weights7 = (nsd == 3)? weights[7] : 0.; 
  real weights8 = (nsd == 3)? weights[8] : 0.;
  real weights9 = (nsd == 3)? weights[9] : 0.;

  int offsets1 = offsets[1];
  int offsets2 = offsets[2];
  int offsets3 = offsets[3];
  int offsets4 = offsets[4];
  int offsets5 = (nsd == 3)? offsets[5] : 0;
  int offsets6 = (nsd == 3)? offsets[6] : 0;

  for(int loop = 0; loop< indexes->nOfLoopIndex(); loop++ )
    {
      int lstart = indexes->loopStart( loop );
      int lstop  = indexes->loopStop( loop );

       cout << "### WaveEqInterior " << loop << " from " << lstart << " to " 
        << lstop << endl;
 

      if(nsd == 3)
      {
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
      }	  
      else if (nsd==2)
      {
	  for(int n = lstart; n<= lstop; n++ )
	  {


// elliptic equation obtained after applying Laplace transform

// equation (1.3) in my remarks or: laplace w - omega^2 c(x) omega = -delta(x-x_0)
	  
	    // elliptic equation after transformation w=ln v
	    
 // y[n] =  
 // (weights0 * x[ n ] +
 // weights1 * x[ n+offsets1 ] +
//  weights2 * x[ n+offsets2 ] +
//  weights3 * x[ n+offsets3 ] +
//  weights4 * x[ n+offsets4 ] +
 //  (weights5*(x[ n+offsets2 ] - x[n]))*(weights5*(x[ n+offsets2 ] - x[n]))
  // +   (weights6*(x[n+offsets4 ] - x[n]))*(weights6*(x[ n+offsets4 ] - x[n])))/(omega_*omega_);
	    


  y[n] =  
  weights0 * x[ n ] +
  weights1 * x[ n+offsets1 ] +
  weights2 * x[ n+offsets2 ] +
  weights3 * x[ n+offsets3 ] +
  weights4 * x[ n+offsets4 ] +
   (omega_*omega_)*(weights5*(x[ n+offsets2 ] - x[n])*weights5*(x[ n+offsets2 ] - x[n])
   +   weights6*(x[n+offsets4 ] - x[n])*weights6*(x[ n+offsets4 ] - x[n]));

	    cout<<" in computing operator ( "<<n<<")="<<y[n]<<endl;

	   
// cout<<"x("<<n<<") = "<<x[n]<<endl;
	//  cout<<"x("<<n+offsets1<<") = "<<x[n+offsets1]<<endl;
	//  cout<<"x("<<n+offsets2<<") = "<<x[n+offsets2]<<endl;
	//  cout<<"x("<<n+offsets3<<") = "<<x[n+offsets3]<<endl;
	//   cout<<"x("<<n+offsets4<<") = "<<x[n+offsets4]<<endl;
//	   cout<<" laplace = "<<weights0 * x[ n ] +weights1 * x[ n+offsets1 ] +
// weights2 * x[ n+offsets2 ] +weights3 * x[ n+offsets3 ] +weights4 * x[ n+offsets4 ]<<endl;
         
//	   cout<<" d/dx "<<weights5*(x[ n+offsets2 ] - x[n])*weights5*(x[ n+offsets2 ] - x[n])<<endl;
//		   cout<<" d/dy "<<weights6*(x[ n+offsets4] - x[n])*weights6*(x[ n+offsets4 ] - x[n])<<endl;

	   
	  }    
      }
    }
 
    return true;
    }


//----------------------------------------------------------------------
//
// class WavesWaveEqInteriorLaplace
//
//----------------------------------------------------------------------


WavesWaveEqInteriorLaplace::WavesWaveEqInteriorLaplace(WavesSDGeometry *sdg)
    : WavesSDOperator( new WavesSDInterior(*sdg) ), freeIndexes(true)
{
  doInitialize();


}

WavesWaveEqInteriorLaplace::WavesWaveEqInteriorLaplace(WavesSDIndexes *sdi)
  : WavesSDOperator( sdi ), freeIndexes(false)
{
  doInitialize();
}

void WavesWaveEqInteriorLaplace::doInitialize() 
{
//  cout<<"inside doInitialize"<<endl;
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 

  // 1) create Laplace operator
  const int nsd = geometry->getNoSpaceDim();

   assert( nsd>=2 && nsd <= 3);

  nrOfWeights = (nsd==2) ? 7 : 10;
  
  weights = new real[nrOfWeights];
  offsets = new int[nrOfWeights];

  real dx2 = geometry->getDx()*geometry->getDx();
  real dy2 = geometry->getDy()*geometry->getDy();
  real dz2 = geometry->getDz()*geometry->getDz();

  real dx = geometry->getDx();
  real dy = geometry->getDy();
  real dz = geometry->getDz();

  if(nsd == 3) {
    weights[0] = -2.0/dx2 - 2.0/dy2 -2.0/dz2 ; 
    weights[1] = weights[2] = 1.0/dx2;
    weights[3] = weights[4] = 1.0/dy2;
    weights[5] = weights[6] = 1.0/dz2;
    weights[7] = 1.0/dx;
    weights[8] = 1.0/dy;
    weights[9] = 1.0/dz;

    offsets[0] = 0;
    offsets[1] = geometry->node( -1, 0, 0 );
    offsets[2] = geometry->node(  1, 0, 0 );
    offsets[3] = geometry->node(  0,-1, 0 );
    offsets[4] = geometry->node(  0, 1, 0 );
    offsets[5] = geometry->node(  0, 0,-1 );
    offsets[6] = geometry->node(  0, 0, 1 );
  }
  else {
    weights[0] = -2.0/dx2 - 2.0/dy2; 
    weights[1] = weights[2] = 1.0/dx2;
    weights[3] = weights[4] = 1.0/dy2;
    weights[5] = 1.0/dx;
    weights[6] = 1.0/dy;


    offsets[0] = 0;
    offsets[1] = geometry->node( -1, 0, 0 );
    offsets[2] = geometry->node(  1, 0, 0 );
    offsets[3] = geometry->node(  0,-1, 0 );
    offsets[4] = geometry->node(  0, 1, 0 );

    cout<<" dx = "<<dx<<"  and dy = "<<dy<<endl;

  }  

}


bool WavesWaveEqInteriorLaplace::doApply(const double *x,  double *y) const
{
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 
  //cout<<"nsd = "<< geometry->getNoSpaceDim()<<endl;
  // indexes->presentSDIndexes();
  const int nsd = geometry->getNoSpaceDim();
  assert( nsd>=2 && nsd <= 3);

  // We apply y <- (dt^2 Laplace + 2) x - y

  // For faster loops:
  real weights0 = weights[0];
  real weights1 = weights[1];
  real weights2 = weights[2];
  real weights3 = weights[3];
  real weights4 = weights[4];
  real weights5 = weights[5];
  real weights6 = weights[6];
  real weights7 = (nsd == 3)? weights[7] : 0.; 
  real weights8 = (nsd == 3)? weights[8] : 0.;
  real weights9 = (nsd == 3)? weights[9] : 0.;

  int offsets1 = offsets[1];
  int offsets2 = offsets[2];
  int offsets3 = offsets[3];
  int offsets4 = offsets[4];
  int offsets5 = (nsd == 3)? offsets[5] : 0;
  int offsets6 = (nsd == 3)? offsets[6] : 0;

  for(int loop = 0; loop< indexes->nOfLoopIndex(); loop++ )
    {
      int lstart = indexes->loopStart( loop );
      int lstop  = indexes->loopStop( loop );

       cout << "### WaveEqInterior " << loop << " from " << lstart << " to " 
        << lstop << endl;
 

      if(nsd == 3)
      {
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
      }	  
      else if (nsd==2)
      {
	  for(int n = lstart; n<= lstop; n++ )
	  {


// elliptic equation obtained after applying Laplace transform

// equation (1.3) in my remarks or: laplace w - omega^2 c(x) omega = -delta(x-x_0)
	    
	    //     y[n] =  (weights0 * x[ n ] +
	    //	       weights1 * x[ n+offsets1 ] +
	    //		       weights2 * x[ n+offsets2 ] +
	    //       weights3 * x[ n+offsets3 ] +
	    //       weights4 * x[ n+offsets4 ])/(omega_*omega_*x[n]);
	    



	    // elliptic equation after transformation w=ln v
	    
	    // y[n] =  
	    // (weights0 * x[ n ] +
	    //  weights1 * x[ n+offsets1 ] +
	    // weights2 * x[ n+offsets2 ] +
	    //  weights3 * x[ n+offsets3 ] +
	    //  weights4 * x[ n+offsets4 ] +
	    //  (weights5*(x[ n+offsets2 ] - x[n]))*(weights5*(x[ n+offsets2 ] - x[n]))
	    // +   (weights6*(x[n+offsets4 ] - x[n]))*(weights6*(x[ n+offsets4 ] - x[n])))/(omega_*omega_);
	    


  y[n] =  
  weights0 * x[ n ] +
  weights1 * x[ n+offsets1 ] +
  weights2 * x[ n+offsets2 ] +
  weights3 * x[ n+offsets3 ] +
  weights4 * x[ n+offsets4 ];

	    cout<<" in computing operator ( "<<n<<")="<<y[n]<<endl;

	    
	    // cout<<"x("<<n<<") = "<<x[n]<<endl;
	    //  cout<<"x("<<n+offsets1<<") = "<<x[n+offsets1]<<endl;
	    // cout<<"x("<<n+offsets2<<") = "<<x[n+offsets2]<<endl;
	    // cout<<"x("<<n+offsets3<<") = "<<x[n+offsets3]<<endl;
	    // cout<<"x("<<n+offsets4<<") = "<<x[n+offsets4]<<endl;
	    // cout<<" laplace = "<<weights0 * x[ n ] +weights1 * x[ n+offsets1 ] + weights2 * x[ n+offsets2 ] +weights3 * x[ n+offsets3 ] +weights4 * x[ n+offsets4 ]<<endl;
	   
	    // cout<<" d/dx "<<weights5*(x[ n+offsets2 ] - x[n])*weights5*(x[ n+offsets2 ] - x[n])<<endl;

	    //   cout<<" d/dy "<<weights6*(x[ n+offsets4] - x[n])*weights6*(x[ n+offsets4 ] - x[n])<<endl;

	   
	  }    
      }
    }
 
    return true;
    }



//----------------------------------------------------------------------
//
// class WavesWaveEqInteriorElliptic_laplace
//
//----------------------------------------------------------------------


WavesWaveEqInteriorElliptic_laplace::WavesWaveEqInteriorElliptic_laplace(WavesSDGeometry *sdg)
    : WavesSDOperator( new WavesSDInterior(*sdg) ), freeIndexes(true)
{
  doInitialize();
}

WavesWaveEqInteriorElliptic_laplace::WavesWaveEqInteriorElliptic_laplace(WavesSDIndexes *sdi)
  : WavesSDOperator( sdi ), freeIndexes(false)
{
  doInitialize();
}

void WavesWaveEqInteriorElliptic_laplace::doInitialize() 
{
//  cout<<"inside doInitialize"<<endl;
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 

  // 1) create Laplace operator
  const int nsd = geometry->getNoSpaceDim();

   assert( nsd>=2 && nsd <= 3);

  nrOfWeights = (nsd==2) ? 7 : 10;
  
  weights = new real[nrOfWeights];
  offsets = new int[nrOfWeights];

  real dx2 = geometry->getDx()*geometry->getDx();
  real dy2 = geometry->getDy()*geometry->getDy();
  real dz2 = geometry->getDz()*geometry->getDz();

  real dx = geometry->getDx();
  real dy = geometry->getDy();
  real dz = geometry->getDz();

  if(nsd == 3) {
    weights[0] = -2.0/dx2 - 2.0/dy2 -2.0/dz2 ; 
    weights[1] = weights[2] = 1.0/dx2;
    weights[3] = weights[4] = 1.0/dy2;
    weights[5] = weights[6] = 1.0/dz2;
    weights[7] = 1.0/dx;
    weights[8] = 1.0/dy;
    weights[9] = 1.0/dz;

    offsets[0] = 0;
    offsets[1] = geometry->node( -1, 0, 0 );
    offsets[2] = geometry->node(  1, 0, 0 );
    offsets[3] = geometry->node(  0,-1, 0 );
    offsets[4] = geometry->node(  0, 1, 0 );
    offsets[5] = geometry->node(  0, 0,-1 );
    offsets[6] = geometry->node(  0, 0, 1 );
  }
  else {
    weights[0] = -2.0/dx2 - 2.0/dy2; 
    weights[1] = weights[2] = 1.0/dx2;
    weights[3] = weights[4] = 1.0/dy2;
    weights[5] = 1.0/dx;
    weights[6] = 1.0/dy;


    offsets[0] = 0;
    offsets[1] = geometry->node( -1, 0, 0 );
    offsets[2] = geometry->node(  1, 0, 0 );
    offsets[3] = geometry->node(  0,-1, 0 );
    offsets[4] = geometry->node(  0, 1, 0 );
  }  

}


bool WavesWaveEqInteriorElliptic_laplace::doApply(const double *x,  double *y) const
{
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 
  //cout<<"nsd = "<< geometry->getNoSpaceDim()<<endl;
  // indexes->presentSDIndexes();
  const int nsd = geometry->getNoSpaceDim();
  assert( nsd>=2 && nsd <= 3);

  // We apply y <- (dt^2 Laplace + 2) x - y

  // For faster loops:
  real weights0 = weights[0];
  real weights1 = weights[1];
  real weights2 = weights[2];
  real weights3 = weights[3];
  real weights4 = weights[4];
  real weights5 = weights[5];
  real weights6 = weights[6];
  real weights7 = (nsd == 3)? weights[7] : 0.; 
  real weights8 = (nsd == 3)? weights[8] : 0.;
  real weights9 = (nsd == 3)? weights[9] : 0.;

  int offsets1 = offsets[1];
  int offsets2 = offsets[2];
  int offsets3 = offsets[3];
  int offsets4 = offsets[4];
  int offsets5 = (nsd == 3)? offsets[5] : 0;
  int offsets6 = (nsd == 3)? offsets[6] : 0;

  for(int loop = 0; loop< indexes->nOfLoopIndex(); loop++ )
    {
      int lstart = indexes->loopStart( loop );
      int lstop  = indexes->loopStop( loop );

   //    cout << "### WaveEqInterior " << loop << " from " << lstart << " to " 
      //  << lstop << endl;
 

      if(nsd == 3)
      {
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
      }	  
      else if (nsd==2)
      {
	  for(int n = lstart; n<= lstop; n++ )
	  {


// elliptic equation obtained after applying Laplace transform

// equation (1.3) in my remarks or: laplace w - omega^2 c(x) omega = -delta(x-x_0)
	    
	    //    y[n] =  (weights0 * x[ n ] +
	    //	       weights1 * x[ n+offsets1 ] +
	    //	       weights2 * x[ n+offsets2 ] +
	    //	       weights3 * x[ n+offsets3 ] +
	    //	       weights4 * x[ n+offsets4 ])/(omega_*omega_*x[n]);
	    
	    // elliptic equation after transformation w=ln v
	    
	    //y[n] =  
	    //(weights0 * x[ n ] +
	    //weights1 * x[ n+offsets1 ] +
	    //weights2 * x[ n+offsets2 ] +
	    //weights3 * x[ n+offsets3 ] +
	    //weights4 * x[ n+offsets4 ] +
	    //(weights5*(x[ n+offsets2 ] - x[n]))*(weights5*(x[ n+offsets2 ] - x[n]))
	    //+   (weights6*(x[n+offsets4 ] - x[n]))*(weights6*(x[ n+offsets4 ] - x[n])))/(omega_*omega_);
	    

	    // first compute parameter c
  y[n] =  
  weights0 * x[ n ] +
  weights1 * x[ n+offsets1 ] +
  weights2 * x[ n+offsets2 ] +
  weights3 * x[ n+offsets3 ] +
  weights4 * x[ n+offsets4 ] +
   (omega_*omega_)*(weights5*(x[ n+offsets2 ] - x[n])*weights5*(x[ n+offsets2 ] - x[n])
   +   weights6*(x[n+offsets4 ] - x[n])*weights6*(x[ n+offsets4 ] - x[n]));

  // then again compute laplacian, second derivatives, only to check convergence for omega

  y[n] -=  (omega_*omega_)*(weights5*(x[ n+offsets2 ] - x[n])*weights5*(x[ n+offsets2 ] - x[n])
		  +   weights6*(x[n+offsets4 ] - x[n])*weights6*(x[ n+offsets4 ] - x[n]));

	    cout<<" for log  y( "<<n<<")="<<y[n]<<endl;
	  }    
      }
    }
 
    return true;
    }

//=============================================================
 
//----------------------------------------------------------------------
//
// class WavesWaveEqBoundaryRight
//
//----------------------------------------------------------------------


WavesWaveEqBoundaryRight::WavesWaveEqBoundaryRight(WavesSDIndexes *sdi)
  : WavesSDOperator( sdi ), freeIndexes(false)
{
  doInitialize();
}

void WavesWaveEqBoundaryRight::doInitialize() 
{
//  cout<<"inside doInitialize"<<endl;
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 

  // 1) create Laplace operator
  const int nsd = geometry->getNoSpaceDim();

   assert( nsd>=2 && nsd <= 3);

  nrOfWeights = (nsd==2) ? 7 : 10;
  
  weights = new real[nrOfWeights];
  offsets = new int[nrOfWeights];

  real dx2 = geometry->getDx()*geometry->getDx();
  real dy2 = geometry->getDy()*geometry->getDy();
  real dz2 = geometry->getDz()*geometry->getDz();

  real dx = geometry->getDx();
  real dy = geometry->getDy();
  real dz = geometry->getDz();

  if(nsd == 3) {
    weights[0] = -2.0/dx2 - 2.0/dy2 -2.0/dz2 ; 
    weights[1] = weights[2] = 1.0/dx2;
    weights[3] = weights[4] = 1.0/dy2;
    weights[5] = weights[6] = 1.0/dz2;
    weights[7] = 1.0/dx;
    weights[8] = 1.0/dy;
    weights[9] = 1.0/dz;

    offsets[0] = 0;
    offsets[1] = geometry->node( -1, 0, 0 );
    offsets[2] = geometry->node(  1, 0, 0 );
    offsets[3] = geometry->node(  0,-1, 0 );
    offsets[4] = geometry->node(  0, 1, 0 );
    offsets[5] = geometry->node(  0, 0,-1 );
    offsets[6] = geometry->node(  0, 0, 1 );
  }
  else {
    weights[0] = 1.0/dx2 - 2.0/dy2; 
    weights[1] = weights[2] = 1.0/dx2;
    weights[3] = weights[4] = 1.0/dy2;
    weights[5] = 1.0/dx;
    weights[6] = 1.0/dy;


    offsets[0] = 0;
    offsets[1] = geometry->node( -1, 0, 0 );
    offsets[2] = geometry->node(  -2, 0, 0 );
    offsets[3] = geometry->node(  0,-1, 0 );
    offsets[4] = geometry->node(  0, 1, 0 );
   

    cout<<" dx = "<<dx<<"  and dy = "<<dy<<endl;

  }  

}


bool WavesWaveEqBoundaryRight::doApply(const double *x,  double *y) const
{
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 
  //cout<<"nsd = "<< geometry->getNoSpaceDim()<<endl;
  // indexes->presentSDIndexes();
  const int nsd = geometry->getNoSpaceDim();
  assert( nsd>=2 && nsd <= 3);

  // We apply y <- (dt^2 Laplace + 2) x - y

  // For faster loops:
  real weights0 = weights[0];
  real weights1 = weights[1];
  real weights2 = weights[2];
  real weights3 = weights[3];
  real weights4 = weights[4];
  real weights5 = weights[5];
  real weights6 = weights[6];
  real weights7 = (nsd == 3)? weights[7] : 0.; 
  real weights8 = (nsd == 3)? weights[8] : 0.;
  real weights9 = (nsd == 3)? weights[9] : 0.;

  int offsets1 = offsets[1];
  int offsets2 = offsets[2];
  int offsets3 = offsets[3];
  int offsets4 = offsets[4];
  int offsets5 = (nsd == 3)? offsets[5] : 0;
  int offsets6 = (nsd == 3)? offsets[6] : 0;

  for(int loop = 0; loop< indexes->nOfLoopIndex(); loop++ )
    {
      int lstart = indexes->loopStart( loop );
      int lstop  = indexes->loopStop( loop );

       cout << "### WavesWaveEqBoundary " << loop << " from " << lstart << " to " 
        << lstop << endl;
 

      if(nsd == 3)
      {
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
      }	  
      else if (nsd==2)
      {
	  for(int n = lstart; n<= lstop; n++ )
	  {


// elliptic equation obtained after applying Laplace transform

// equation (1.3) in my remarks or: laplace w - omega^2 c(x) omega = -delta(x-x_0)
	   

	    // elliptic equation after transformation w=ln v



  y[n] =  
  weights0 * x[ n ] -
  2*weights1 * x[ n+offsets1 ] +
   weights2 * x[ n+offsets2 ] +
  weights3 * x[ n+offsets3 ] +
  weights4 * x[ n+offsets4 ] +
   (omega_*omega_)*(weights5*(x[n] - x[ n+offsets1])*weights5*(x[n] - x[ n+offsets1])
   +   weights6*(x[n+offsets4 ] - x[n])*weights6*(x[ n+offsets4 ] - x[n]));

	    cout<<" in computing operator ( "<<n<<")="<<y[n]<<endl;

	  
	  }    
      }
    }
 
    return true;
    }

//====================================================================================================
//----------------------------------------------------------------------
//
// class WavesWaveEqBoundaryBot
//
//----------------------------------------------------------------------


WavesWaveEqBoundaryBot::WavesWaveEqBoundaryBot(WavesSDIndexes *sdi)
  : WavesSDOperator( sdi ), freeIndexes(false)
{
  doInitialize();
}

void WavesWaveEqBoundaryBot::doInitialize() 
{
//  cout<<"inside doInitialize"<<endl;
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 

  // 1) create Laplace operator
  const int nsd = geometry->getNoSpaceDim();

   assert( nsd>=2 && nsd <= 3);

  nrOfWeights = (nsd==2) ? 7 : 10;
  
  weights = new real[nrOfWeights];
  offsets = new int[nrOfWeights];

  real dx2 = geometry->getDx()*geometry->getDx();
  real dy2 = geometry->getDy()*geometry->getDy();
  real dz2 = geometry->getDz()*geometry->getDz();

  real dx = geometry->getDx();
  real dy = geometry->getDy();
  real dz = geometry->getDz();

  if(nsd == 3) {
    weights[0] = -2.0/dx2 - 2.0/dy2 -2.0/dz2 ; 
    weights[1] = weights[2] = 1.0/dx2;
    weights[3] = weights[4] = 1.0/dy2;
    weights[5] = weights[6] = 1.0/dz2;
    weights[7] = 1.0/dx;
    weights[8] = 1.0/dy;
    weights[9] = 1.0/dz;

    offsets[0] = 0;
    offsets[1] = geometry->node( -1, 0, 0 );
    offsets[2] = geometry->node(  1, 0, 0 );
    offsets[3] = geometry->node(  0,-1, 0 );
    offsets[4] = geometry->node(  0, 1, 0 );
    offsets[5] = geometry->node(  0, 0,-1 );
    offsets[6] = geometry->node(  0, 0, 1 );
  }
  else {
    weights[0] = -2.0/dx2 + 1.0/dy2; 
    weights[1] = weights[2] = 1.0/dx2;
    weights[3] = weights[4] = 1.0/dy2;
    weights[5] = 1.0/dx;
    weights[6] = 1.0/dy;


    offsets[0] = 0;
    offsets[1] = geometry->node( 1, 0, 0 );
    offsets[2] = geometry->node( -1, 0, 0 );
    offsets[3] = geometry->node(  0, 2, 0 );
    offsets[4] = geometry->node(  0, 1, 0 );
   

    cout<<" dx = "<<dx<<"  and dy = "<<dy<<endl;

  }  

}


bool WavesWaveEqBoundaryBot::doApply(const double *x,  double *y) const
{
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 
  //cout<<"nsd = "<< geometry->getNoSpaceDim()<<endl;
  // indexes->presentSDIndexes();
  const int nsd = geometry->getNoSpaceDim();
  assert( nsd>=2 && nsd <= 3);

  // We apply y <- (dt^2 Laplace + 2) x - y

  // For faster loops:
  real weights0 = weights[0];
  real weights1 = weights[1];
  real weights2 = weights[2];
  real weights3 = weights[3];
  real weights4 = weights[4];
  real weights5 = weights[5];
  real weights6 = weights[6];
  real weights7 = (nsd == 3)? weights[7] : 0.; 
  real weights8 = (nsd == 3)? weights[8] : 0.;
  real weights9 = (nsd == 3)? weights[9] : 0.;

  int offsets1 = offsets[1];
  int offsets2 = offsets[2];
  int offsets3 = offsets[3];
  int offsets4 = offsets[4];
  int offsets5 = (nsd == 3)? offsets[5] : 0;
  int offsets6 = (nsd == 3)? offsets[6] : 0;

  for(int loop = 0; loop< indexes->nOfLoopIndex(); loop++ )
    {
      int lstart = indexes->loopStart( loop );
      int lstop  = indexes->loopStop( loop );

       cout << "### WavesWaveEqBoundary " << loop << " from " << lstart << " to " 
        << lstop << endl;
 

      if(nsd == 3)
      {
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
      }	  
      else if (nsd==2)
      {
	  for(int n = lstart; n<= lstop; n++ )
	  {


// elliptic equation obtained after applying Laplace transform

// equation (1.3) in my remarks or: laplace w - omega^2 c(x) omega = -delta(x-x_0)
	

	    // elliptic equation after transformation w=ln v
	 
  y[n] =  
  weights0 * x[ n ] +
 weights1 * x[ n+offsets1 ] +
 weights2 * x[ n+offsets2 ] +
 weights3 * x[ n+offsets3 ] -
  2*weights4 * x[ n+offsets4 ] +
   (omega_*omega_)*(weights5*( x[ n+offsets1]  -  x[n])*weights5*(x[ n+offsets1]  - x[n] )
   +   weights6*(x[n+offsets4 ] - x[n])*weights6*(x[ n+offsets4 ] - x[n]));

	    cout<<" in computing bottom boundary( "<<n<<")="<<y[n]<<endl;

	  }    
      }
    }
 
    return true;
    }

//======================================================================================================
//----------------------------------------------------------------------


WavesWaveEqBoundaryLeft::WavesWaveEqBoundaryLeft(WavesSDIndexes *sdi)
  : WavesSDOperator( sdi ), freeIndexes(false)
{
  doInitialize();
}

void WavesWaveEqBoundaryLeft::doInitialize() 
{
//  cout<<"inside doInitialize"<<endl;
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 

  // 1) create Laplace operator
  const int nsd = geometry->getNoSpaceDim();

   assert( nsd>=2 && nsd <= 3);

  nrOfWeights = (nsd==2) ? 7 : 10;
  
  weights = new real[nrOfWeights];
  offsets = new int[nrOfWeights];

  real dx2 = geometry->getDx()*geometry->getDx();
  real dy2 = geometry->getDy()*geometry->getDy();
  real dz2 = geometry->getDz()*geometry->getDz();

  real dx = geometry->getDx();
  real dy = geometry->getDy();
  real dz = geometry->getDz();

  if(nsd == 3) {
    weights[0] = -2.0/dx2 - 2.0/dy2 -2.0/dz2 ; 
    weights[1] = weights[2] = 1.0/dx2;
    weights[3] = weights[4] = 1.0/dy2;
    weights[5] = weights[6] = 1.0/dz2;
    weights[7] = 1.0/dx;
    weights[8] = 1.0/dy;
    weights[9] = 1.0/dz;

    offsets[0] = 0;
    offsets[1] = geometry->node( -1, 0, 0 );
    offsets[2] = geometry->node(  1, 0, 0 );
    offsets[3] = geometry->node(  0,-1, 0 );
    offsets[4] = geometry->node(  0, 1, 0 );
    offsets[5] = geometry->node(  0, 0,-1 );
    offsets[6] = geometry->node(  0, 0, 1 );
  }
  else {
    weights[0] = 1.0/dx2 - 2.0/dy2; 
    weights[1] = weights[2] = 1.0/dx2;
    weights[3] = weights[4] = 1.0/dy2;
    weights[5] = 1.0/dx;
    weights[6] = 1.0/dy;


    offsets[0] = 0;
    offsets[1] = geometry->node( 1, 0, 0 );
    offsets[2] = geometry->node( 2, 0, 0 );
    offsets[3] = geometry->node(  0, 1, 0 );
    offsets[4] = geometry->node(  0, -1, 0 );
   

    cout<<" dx = "<<dx<<"  and dy = "<<dy<<endl;

  }  

}


bool WavesWaveEqBoundaryLeft::doApply(const double *x,  double *y) const
{
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 
  //cout<<"nsd = "<< geometry->getNoSpaceDim()<<endl;
  // indexes->presentSDIndexes();
  const int nsd = geometry->getNoSpaceDim();
  assert( nsd>=2 && nsd <= 3);

  // We apply y <- (dt^2 Laplace + 2) x - y

  // For faster loops:
  real weights0 = weights[0];
  real weights1 = weights[1];
  real weights2 = weights[2];
  real weights3 = weights[3];
  real weights4 = weights[4];
  real weights5 = weights[5];
  real weights6 = weights[6];
  real weights7 = (nsd == 3)? weights[7] : 0.; 
  real weights8 = (nsd == 3)? weights[8] : 0.;
  real weights9 = (nsd == 3)? weights[9] : 0.;

  int offsets1 = offsets[1];
  int offsets2 = offsets[2];
  int offsets3 = offsets[3];
  int offsets4 = offsets[4];
  int offsets5 = (nsd == 3)? offsets[5] : 0;
  int offsets6 = (nsd == 3)? offsets[6] : 0;

  for(int loop = 0; loop< indexes->nOfLoopIndex(); loop++ )
    {
      int lstart = indexes->loopStart( loop );
      int lstop  = indexes->loopStop( loop );

       cout << "### WavesWaveEqBoundary " << loop << " from " << lstart << " to " 
        << lstop << endl;
 

      if(nsd == 3)
      {
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
      }	  
      else if (nsd==2)
      {
	  for(int n = lstart; n<= lstop; n++ )
	  {


// elliptic equation obtained after applying Laplace transform

// equation (1.3) in my remarks or: laplace w - omega^2 c(x) omega = -delta(x-x_0)



  y[n] =  
  weights0 * x[ n ] -
  2*weights1 * x[ n+offsets1 ] +
   weights2 * x[ n+offsets2 ] +
  weights3 * x[ n+offsets3 ] +
  weights4 * x[ n+offsets4 ] +
   (omega_*omega_)*(weights5*(x[n] - x[ n+offsets1])*weights5*(x[n] - x[ n+offsets1])
   +   weights6*(x[n+offsets4 ] - x[n])*weights6*(x[ n+offsets4 ] - x[n]));

	    cout<<" in computing operator ( "<<n<<")="<<y[n]<<endl;

	  
	  }    
      }
    }
 
    return true;
    }

//=====================================================================================
//----------------------------------------------------------------------


WavesWaveEqBoundaryTop::WavesWaveEqBoundaryTop(WavesSDIndexes *sdi)
  : WavesSDOperator( sdi ), freeIndexes(false)
{
  doInitialize();
}

void WavesWaveEqBoundaryTop::doInitialize() 
{
 // cout<<"inside doInitialize"<<endl;
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 

  // 1) create Laplace operator
  const int nsd = geometry->getNoSpaceDim();

   assert( nsd>=2 && nsd <= 3);

  nrOfWeights = (nsd==2) ? 7 : 10;
  
  weights = new real[nrOfWeights];
  offsets = new int[nrOfWeights];

  real dx2 = geometry->getDx()*geometry->getDx();
  real dy2 = geometry->getDy()*geometry->getDy();
  real dz2 = geometry->getDz()*geometry->getDz();

  real dx = geometry->getDx();
  real dy = geometry->getDy();
  real dz = geometry->getDz();

  if(nsd == 3) {
    weights[0] = -2.0/dx2 - 2.0/dy2 -2.0/dz2 ; 
    weights[1] = weights[2] = 1.0/dx2;
    weights[3] = weights[4] = 1.0/dy2;
    weights[5] = weights[6] = 1.0/dz2;
    weights[7] = 1.0/dx;
    weights[8] = 1.0/dy;
    weights[9] = 1.0/dz;

    offsets[0] = 0;
    offsets[1] = geometry->node( -1, 0, 0 );
    offsets[2] = geometry->node(  1, 0, 0 );
    offsets[3] = geometry->node(  0,-1, 0 );
    offsets[4] = geometry->node(  0, 1, 0 );
    offsets[5] = geometry->node(  0, 0,-1 );
    offsets[6] = geometry->node(  0, 0, 1 );
  }
  else {
    weights[0] = -2.0/dx2 + 1.0/dy2; 
    weights[1] = weights[2] = 1.0/dx2;
    weights[3] = weights[4] = 1.0/dy2;
    weights[5] = 1.0/dx;
    weights[6] = 1.0/dy;


    offsets[0] = 0;
    offsets[1] = geometry->node( 1, 0, 0 );
    offsets[2] = geometry->node( -1, 0, 0 );
    offsets[3] = geometry->node(  0, -1, 0 );
    offsets[4] = geometry->node(  0, -2, 0 );
   

    cout<<" dx = "<<dx<<"  and dy = "<<dy<<endl;

  }  

}


bool WavesWaveEqBoundaryTop::doApply(const double *x,  double *y) const
{
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 
  //cout<<"nsd = "<< geometry->getNoSpaceDim()<<endl;
  // indexes->presentSDIndexes();
  const int nsd = geometry->getNoSpaceDim();
  assert( nsd>=2 && nsd <= 3);

  // We apply y <- (dt^2 Laplace + 2) x - y

  // For faster loops:
  real weights0 = weights[0];
  real weights1 = weights[1];
  real weights2 = weights[2];
  real weights3 = weights[3];
  real weights4 = weights[4];
  real weights5 = weights[5];
  real weights6 = weights[6];
  real weights7 = (nsd == 3)? weights[7] : 0.; 
  real weights8 = (nsd == 3)? weights[8] : 0.;
  real weights9 = (nsd == 3)? weights[9] : 0.;

  int offsets1 = offsets[1];
  int offsets2 = offsets[2];
  int offsets3 = offsets[3];
  int offsets4 = offsets[4];
  int offsets5 = (nsd == 3)? offsets[5] : 0;
  int offsets6 = (nsd == 3)? offsets[6] : 0;

  for(int loop = 0; loop< indexes->nOfLoopIndex(); loop++ )
    {
      int lstart = indexes->loopStart( loop );
      int lstop  = indexes->loopStop( loop );

       cout << "### WavesWaveEqBoundary " << loop << " from " << lstart << " to " 
        << lstop << endl;
 

      if(nsd == 3)
      {
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
      }	  
      else if (nsd==2)
      {
	  for(int n = lstart; n<= lstop; n++ )
	  {


// elliptic equation obtained after applying Laplace transform

// equation (1.3) in my remarks or: laplace w - omega^2 c(x) omega = -delta(x-x_0)
	

  y[n] =  
  weights0 * x[ n ] +
  weights1 * x[ n+offsets1 ] +
   weights2 * x[ n+offsets2 ] -
  2*weights3 * x[ n+offsets3 ] +
  weights4 * x[ n+offsets4 ] +
   (omega_*omega_)*(weights5*(x[n] - x[ n+offsets1])*weights5*(x[n] - x[ n+offsets1])
   +   weights6*(x[n] - x[n+offsets3] )*weights6*(x[n] -  x[n+offsets3] ));

  //	    cout<<" in computing top boundary ( "<<n<<")="<<y[n]<<endl;

	  }    
      }
    }
 
    return true;
    }
