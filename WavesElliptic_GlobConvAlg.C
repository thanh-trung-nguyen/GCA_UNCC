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
#include "include/wavesElliptic_GlobConvAlg.h"

//----------------------------------------------------------------------
//
// class WavesElliptic_GlobConvAlgI3
//
//----------------------------------------------------------------------


WavesElliptic_GlobConvAlgI3::WavesElliptic_GlobConvAlgI3(WavesSDGeometry *sdg)
    : WavesSDOperator( new WavesSDInterior(*sdg) ), freeIndexes(true)
{
  doInitialize();
}

WavesElliptic_GlobConvAlgI3::WavesElliptic_GlobConvAlgI3(WavesSDIndexes *sdi)
  : WavesSDOperator( sdi ), freeIndexes(false)
{
  doInitialize();
}

void WavesElliptic_GlobConvAlgI3::doInitialize() 
{
  //cout<<"inside doInitialize"<<endl;
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 

  // 1) create Laplace operator
  const int nsd = geometry->getNoSpaceDim();

   assert( nsd>=2 && nsd <= 3);

  nrOfWeights = (nsd==2) ? 8 : 10;
  
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
    weights[7] = dx;

    offsets[0] = 0;
    offsets[1] = geometry->node( -1, 0, 0 );
    offsets[2] = geometry->node(  1, 0, 0 );
    offsets[3] = geometry->node(  0,-1, 0 );
    offsets[4] = geometry->node(  0, 1, 0 );
  }  

}


bool WavesElliptic_GlobConvAlgI3::doApply(const double *x,  double *y) const
{
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 
  //cout<<"nsd = "<< geometry->getNoSpaceDim()<<endl;
  // indexes->presentSDIndexes();
  const int nsd = geometry->getNoSpaceDim();
  double  del=0.0;
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
  real weights7 = weights[7];
 
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
	
	    //==================0
  for(int n = lstart; n<= lstop; n++ )
	      {
  
	    // case for q_1
	    //   (i-1,j) =  offsets1 
		// (k,s) = s 
		// (i+1, j) = offsets2
		// (i, j-1) = offsets3
		// (i, j+1) = offsets4 

  for(int loop1 = 0; loop1< indexes->nOfLoopIndex(); loop1++ )
    {
      int lstart1 = indexes->loopStart( loop1 );
      int lstop1  = indexes->loopStop( loop1 );

 cout << "### GlobConvAlg inner loop s, number of loop is " << loop1 << " from " << lstart1 << " to " 
       << lstop1 << endl;
 
    
		for(int s = lstart1; s <= lstop1; s++ )
		  {		

		    if (s == n+offsets1 || s == n+offsets2 || s == n+offsets3 || s==n+offsets4 || s==n)
		      { 



			// grad_q_i is vector with values < 0  which should be prepaired 

			 del = weights[5]*(delta(s, n + offsets1) - 4*delta(s,n) + delta(s, n+offsets2) + delta(s, n+offsets3) + delta(s, n+offsets4)) + grad_q_i_[n]*(delta(s,n+offsets2) + delta(s,n+offsets4) + 2*delta(s,n));


			/*
			if (del != 0.0)
			  cout<<" node s = "<<s<<"  and n = "<<n<<"  delta(s, n+offsets1)="<<delta(s, n + offsets1)<<" 4*delta(s,n) ="<<4*delta(s,n)<<"   delta(s, n+offsets2)="<<delta(s, n+offsets2)<<"     delta(s, n+offsets3)= "<<delta(s, n+offsets3)<<"     delta(s, n+offsets4) ="<<delta(s, n+offsets4)<<endl;
				
			*/
cout<<"  previous value y(s) = y("<<s<<") = "<<y[s]<<endl;
			

 y[s] +=   grad_q_square_[n]*del;   

	cout<<"  y(s) = y("<<s<<") = "<<y[s]<<endl;
			
			if (s==n)
			  cout<<" ***node s= "<<s<<"  and node n="<<n<<" 4*delta(s,n)="<<4*delta(s,n)<<"  del "<<del<<endl;		


      }
		    
  //    func_q_i_[n]*(weights5*(x[ n+offsets2 ] - x[n]) + weights6*(x[n+offsets4 ] - x[n]));
		
		  }  // for s
    }     // for loop 1
	      } // for n
      } // for nsd==2
    } // for outer loop
 //==============================================0

 
    return true;
    }


//----------------------------------------------------------------------
//
// class WavesElliptic_GlobConvAlgI2
//
//----------------------------------------------------------------------


WavesElliptic_GlobConvAlgI2::WavesElliptic_GlobConvAlgI2(WavesSDGeometry *sdg)
    : WavesSDOperator( new WavesSDInterior(*sdg) ), freeIndexes(true)
{
  doInitialize();
}

WavesElliptic_GlobConvAlgI2::WavesElliptic_GlobConvAlgI2(WavesSDIndexes *sdi)
  : WavesSDOperator( sdi ), freeIndexes(false)
{
  doInitialize();
}

void WavesElliptic_GlobConvAlgI2::doInitialize() 
{
  //cout<<"inside doInitialize"<<endl;
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 

  // 1) create Laplace operator
  const int nsd = geometry->getNoSpaceDim();

   assert( nsd>=2 && nsd <= 3);

  nrOfWeights = (nsd==2) ? 8 : 10;
  
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
    weights[7] = dx;

    offsets[0] = 0;
    offsets[1] = geometry->node( -1, 0, 0 );
    offsets[2] = geometry->node(  1, 0, 0 );
    offsets[3] = geometry->node(  0,-1, 0 );
    offsets[4] = geometry->node(  0, 1, 0 );
  }  

}


bool WavesElliptic_GlobConvAlgI2::doApply(const double *x,  double *y) const
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
  real weights7 = weights[7];
 
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
	
	    //==================0
  for(int n = lstart; n<= lstop; n++ )
	      {

	    // case for q_1
	    //   (i-1,j) =  offsets1 
		// (k,s) = s 
		// (i+1, j) = offsets2
		// (i, j-1) = offsets3
		// (i, j+1) = offsets4 

  for(int loop1 = 0; loop1< indexes->nOfLoopIndex(); loop1++ )
    {
      int lstart1 = indexes->loopStart( loop1 );
      int lstop1  = indexes->loopStop( loop1 );

 cout << "### GlobConvAlg inner loop s, number of loop is " << loop1 << " from " << lstart1 << " to " 
       << lstop1 << endl;
 
    
		for(int s = lstart1; s <= lstop1; s++ )
		  {		

		    if (s == n+offsets1 || s == n+offsets2 || s == n+offsets3 || s==n+offsets4 || s==n)
		      { 



			// grad_q_i is vector with values < 0  which should be prepaired 

			double del = weights5*(delta(s, n + offsets1) - 4*delta(s,n) + delta(s, n+offsets2) + delta(s, n+offsets3) + delta(s, n+offsets4)) + grad_q_i_[n]*(delta(s,n+offsets2) + delta(s,n+offsets4) + 2*delta(s,n));


			/*
			if (del != 0.0)
			  cout<<" node s = "<<s<<"  and n = "<<n<<"  delta(s, n+offsets1)="<<delta(s, n + offsets1)<<" 4*delta(s,n) ="<<4*delta(s,n)<<"   delta(s, n+offsets2)="<<delta(s, n+offsets2)<<"     delta(s, n+offsets3)= "<<delta(s, n+offsets3)<<"     delta(s, n+offsets4) ="<<delta(s, n+offsets4)<<endl;
				
			*/
cout<<"  previous value y(s) = y("<<s<<") = "<<y[s]<<endl;
			

 y[s] +=  grad_q_i_CWF_[n]*(weights5*(x[ n+offsets2 ] - x[n]) + weights6*(x[n+offsets4 ] - x[n]))*del;
	     

	cout<<"  y(s) = y("<<s<<") = "<<y[s]<<endl;
			
			if (s==n)
			  cout<<" ***node s= "<<s<<"  and node n="<<n<<" 4*delta(s,n)="<<4*delta(s,n)<<"  del "<<del<<endl;		


      }
		    
  //    func_q_i_[n]*(weights5*(x[ n+offsets2 ] - x[n]) + weights6*(x[n+offsets4 ] - x[n]));
		
		  }  // for s
    }     // for loop 1
	      } // for n
      } // for nsd==2
    } // for outer loop
 //==============================================0

 
    return true;
    }


//----------------------------------------------------------------------
//
// class WavesElliptic_GlobConvAlgI1
//
//----------------------------------------------------------------------


WavesElliptic_GlobConvAlgI1::WavesElliptic_GlobConvAlgI1(WavesSDGeometry *sdg)
    : WavesSDOperator( new WavesSDInterior(*sdg) ), freeIndexes(true)
{
  doInitialize();
}

WavesElliptic_GlobConvAlgI1::WavesElliptic_GlobConvAlgI1(WavesSDIndexes *sdi)
  : WavesSDOperator( sdi ), freeIndexes(false)
{
  doInitialize();
}

void WavesElliptic_GlobConvAlgI1::doInitialize() 
{
  //cout<<"inside doInitialize"<<endl;
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 

  // 1) create Laplace operator
  const int nsd = geometry->getNoSpaceDim();

   assert( nsd>=2 && nsd <= 3);

  nrOfWeights = (nsd==2) ? 8 : 10;
  
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
    weights[7] = dx;

    offsets[0] = 0;
    offsets[1] = geometry->node( -1, 0, 0 );
    offsets[2] = geometry->node(  1, 0, 0 );
    offsets[3] = geometry->node(  0,-1, 0 );
    offsets[4] = geometry->node(  0, 1, 0 );
  }  

}


bool WavesElliptic_GlobConvAlgI1::doApply(const double *x,  double *y) const
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
  real weights7 = weights[7];
 
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
	
	    //==================0
  for(int n = lstart; n<= lstop; n++ )
	      {

	    // case for q_1
	    //   (i-1,j) =  offsets1 
		// (k,s) = s 
		// (i+1, j) = offsets2
		// (i, j-1) = offsets3
		// (i, j+1) = offsets4 

  for(int loop1 = 0; loop1< indexes->nOfLoopIndex(); loop1++ )
    {
      int lstart1 = indexes->loopStart( loop1 );
      int lstop1  = indexes->loopStop( loop1 );

 cout << "### GlobConvAlg inner loop s, number of loop is " << loop1 << " from " << lstart1 << " to " 
       << lstop1 << endl;
 
    
		for(int s = lstart1; s <= lstop1; s++ )
		  {		

		    if (s == n+offsets1 || s == n+offsets2 || s == n+offsets3 || s==n+offsets4 || s==n)
		      { 



			// grad_q_i is vector with values < 0  which should be prepaired 

			double del = delta(s, n + offsets1) - 4*delta(s,n) + delta(s, n+offsets2) + delta(s, n+offsets3) + delta(s, n+offsets4) + grad_q_i_[n]*(delta(s,n+offsets2) + delta(s,n+offsets4) + 2*delta(s,n));
			

			/*
			if (del != 0.0)
			  cout<<" node s = "<<s<<"  and n = "<<n<<"  delta(s, n+offsets1)="<<delta(s, n + offsets1)<<" 4*delta(s,n) ="<<4*delta(s,n)<<"   delta(s, n+offsets2)="<<delta(s, n+offsets2)<<"     delta(s, n+offsets3)= "<<delta(s, n+offsets3)<<"     delta(s, n+offsets4) ="<<delta(s, n+offsets4)<<endl;
				
			*/
			cout<<"  previous value y(s) = y("<<s<<") = "<<y[s]<<"  grad_q_i "<< grad_q_i_[n]<<"  n  "<<n<<endl;
			

 y[s] +=  
   (weights0 * x[ n ] +
    weights1 * x[ n+offsets1 ] +
    weights2 * x[ n+offsets2 ] +
    weights3 * x[ n+offsets3 ] +
    weights4 * x[ n+offsets4 ] )*del;
    	     

	cout<<"  y(s) = y("<<s<<") = "<<y[s]<<endl;
			
			if (s==n)
			  cout<<" ***node s= "<<s<<"  and node n="<<n<<" 4*delta(s,n)="<<4*delta(s,n)<<"  del "<<del<<endl;		


      }
		    
  //    func_q_i_[n]*(weights5*(x[ n+offsets2 ] - x[n]) + weights6*(x[n+offsets4 ] - x[n]));
		
		  }  // for s
    }     // for loop 1
	      } // for n
      } // for nsd==2
    } // for outer loop
 //==============================================0

 
    return true;
    }



//----------------------------------------------------------------------
//
// class WavesElliptic_GlobConvAlg_CWF
//
//----------------------------------------------------------------------


WavesElliptic_GlobConvAlg_CWF::WavesElliptic_GlobConvAlg_CWF(WavesSDGeometry *sdg, int code)
    : WavesSDOperator( new WavesSDInterior(*sdg, code) ), freeIndexes(true)
{
  doInitialize();
}

WavesElliptic_GlobConvAlg_CWF::WavesElliptic_GlobConvAlg_CWF(WavesSDGeometry *sdg)
    : WavesSDOperator( new WavesSDInterior(*sdg) ), freeIndexes(true)
{
  doInitialize();
}


WavesElliptic_GlobConvAlg_CWF::WavesElliptic_GlobConvAlg_CWF(WavesSDIndexes *sdi)
  : WavesSDOperator( sdi ), freeIndexes(false)
{
  doInitialize();
}

void WavesElliptic_GlobConvAlg_CWF::doInitialize() 
{
  //cout<<"inside doInitialize"<<endl;
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


bool WavesElliptic_GlobConvAlg_CWF::doApply(const double *x,  double *y) const
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


     
      //  cout << "### GlobConvAlg outer loop n, number of loop is " << loop << " from " << lstart << " to " 
      // << lstop << endl;
 
 
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
	                                     
	  }
      }	  
      else if (nsd==2)
      {
	//	for(int s = lstart1; s <= lstop1; s++ )
	// {
	    for(int n = lstart; n<= lstop; n++ )
	      {

	    // case for q_1
	    //   (i-1,j) =  offsets1 
		// (k,s) = s 
		// (i+1, j) = offsets2
		// (i, j-1) = offsets3
		// (i, j+1) = offsets4 

  for(int loop1 = 0; loop1< indexes->nOfLoopIndex(); loop1++ )
    {
      int lstart1 = indexes->loopStart( loop1 );
      int lstop1  = indexes->loopStop( loop1 );

      // cout << "### GlobConvAlg inner loop s, number of loop is " << loop1 << " from " << lstart1 << " to " 
      // << lstop1 << endl;
 
    
		for(int s = lstart1; s <= lstop1; s++ )
		  {		

		    if (s == n+offsets1 || s == n+offsets2 || s == n+offsets3 || s==n+offsets4 || s==n)
		      { 

			double del = delta(s, n + offsets1) - 4*delta(s,n) + delta(s, n+offsets2) + delta(s, n+offsets3) + delta(s, n+offsets4);


			/*
			if (del != 0.0)
			  cout<<" node s = "<<s<<"  and n = "<<n<<"  delta(s, n+offsets1)="<<delta(s, n + offsets1)<<" 4*delta(s,n) ="<<4*delta(s,n)<<"   delta(s, n+offsets2)="<<delta(s, n+offsets2)<<"     delta(s, n+offsets3)= "<<delta(s, n+offsets3)<<"     delta(s, n+offsets4) ="<<delta(s, n+offsets4)<<endl;
				
			*/
cout<<"  previous value y(s) = y("<<s<<") = "<<y[s]<<endl;
			
			y[s] +=   
			  (weights0 * x[ n ] +
			   weights1 * x[ n+offsets1 ] +
			   weights2 * x[ n+offsets2 ] +
			   weights3 * x[ n+offsets3 ] +
			   weights4 * x[ n+offsets4 ])*del;

			//	cout<<"  y(s) = y("<<s<<") = "<<y[s]<<endl;
			
			//	if (s==n)
			//  cout<<" ***node s= "<<s<<"  and node n="<<n<<" 4*delta(s,n)="<<4*delta(s,n)<<"  del "<<del<<endl;		


      }
		    
  //    func_q_i_[n]*(weights5*(x[ n+offsets2 ] - x[n]) + weights6*(x[n+offsets4 ] - x[n]));
		
		  }  // for s
    }     // for loop 1
	      } // for n
      } // for nsd==2
    } // for outer loop
  
  return true;
    }


// function assign value at the left boundary to the all nodes in the domain such that
// we initialize like plane wave with const material type inside domain

void WavesElliptic_GlobConvAlg_CWF::doInitGuess(const double *x,  double *y) 
{

 const WavesSDGeometry *geometry = &indexes->getGeometry(); 
  //cout<<"nsd = "<< geometry->getNoSpaceDim()<<endl;
  // indexes->presentSDIndexes();
  const int nsd = geometry->getNoSpaceDim();
  assert( nsd>=2 && nsd <= 3);

  bool doIncludeCorners=1;
  WavesSDBoundary leftBoundary(*geometry, SDiLow, doIncludeCorners);
 
  for(int loop = 0; loop< indexes->nOfLoopIndex(); loop++ )
       {
	 int lstart = indexes->loopStart( loop );
	 int lstop  = indexes->loopStop ( loop );
	 

for(int n = lstart; n<= lstop; n++ )
	   {

	     /*
	       for(int loop_left = 0; loop_left<  leftBoundary.nOfLoopIndex(); loop_left++ )
	       {
	       int left_start = leftBoundary.loopStart( loop_left );
	       int left_stop  = leftBoundary.loopStop( loop_left );
	       
	       cout<<" left boundary  left_start"<<left_start<<"  lstop"<<left_stop<<endl;
	       for(int l = left_start;  l<= left_stop; l++ )
	       {
	       
	       y[n] =   x[l]; 
	       cout<<" y("<<n<<"), node on the left boundary "<<l<<"  node inside domain"<<n<<" value"<<x[l]<<endl;
	       }   
	       }
	     */

	     // assign value at the left boundary
	     y[n] =   x[lstart+offsets[1]]; 
	     cout<<" y("<<n<<"), node on the left boundary "<<lstart<<"  node inside domain"<<n<<" value"<<x[lstart]<<endl;

	   }
       }
  
}

//==================================================================
//----------------------------------------------------------------------
//
// class WavesElliptic_GlobConvAlg
//
//----------------------------------------------------------------------


WavesElliptic_GlobConvAlg::WavesElliptic_GlobConvAlg(WavesSDGeometry *sdg, int code)
    : WavesSDOperator( new WavesSDInterior(*sdg, code) ), freeIndexes(true)
{
  doInitialize();
}

WavesElliptic_GlobConvAlg::WavesElliptic_GlobConvAlg(WavesSDGeometry *sdg)
    : WavesSDOperator( new WavesSDInterior(*sdg) ), freeIndexes(true)
{
  doInitialize();
}


WavesElliptic_GlobConvAlg::WavesElliptic_GlobConvAlg(WavesSDIndexes *sdi)
  : WavesSDOperator( sdi ), freeIndexes(false)
{
  doInitialize();
}

void WavesElliptic_GlobConvAlg::doInitialize() 
{
  //cout<<"inside doInitialize"<<endl;
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


bool WavesElliptic_GlobConvAlg::doApply(const double *x,  double *y) const
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

 
     
       cout << "### GlobConvAlg outer loop n, number of loop is " << loop << " from " << lstart << " to " 
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
	                                     
	  }
      }	  
      else if (nsd==2)
      {
	//	for(int s = lstart1; s <= lstop1; s++ )
	// {
	    for(int n = lstart; n<= lstop; n++ )
	      {

	    // case for q_1
	    //   (i-1,j) =  offsets1 
		// (k,s) = s 
		// (i+1, j) = offsets2
		// (i, j-1) = offsets3
		// (i, j+1) = offsets4 

  for(int loop1 = 0; loop1< indexes->nOfLoopIndex(); loop1++ )
    {
      int lstart1 = indexes->loopStart( loop1 );
      int lstop1  = indexes->loopStop( loop1 );

 cout << "### GlobConvAlg inner loop s, number of loop is " << loop1 << " from " << lstart1 << " to " 
       << lstop1 << endl;
 
    
		for(int s = lstart1; s <= lstop1; s++ )
		  {		

		    if (s == n+offsets1 || s == n+offsets2 || s == n+offsets3 || s==n+offsets4 || s==n)
		      { 



		

			double del =  weights1*delta(s, n + offsets1) - 4*weights0*delta(s,n) +  weights2*delta(s, n+offsets2) + weights3*delta(s, n+offsets3) +  weights4*delta(s, n+offsets4);


			/*
			if (del != 0.0)
			  cout<<" node s = "<<s<<"  and n = "<<n<<"  delta(s, n+offsets1)="<<delta(s, n + offsets1)<<" 4*delta(s,n) ="<<4*delta(s,n)<<"   delta(s, n+offsets2)="<<delta(s, n+offsets2)<<"     delta(s, n+offsets3)= "<<delta(s, n+offsets3)<<"     delta(s, n+offsets4) ="<<delta(s, n+offsets4)<<endl;
				
			*/
cout<<"  previous value y(s) = y("<<s<<") = "<<y[s]<<endl;
			
			y[n] =   
			  (weights0 * x[ n ] +
			   weights1 * x[ n+offsets1 ] +
			   weights2 * x[ n+offsets2 ] +
			   weights3 * x[ n+offsets3 ] +
			   weights4 * x[ n+offsets4 ])*del;

			cout<<"  y(s) = y("<<s<<") = "<<y[s]<<endl;
			
			if (s==n)
			  cout<<" ***node s= "<<s<<"  and node n="<<n<<" 4*delta(s,n)="<<4*delta(s,n)<<"  del "<<del<<endl;		


      }
		    
  //    func_q_i_[n]*(weights5*(x[ n+offsets2 ] - x[n]) + weights6*(x[n+offsets4 ] - x[n]));
		
		  }  // for s
    }     // for loop 1
	      } // for n
      } // for nsd==2
    } // for outer loop
  
  return true;
    }


// function assign value at the left boundary to the all nodes in the domain such that
// we initialize like plane wave with const material type inside domain

void WavesElliptic_GlobConvAlg::doInitGuess(const double *x,  double *y) 
{

 const WavesSDGeometry *geometry = &indexes->getGeometry(); 
  //cout<<"nsd = "<< geometry->getNoSpaceDim()<<endl;
  // indexes->presentSDIndexes();
  const int nsd = geometry->getNoSpaceDim();
  assert( nsd>=2 && nsd <= 3);

  bool doIncludeCorners=1;
  WavesSDBoundary leftBoundary(*geometry, SDiLow, doIncludeCorners);
 
  for(int loop = 0; loop< indexes->nOfLoopIndex(); loop++ )
       {
	 int lstart = indexes->loopStart( loop );
	 int lstop  = indexes->loopStop ( loop );
	 

for(int n = lstart; n<= lstop; n++ )
	   {

	     /*
	       for(int loop_left = 0; loop_left<  leftBoundary.nOfLoopIndex(); loop_left++ )
	       {
	       int left_start = leftBoundary.loopStart( loop_left );
	       int left_stop  = leftBoundary.loopStop( loop_left );
	       
	       cout<<" left boundary  left_start"<<left_start<<"  lstop"<<left_stop<<endl;
	       for(int l = left_start;  l<= left_stop; l++ )
	       {
	       
	       y[n] =   x[l]; 
	       cout<<" y("<<n<<"), node on the left boundary "<<l<<"  node inside domain"<<n<<" value"<<x[l]<<endl;
	       }   
	       }
	     */

	     // assign value at the left boundary
	     y[n] =   x[lstart+offsets[1]]; 
	     cout<<" y("<<n<<"), node on the left boundary "<<lstart<<"  node inside domain"<<n<<" value"<<x[lstart]<<endl;

	   }
       }
  
}

//*************************************************************************

//----------------------------------------------------------------------
//
// class WavesElliptic_gradient
//
//----------------------------------------------------------------------


WavesElliptic_gradient::WavesElliptic_gradient(WavesSDGeometry *sdg, int code)
    : WavesSDOperator( new WavesSDInterior(*sdg, code) ), freeIndexes(true)
{
  doInitialize();
}

WavesElliptic_gradient::WavesElliptic_gradient(WavesSDGeometry *sdg)
    : WavesSDOperator( new WavesSDInterior(*sdg) ), freeIndexes(true)
{
  doInitialize();
}


WavesElliptic_gradient::WavesElliptic_gradient(WavesSDIndexes *sdi)
  : WavesSDOperator( sdi ), freeIndexes(false)
{
  doInitialize();
}

void WavesElliptic_gradient::doInitialize() 
{
  //cout<<"inside doInitialize"<<endl;
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


bool WavesElliptic_gradient::doApply(const double *x,  double *y) const
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


     
      //  cout << "### GlobConvAlg outer loop n, number of loop is " << loop << " from " << lstart << " to " 
      // << lstop << endl;
 
 
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
	                                     
	  }
      }	  
      else if (nsd==2)
      {
	//	for(int s = lstart1; s <= lstop1; s++ )
	// {
	    for(int n = lstart; n<= lstop; n++ )
	      {

	    // case for q_1
	    //   (i-1,j) =  offsets1 
		// (k,s) = s 
		// (i+1, j) = offsets2
		// (i, j-1) = offsets3
		// (i, j+1) = offsets4 

  for(int loop1 = 0; loop1< indexes->nOfLoopIndex(); loop1++ )
    {
      int lstart1 = indexes->loopStart( loop1 );
      int lstop1  = indexes->loopStop( loop1 );

      // cout << "### GlobConvAlg inner loop s, number of loop is " << loop1 << " from " << lstart1 << " to " 
      // << lstop1 << endl;
 
    
		for(int s = lstart1; s <= lstop1; s++ )
		  {		

		    if ( s == n+offsets2  || s==n+offsets4 || s==n)
		      { 
		
			double del =	weights5*delta(s, n+offsets2) - weights5*delta(s,n) - weights6*delta(s,n) +  weights6*delta(s, n+offsets4);

		
			/*
			if (del != 0.0)
			  cout<<" node s = "<<s<<"  and n = "<<n<<"  delta(s, n+offsets1)="<<delta(s, n + offsets1)<<" 4*delta(s,n) ="<<4*delta(s,n)<<"   delta(s, n+offsets2)="<<delta(s, n+offsets2)<<"     delta(s, n+offsets3)= "<<delta(s, n+offsets3)<<"     delta(s, n+offsets4) ="<<delta(s, n+offsets4)<<endl;
				
			*/
cout<<"  previous value y(s) = y("<<s<<") = "<<y[s]<<endl;
 
                       y[s] +=   (weights5*(x[ n+offsets2 ] - x[n]) +
	                          weights6*(x[ n+offsets4 ] - x[n]))*del;
			//	cout<<"  y(s) = y("<<s<<") = "<<y[s]<<endl;
			
			//	if (s==n)
			//  cout<<" ***node s= "<<s<<"  and node n="<<n<<" 4*delta(s,n)="<<4*delta(s,n)<<"  del "<<del<<endl;		


      }
		    
  //    func_q_i_[n]*(weights5*(x[ n+offsets2 ] - x[n]) + weights6*(x[n+offsets4 ] - x[n]));
		
		  }  // for s
    }     // for loop 1
	      } // for n
      } // for nsd==2
    } // for outer loop
  
  return true;
    }


//  compute square of the gradient of the function x  and assign it to y
void WavesElliptic_gradient::compute_square_gradient(const double *x,  double *y) 
{

 const WavesSDGeometry *geometry = &indexes->getGeometry(); 
  //cout<<"nsd = "<< geometry->getNoSpaceDim()<<endl;
  // indexes->presentSDIndexes();
  const int nsd = geometry->getNoSpaceDim();
  assert( nsd>=2 && nsd <= 3);
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
	 int lstop  = indexes->loopStop ( loop );
	 

for(int n = lstart; n<= lstop; n++ )
	   {

	   
	     // compute gradient of function
	     y[n] =   (weights5*(x[ n+offsets2 ] - x[n])*weights5*(x[ n+offsets2 ] - x[n])
		       +   weights6*(x[n+offsets4 ] - x[n])*weights6*(x[ n+offsets4 ] - x[n]));

	     cout<<" y("<<n<<"), node on the left boundary "<<lstart<<"  node inside domain"<<n<<" value"<<x[lstart]<<endl;

	   }
       }
  
}

//  compute   gradient of the function x  and assign it to y
void WavesElliptic_gradient::compute_gradient(const double *x,  double *y) 
{

 const WavesSDGeometry *geometry = &indexes->getGeometry(); 
  //cout<<"nsd = "<< geometry->getNoSpaceDim()<<endl;
  // indexes->presentSDIndexes();
  const int nsd = geometry->getNoSpaceDim();
  assert( nsd>=2 && nsd <= 3);
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
	 int lstop  = indexes->loopStop ( loop );
	 

for(int n = lstart; n<= lstop; n++ )
	   {

	   
	     // compute gradient of function
	     y[n] =   (weights5*(x[ n+offsets2 ] - x[n]) +
		       weights6*(x[ n+offsets4 ] - x[n]));

	     cout<<" y("<<n<<"), node on the left boundary "<<lstart<<"  node inside domain"<<n<<" value"<<x[lstart]<<endl;

	   }
       }
  
}

//  compute   gradient of the function x  and assign it to y
void WavesElliptic_gradient::compute_gradient(const double *x,
					      double *y_1,
					      double *y_2,
					      double *y_3  ) 
{
  
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 
  //cout<<"nsd = "<< geometry->getNoSpaceDim()<<endl;
  // indexes->presentSDIndexes();
  const int nsd = geometry->getNoSpaceDim();
  assert( nsd>=2 && nsd <= 3);
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
      int lstop  = indexes->loopStop ( loop );
      
      
      for(int n = lstart; n<= lstop; n++ )
	{
	  
	  if (nsd ==2)
	    { 
	      // compute gradient of function
	      y_1[n] =   weights5*(x[ n+offsets2 ] - x[n]);
	      
	      y_2[n] =  weights6*(x[ n+offsets4 ] - x[n]);
	      
	      
	    }
	  else if (nsd == 3)
	    {
	      // compute gradient of function
	      y_1[n] = weights7*(x[n+offsets2] - x[n]);
	      y_2[n] =   weights8*(x[n+offsets4] - x[n]);
	      y_3[n] =  weights9*(x[n+offsets6] - x[n]);	 
	      
	    }
	  
	  cout<<" y("<<n<<"), node on the left boundary "<<lstart<<"  node inside domain"<<n<<" value"<<x[lstart]<<endl;
	  
	}
    }
  
}


//  compute   gradient of the function x  and assign it to y_1, y_2, y_3 using MV_Vector type. added by Thanh
void WavesElliptic_gradient::compute_gradient(MV_Vector<double> x, MV_Vector<double>& y_1, MV_Vector<double>& y_2, MV_Vector<double>& y_3) 
{
  
  const WavesSDGeometry *geometry = &indexes->getGeometry(); 
  //cout<<"nsd = "<< geometry->getNoSpaceDim()<<endl;
  // indexes->presentSDIndexes();
  const int nsd = geometry->getNoSpaceDim();
  assert( nsd>=2 && nsd <= 3);
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
      int lstop  = indexes->loopStop ( loop );
      
      
      for(int n = lstart; n<= lstop; n++ )
	{
	  
	  if (nsd ==2)
	    { 
	      // compute gradient of function
	      y_1(n) =  weights5*(x(n+offsets2) - x(n));
	      
	      y_2(n) =  weights6*(x(n+offsets4) - x(n));
	      
	      
	    }
	  else if (nsd == 3)
	    {
	      // compute gradient of function
	      y_1(n) =  weights7*(x(n+offsets2) - x(n));
	      y_2(n) =  weights8*(x(n+offsets4) - x(n));
	      y_3(n) =  weights9*(x(n+offsets6) - x(n));	 
	      
	    }	  
	}
    }
  
}


//==================================================================
