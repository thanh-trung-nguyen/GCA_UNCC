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



#include <iostream>
#include <string.h>
#include <fstream>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "include/wavesPoissonEqOp.h"
#include "include/wavesFEcomp.h"
#include"include/waveskraftwerk.hh"
#include <sys/times.h>
#include "include/wavesSDIndexes.h"
#include "include/wavesFuncDefs.h"
//#include "NeumannBC1.h"
#include "include/wavesWaveEquationOp.h"
#include <petscis.h>

#if PETSC_USE_DEBUG
#define CHKERRA(e) if(e){PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,e,0,0);exit(1);}
#else
#define CHKERRA(e); /* obsolete in petsc 2.1 */
#endif 
#define SETERRA(n,p,s)  SETERRQ(n,s)  /* obsolete in petsc 2.1 */


typedef double real;
typedef MV_ColMat<real> Mat_real;

// two functions to use in the debugger 
static real get_ij(Mat_real &A, int i, int j){
  return A(i,j);
}

static real get_i(MV_Vector <real> &A, int i)
{
  return A(i);
}

WavesPoissonEqOp::WavesPoissonEqOp(Grid &gg_,
			 bool& USE_FEM_,
			 bool &USE_RHS_,
			 bool &USE_DIRICHLET_FEM_,
			 bool &PRINT_FILES_,
			 double &EVALPT_X_,
			 double &EVALPT_Y_,double &EVALPT_Z_,
			 double &rhs_, int& noQpts_) 
  :  opt(gg_)
{
 
  gg = gg_;
  //cout<<gg.getNoNodes()<<endl;
  
  USE_FEM =  USE_FEM_;

  USE_RHS = USE_RHS_;
 
  USE_DIRICHLET_FEM = USE_DIRICHLET_FEM_;

  PRINT_FILES = PRINT_FILES_;
 
  EVALPT_X = EVALPT_X_;
  EVALPT_Y = EVALPT_Y_;
  EVALPT_Z = EVALPT_Z_;
 
  rhs = rhs_;
  noQpts = noQpts_;
 // cout<<"noqpts "<<noQpts<<endl;

  bvalues_initialised=0;
  Poisson_FEM_initialized = 0;
  nbn = 0;
  // nsd = gg.getNoSpaceDim();


 
} 


  
WavesPoissonEqOp::~WavesPoissonEqOp()
{
  if(bvalues_initialised){
    //  ierr = PetscFree( boundindex );
    delete[] boundindex;
    delete[] bound_values;
    delete[] bvalues;
    // ierr = PetscFree( bvalues );
  }

 if(Poisson_FEM_initialized){
   ierr = VecDestroy(puas_sol); CHKERRA(ierr);
   ierr = MatDestroy(A); CHKERRA(ierr);
   //delete[] boundindex;
   //delete[]  bound_values;
   //delete[] bvalues;
 }

}


const int DEBUG=1;



void WavesPoissonEqOp::Configurate(char *fname)

{
 const int maxlen = 30;
 char var[maxlen];
 char fname1[maxlen];

 double dx,dy,dz,x,y,z;

 char text1[20],text2[20],text3[20],text4[20],text5[20],text6[20],text7[20],text8[20],text9[20],text10[20],text11[20],text12[20];
 int dbg=1;
 ifstream inp;
 inp.open(fname);   
 inp>>text1>>nsd;
 inp>>fname1;
 inp>>text6>>USE_RHS;
 inp>>text8>>USE_DIRICHLET_FEM;
 inp>>text10>>PRINT_FILES;
 inp>>text11>>EVALPT_X;
 inp>>text12>>EVALPT_Y;

 //if(dbg)cout<<"Dimensions:"<<nsd<<endl;   


 gg.scan(fname1);
 inp.close();
}



int WavesPoissonEqOp::dirichletBCxxx(Mat&  Matrice )
{
  PetscFunctionBegin;
  int i,j,neq;
  int *bndNodesfem;  
  double *diriVal;
  PetscScalar *tmpvec, one = 1.0, zero = 0.0, minusone = -1.0;
  IS bndRows;

  MV_Vector<int> markNodes;
  opt.findBoundaryNodes(gg,markNodes); 
  MV_Vector<int> bndNodes;
  opt.makeSortedNodesIndex(markNodes,bndNodes,1);
  int nbn = bndNodes.size();
 
  PetscInt nbn_petsc=nbn;

  boundindex = new int[nbn];
  bound_values = new double[nbn];
  
//initialize initial values = 0.0 in the all points

  for( int i=0;i<nbn;i++)
 {
   boundindex[i] = bndNodes(i);
   if (USE_RHS)
   bound_values[i] = 0.0;
   else
     bound_values[i] = 1.0;

  } 
 
  PetscMalloc((nbn_petsc)*sizeof(PetscInt), &bndNodesfem);

  // (int *) PetscMalloc((gg.getNoNodes()+1)*sizeof(int), bndNodesfem);


  for( i = 0; i < nbn; i++)
    {
      bndNodesfem[ i] = boundindex[i];
      //      cout<<"bndNodesfem = "<<bndNodesfem[i]<<endl;
    }      
  
  // create general index set
  ierr = ISCreateGeneral(PETSC_COMM_SELF, nbn_petsc, bndNodesfem,  &bndRows);
  
//  for( i = 0; i < nbn; i++)
//      cout<<"bndNodesfem = "<<bndNodesfem[i]<<endl;
    
 
  ierr = PetscFree( bndNodesfem );

//  cout<<"before matzerorows"<<endl;

  ierr = MatZeroRowsIS(Matrice, bndRows, one); CHKERRQ(ierr);

  ierr = ISDestroy( bndRows ); 
 
  PetscFunctionReturn(0);
}




void WavesPoissonEqOp::getTriangleQuadrature(const int noQpts, 
						   MV_ColMat<real>& points, 
						   MV_Vector<real>& weights)
{
  // should also check that the size is correct??
  if ( noQpts==3 ) {
    weights=1.0/3.0;
    /*
    points(0,0)=0.66666666666666667;
    points(0,1)=0.16666666666666667;
    points(0,2)=0.16666666666666667;
      
    points(1,0)=0.16666666666666667;
    points(1,1)=0.66666666666666667;
    points(1,2)=0.16666666666666667;
      
    points(2,0)=0.16666666666666667;
    points(2,1)=0.16666666666666667;
    points(2,2)=0.66666666666666667;
    */
     
    points(0,0)=1.0;
    points(0,1)=0.0;
    points(0,2)=0.0;
      
    points(1,0)=0.0;
    points(1,1)=1.0;
    points(1,2)=0.0;
      
    points(2,0)=0.0;
    points(2,1)=0.0;
    points(2,2)=1.0;
    


  }
  else if (noQpts==6){
    weights(0) = 0.223381589678011;
    weights(1) = 0.223381589678011;
    weights(2) = 0.223381589678011;
    weights(3) = 0.109951743655322;
    weights(4) = 0.109951743655322;
    weights(5) = 0.109951743655322;

    points(0,0)=0.108103018168070;
    points(0,1)=0.445948490915965;
    points(0,2)=0.445948490915965;
    
    points(1,0)=0.445948490915965;
    points(1,1)=0.445948490915965;
    points(1,2)=0.108103018168070;

    points(2,0)=0.445948490915965;
    points(2,1)=0.108103018168070;
    points(2,2)=0.445948490915965;
  
    points(3,0)=0.816847572980459;
    points(3,1)=0.091576213509771;
    points(3,2)=0.091576213509771;
  
    points(4,0)=0.091576213509771;
    points(4,1)=0.091576213509771;
    points(4,2)=0.816847572980459;
  
    points(5,0)=0.091576213509771;
    points(5,1)=0.816847572980459;
    points(5,2)=0.091576213509771;
  }
  else if (noQpts==12){

    weights(0) = 0.116786275726379; 
    weights(1) = 0.116786275726379; 
    weights(2) = 0.116786275726379; 
    weights(3) = 0.050844906370207; 
    weights(4) = 0.050844906370207; 
    weights(5) = 0.050844906370207; 
    weights(6) = 0.082851075618374; 
    weights(7) = 0.082851075618374; 
    weights(8) = 0.082851075618374; 
    weights(9) = 0.082851075618374; 
    weights(10) = 0.082851075618374;
    weights(11) = 0.082851075618374;
    
    points(0,0)=0.501426509658179;  
    points(0,1)=0.249286745170910;  
    points(0,2)=0.249286745170910;
    
    points(1,0)=0.249286745170910;
    points(1,1)=0.249286745170910;
    points(1,2)=0.501426509658179;

    points(2,0)=0.249286745170910; 
    points(2,1)=0.501426509658179;
    points(2,2)=0.249286745170910;

    points(3,0)=0.873821971016996;
    points(3,1)=0.063089014491502;
    points(3,2)=0.063089014491502;
    
    points(4,0)=0.063089014491502;
    points(4,1)=0.063089014491502;
    points(4,2)=0.873821971016996;

    points(5,0)=0.063089014491502;
    points(5,1)=0.873821971016996;
    points(5,2)=0.063089014491502;
    
    points(6,0)=0.053145049844817;
    points(6,1)=0.310352451033784;
    points(6,2)=0.636502499121399;

    points(7,0)=0.310352451033784;
    points(7,1)=0.636502499121399;
    points(7,2)=0.053145049844817;

    points(8,0)=0.636502499121399;
    points(8,1)=0.053145049844817;
    points(8,2)=0.310352451033784;
    
    points(9,0)=0.636502499121399;
    points(9,1)=0.310352451033784;
    points(9,2)=0.053145049844817;

    points(10,0)=0.053145049844817;
    points(10,1)=0.636502499121399;
    points(10,2)=0.310352451033784;
    
    points(11,0)=0.310352451033784;
    points(11,1)=0.053145049844817;
    points(11,2)=0.636502499121399;
  }
}


void WavesPoissonEqOp::getTetraQuadrature2(const int noQpts,
						 MV_ColMat < real > &points,
						 MV_Vector < real > &weights)
{

  if (noQpts == 1) {
    points = 0.25;
    weights = 1.0 / 6.0;
  }
  else if (noQpts == 4) {


    const real g1 = 0.05 * (5.0 - sqrt(5.0));
    const real g2 = 0.05 * (5.0 + 3.0 * sqrt(5.0));


    points(0, 0) = g2;
    points(0, 1) = g1;
    points(0, 2) = g1;
    points(0, 3) = g1;
    points(1, 0) = g1;
    points(1, 1) = g2;
    points(1, 2) = g1;
    points(1, 3) = g1;
    points(2, 0) = g1;
    points(2, 1) = g1;
    points(2, 2) = g2;
    points(2, 3) = g1;
    points(3, 0) = g1;
    points(3, 1) = g1;
    points(3, 2) = g1;
    points(3, 3) = g2;

    weights = 1.0 / 24.0;


  }
  else if (noQpts == 11) {

    int noquadpnts = 11;
    points.newsize(noquadpnts, nsd+1);
    weights.newsize(noquadpnts);

    long double w1 = -0.0131555555555555555555555555555555;
    long double w2 = 0.00762222222222222222222222222222222;
    long double w3 = 0.0248888888888888888888888888888888;

    long double p1 = 0.25;
    long double p2 = 0.0714285714285714285714285714285714;
    long double p3 = 0.399403576166799204996102147461640;
    long double p4 = 0.100596423833200795003897852538359;
    long double p5 = 1.0 - 3.0 * p2;

    int c = -1;
    points(++c, 0) = p1;
    points(c, 1) = p1;
    points(c, 2) = p1;
    points(c, 3) = p1;

    points(++c, 0) = p5;
    points(c, 1) = p2;
    points(c, 2) = p2;
    points(c, 3) = p2;
    points(++c, 0) = p2;
    points(c, 1) = p5;
    points(c, 2) = p2;
    points(c, 3) = p2;
    points(++c, 0) = p2;
    points(c, 1) = p2;
    points(c, 2) = p5;
    points(c, 3) = p2;
    points(++c, 0) = p2;
    points(c, 1) = p2;
    points(c, 2) = p2;
    points(c, 3) = p5;

    points(++c, 0) = p3;
    points(c, 1) = p3;
    points(c, 2) = p4;
    points(c, 3) = p4;
    points(++c, 0) = p3;
    points(c, 1) = p4;
    points(c, 2) = p3;
    points(c, 3) = p4;
    points(++c, 0) = p3;
    points(c, 1) = p4;
    points(c, 2) = p4;
    points(c, 3) = p3;
    points(++c, 0) = p4;
    points(c, 1) = p3;
    points(c, 2) = p3;
    points(c, 3) = p4;
    points(++c, 0) = p4;
    points(c, 1) = p3;
    points(c, 2) = p4;
    points(c, 3) = p3;
    points(++c, 0) = p4;
    points(c, 1) = p4;
    points(c, 2) = p3;
    points(c, 3) = p3;

    weights(0) = w1;
    weights(1) = weights(2) = weights(3) = weights(4) = w2;
    weights(5) = weights(6) = weights(7) = w3;
    weights(8) = weights(9) = weights(10) = w3;
  }
  else if (noQpts == 24) {	// all weights positive and the points in the interior
    int noquadpnts = 24;
    points.newsize(noquadpnts, nsd+1);
    weights.newsize(noquadpnts);

    long double w1 = 0.00665379170969458201661510459291332;
    long double w2 = 0.00167953517588677382466887290765614;
    long double w3 = 0.00922619692394245368252554630895433;
    long double w4 = 0.00803571428571428571428571428571428;

    long double p1 = 0.214602871259152029288839219386284;
    long double p2 = 0.0406739585346113531155794489564100;
    long double p3 = 0.322337890142275510343994470762492;
    long double p4 = 0.0636610018750175252992355276057269;
    long double p5 = 0.269672331458315808034097805727606;
    long double p6 = 1.0 - 3.0 * p1;
    long double p7 = 1.0 - 3.0 * p2;
    long double p8 = 1.0 - 3.0 * p3;
    long double p9 = 1.0 - (2.0 * p4 + p5);

    int c = -1;
    points(++c, 0) = p1;
    points(c, 1) = p1;
    points(c, 2) = p1;
    points(c, 3) = p6;
    points(++c, 0) = p6;
    points(c, 1) = p1;
    points(c, 2) = p1;
    points(c, 3) = p1;
    points(++c, 0) = p1;
    points(c, 1) = p6;
    points(c, 2) = p1;
    points(c, 3) = p1;
    points(++c, 0) = p1;
    points(c, 1) = p1;
    points(c, 2) = p6;
    points(c, 3) = p1;

    points(++c, 0) = p2;
    points(c, 1) = p2;
    points(c, 2) = p2;
    points(c, 3) = p7;
    points(++c, 0) = p7;
    points(c, 1) = p2;
    points(c, 2) = p2;
    points(c, 3) = p2;
    points(++c, 0) = p2;
    points(c, 1) = p7;
    points(c, 2) = p2;
    points(c, 3) = p2;
    points(++c, 0) = p2;
    points(c, 1) = p2;
    points(c, 2) = p7;
    points(c, 3) = p2;

    points(++c, 0) = p3;
    points(c, 1) = p3;
    points(c, 2) = p3;
    points(c, 3) = p8;
    points(++c, 0) = p8;
    points(c, 1) = p3;
    points(c, 2) = p3;
    points(c, 3) = p3;
    points(++c, 0) = p3;
    points(c, 1) = p8;
    points(c, 2) = p3;
    points(c, 3) = p3;
    points(++c, 0) = p3;
    points(c, 1) = p3;
    points(c, 2) = p8;
    points(c, 3) = p3;

    points(++c, 0) = p5;
    points(c, 1) = p4;
    points(c, 2) = p9;
    points(c, 3) = p4;

    points(++c, 0) = p5;
    points(c, 1) = p4;
    points(c, 2) = p4;
    points(c, 3) = p9;
    points(++c, 0) = p4;
    points(c, 1) = p5;
    points(c, 2) = p9;
    points(c, 3) = p4;
    points(++c, 0) = p4;
    points(c, 1) = p5;
    points(c, 2) = p4;
    points(c, 3) = p9;
    points(++c, 0) = p4;
    points(c, 1) = p4;
    points(c, 2) = p5;
    points(c, 3) = p9;
    points(++c, 0) = p5;
    points(c, 1) = p9;
    points(c, 2) = p4;
    points(c, 3) = p4;

    points(++c, 0) = p9;
    points(c, 1) = p4;
    points(c, 2) = p5;
    points(c, 3) = p4;
    points(++c, 0) = p9;
    points(c, 1) = p4;
    points(c, 2) = p4;
    points(c, 3) = p5;
    points(++c, 0) = p4;
    points(c, 1) = p9;
    points(c, 2) = p5;
    points(c, 3) = p4;
    points(++c, 0) = p4;
    points(c, 1) = p9;
    points(c, 2) = p4;
    points(c, 3) = p5;
    points(++c, 0) = p4;
    points(c, 1) = p4;
    points(c, 2) = p9;
    points(c, 3) = p5;
    points(++c, 0) = p9;
    points(c, 1) = p5;
    points(c, 2) = p4;
    points(c, 3) = p4;

    weights(0) = weights(1) = weights(2) = weights(3) = w1;
    weights(4) = weights(5) = weights(6) = weights(7) = w2;
    weights(8) = weights(9) = weights(10) = weights(11) = w3;
    int i;
    for (i = 12; i < 24; i++)
      weights(i) = w4;
  }
  else
    cout << "No quadrature points = " << noQpts << " not available!\n" <<
      flush;
}




void WavesPoissonEqOp::constructGrad_x(Mat& Gradient_x)
{
  int i,n,el;
  int ierr;
  nsd = gg.getNoSpaceDim();
  int nno = gg.getNoNodes();
  int nel = gg.getNoElms();

  
  PetscScalar zero=0.0;
  // cout<<"before vecset"<<endl;
  
 
  WavesNeighborFE& gridNeighbors = gg.getNeighbor ();
  gridNeighbors.init (gg,false,true,false);
  MV_Vector<int> number_of_nbs(nno);
  int count=0;
  for(i=0;i<nno;i++) // Node connecticivity
    number_of_nbs(count++) = (gridNeighbors.couplingsIrow (i+1) -
			      gridNeighbors.couplingsIrow (i));
  gridNeighbors.remove();

  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nno,nno,0, &number_of_nbs(0), &Gradient_x);
  CHKERRA(ierr);
  
 // cout << "Assembling matrix 01 in constructGrad_x" << endl;
  int nobasenodes = nsd+1;// for triangles and tetras
  Array1dReal resultElMat(nobasenodes*nobasenodes);  
  Array1dReal resultElMat_B(nobasenodes*nobasenodes); 
  Array1dInt rowAndColIdx(nobasenodes); 
  
//  cout<<"loop"<<endl;
  if(nsd==2){
    FET3n2D F(&gg);
    for( el = 0; el < nel; el++ ){
      if(gg.getElementType(el)==ELMTRI1) {
	F.refill(el);
	real volume = F.area();

	rowAndColIdx(0) = F.n1();
	rowAndColIdx(1) = F.n2();
	rowAndColIdx(2) = F.n3();
	/*
	resultElMat_B(0) = volume*F.grad_x_a11();
	resultElMat_B(4) = volume*F.grad_x_a22();
	resultElMat_B(8) = volume*F.grad_x_a33();
	*/
        
  	volume *= 0.3333333333333333;

        //volume /= -F.det();

	resultElMat_B(0) = volume*F.dN1x();
	resultElMat_B(4) = volume*F.dN2x();
	resultElMat_B(8) = volume*F.dN3x();
	
	/*
	resultElMat_B(2) = resultElMat_B(1) = resultElMat_B(0);
	resultElMat_B(3) = resultElMat_B(5) = resultElMat_B(4);
	resultElMat_B(6) = resultElMat_B(7) = resultElMat_B(8);
	*/
	resultElMat_B(2) = resultElMat_B(1) = 0.0;
	resultElMat_B(3) = resultElMat_B(5) = 0.0;
	resultElMat_B(6) = resultElMat_B(7) = 0.0;
	
	ierr = MatSetValues(Gradient_x,nobasenodes,&rowAndColIdx(0),
			    nobasenodes,&rowAndColIdx(0),&resultElMat_B(0),
			    ADD_VALUES);  CHKERRA(ierr);


	//	cout<<" Assembl. matrix for grad x, element"<<el<<" area is "<<volume<<endl;
	//  cout<<resultElMat_B(0)<<" ; "<<resultElMat_B(1)<<"  ; "<<resultElMat_B(2)<<endl;
	//	cout<<resultElMat_B(3)<<" ; "<<resultElMat_B(4)<<"  ; "<<resultElMat_B(5)<<endl;
	//cout<<resultElMat_B(6)<<" ; "<<resultElMat_B(7)<<"  ; "<<resultElMat_B(8)<<endl;



      }
    }
  } else if(nsd==3){
    FET4n3D F(&gg);
    for( el = 0; el < nel; el++ ){
      if(gg.getElementType(el)==ELMTET1) {
	F.refill(el);
	real volume = fabs(F.volume());

        if (volume < 0.0)
	  {
	    cout<<"negative volume !!! "<<volume<<" in element "<<el<<endl;
            exit(1);

	  }

	volume *= 0.25;

        volume /= -F.det();

	rowAndColIdx(0) = F.n1();
	rowAndColIdx(1) = F.n2();
	rowAndColIdx(2) = F.n3();
	rowAndColIdx(3) = F.n4();


	resultElMat_B(0) =  volume*F.dN1x();
	resultElMat_B(5) =  volume*F.dN2x();
	resultElMat_B(10) = volume*F.dN3x();
	resultElMat_B(15) = volume*F.dN4x();

  	resultElMat_B(1) = resultElMat_B(2) = resultElMat_B(3) = 0.0;
	resultElMat_B(4) = resultElMat_B(6) = resultElMat_B(7) = 0.0;
	resultElMat_B(8) = resultElMat_B(9) = resultElMat_B(11) = 0.0;
	resultElMat_B(12) = resultElMat_B(13) = resultElMat_B(14) = 0.0;


	ierr = MatSetValues(Gradient_x,nobasenodes,&rowAndColIdx(0),
			    nobasenodes,&rowAndColIdx(0),&resultElMat_B(0),
			    ADD_VALUES);  
	CHKERRA(ierr);

      }
    }
  }

  ierr = MatAssemblyBegin(Gradient_x,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);
  ierr = MatAssemblyEnd(Gradient_x,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);

 // ierr=dirichletBCxxx(A); CHKERRA(ierr);
  //ierr=dirichletBCxxx(B); CHKERRA(ierr);
}
  
void WavesPoissonEqOp::constructGrad_y(Mat& Gradient_y)
{
  int i,n,el;
  int ierr;
  nsd = gg.getNoSpaceDim();
  int nno = gg.getNoNodes();
  int nel = gg.getNoElms();

  
  PetscScalar zero=0.0;
  // cout<<"before vecset"<<endl;
  
  WavesNeighborFE& gridNeighbors = gg.getNeighbor ();
  gridNeighbors.init (gg,false,true,false);
  MV_Vector<int> number_of_nbs(nno);
  int count=0;
  for(i=0;i<nno;i++) // Node connecticivity
    number_of_nbs(count++) = (gridNeighbors.couplingsIrow (i+1) -
			      gridNeighbors.couplingsIrow (i));
  gridNeighbors.remove();
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nno,nno,0, &number_of_nbs(0), & Gradient_y);
  CHKERRA(ierr);
  
 // cout << "Assembling matrix 01 in constructGrad_y" << endl;
  int nobasenodes = nsd+1;// for triangles and tetras
  Array1dReal resultElMat(nobasenodes*nobasenodes);  
  Array1dReal resultElMat_B(nobasenodes*nobasenodes); 
  Array1dInt rowAndColIdx(nobasenodes); 
  
 // cout<<"loop"<<endl;
  if(nsd==2){
    FET3n2D F(&gg);
    for( el = 0; el < nel; el++ ){
      if(gg.getElementType(el)==ELMTRI1) {
	F.refill(el);
	real volume = F.area();

	rowAndColIdx(0) = F.n1();
	rowAndColIdx(1) = F.n2();
	rowAndColIdx(2) = F.n3();
	/*
	resultElMat_B(0) = volume*F.grad_x_a11();
	resultElMat_B(4) = volume*F.grad_x_a22();
	resultElMat_B(8) = volume*F.grad_x_a33();
	*/
        
  	volume *= 0.3333333333333333;

	//        volume /= -F.det();

	resultElMat_B(0) = volume*F.dN1y();
	resultElMat_B(4) = volume*F.dN2y();
	resultElMat_B(8) = volume*F.dN3y();

	/*
	
	resultElMat_B(2) = resultElMat_B(1) = resultElMat_B(0);
	resultElMat_B(3) = resultElMat_B(5) = resultElMat_B(4);
	resultElMat_B(6) = resultElMat_B(7) = resultElMat_B(8);
	*/
	

	resultElMat_B(2) = resultElMat_B(1) = 0.0;
	resultElMat_B(3) = resultElMat_B(5) = 0.0;
	resultElMat_B(6) = resultElMat_B(7) = 0.0;
	
	ierr = MatSetValues(Gradient_y,nobasenodes,&rowAndColIdx(0),
			    nobasenodes,&rowAndColIdx(0),&resultElMat_B(0),
			    ADD_VALUES);  CHKERRA(ierr);


	//	cout<<" Assembl. matrix for grad y, element"<<el<<" area is "<<volume<<endl;
  //      cout<<resultElMat_B(0)<<" ; "<<resultElMat_B(1)<<"  ; "<<resultElMat_B(2)<<endl;
//	cout<<resultElMat_B(3)<<" ; "<<resultElMat_B(4)<<"  ; "<<resultElMat_B(5)<<endl;
//	cout<<resultElMat_B(6)<<" ; "<<resultElMat_B(7)<<"  ; "<<resultElMat_B(8)<<endl;



      }
    }
  } else if(nsd==3){
    FET4n3D F(&gg);
    for( el = 0; el < nel; el++ ){
      if(gg.getElementType(el)==ELMTET1) {
	F.refill(el);
	real volume = fabs(F.volume());

        if (volume < 0.0)
	  {
	    cout<<"negative volume !!! "<<volume<<" in element "<<el<<endl;
            exit(1);

	  }

	volume *= 0.25;

	//  volume /= -F.det();

	rowAndColIdx(0) = F.n1();
	rowAndColIdx(1) = F.n2();
	rowAndColIdx(2) = F.n3();
	rowAndColIdx(3) = F.n4();


	resultElMat_B(0) =  volume*F.dN1y();
	resultElMat_B(5) =  volume*F.dN2y();
	resultElMat_B(10) = volume*F.dN3y();
	resultElMat_B(15) = volume*F.dN4y();

  	resultElMat_B(1) = resultElMat_B(2) = resultElMat_B(3) = 0.0;
	resultElMat_B(4) = resultElMat_B(6) = resultElMat_B(7) = 0.0;
	resultElMat_B(8) = resultElMat_B(9) = resultElMat_B(11) = 0.0;
	resultElMat_B(12) = resultElMat_B(13) = resultElMat_B(14) = 0.0;


	ierr = MatSetValues(Gradient_y,nobasenodes,&rowAndColIdx(0),
			    nobasenodes,&rowAndColIdx(0),&resultElMat_B(0),
			    ADD_VALUES);  
	CHKERRA(ierr);

      }
    }
  }

  ierr = MatAssemblyBegin(Gradient_y,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);
  ierr = MatAssemblyEnd(Gradient_y,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);

 // ierr=dirichletBCxxx(A); CHKERRA(ierr);
  //ierr=dirichletBCxxx(B); CHKERRA(ierr);
}
  


void WavesPoissonEqOp::constructGrad_z(Mat&  Gradient_z)
{
  int i,n,el;
  int ierr;
  nsd = gg.getNoSpaceDim();
  int nno = gg.getNoNodes();
  int nel = gg.getNoElms();
  
  PetscScalar zero=0.0;
  // cout<<"before vecset"<<endl;
  
  WavesNeighborFE& gridNeighbors = gg.getNeighbor ();
  gridNeighbors.init (gg,false,true,false);
  MV_Vector<int> number_of_nbs(nno);
  int count=0;
  for(i=0;i<nno;i++) // Node connecticivity
    number_of_nbs(count++) = (gridNeighbors.couplingsIrow (i+1) -
			      gridNeighbors.couplingsIrow (i));
  gridNeighbors.remove();
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nno,nno,0, &number_of_nbs(0), &Gradient_z);
  CHKERRA(ierr);
  
//  cout << "Assembling matrix 01 in constructGrad_z" << endl;
  int nobasenodes = nsd+1;// for triangles and tetras
  Array1dReal resultElMat(nobasenodes*nobasenodes);  
  Array1dReal resultElMat_Z(nobasenodes*nobasenodes); 
  Array1dInt rowAndColIdx(nobasenodes); 
  
  //  cout<<"loop"<<endl;
 if(nsd==3){
    FET4n3D F(&gg);
    for( el = 0; el < nel; el++ ){
      if(gg.getElementType(el)==ELMTET1) {
	F.refill(el);
	real volume = fabs(F.volume());

        if (volume < 0.0)
	  {
	    cout<<"negative volume !!! "<<volume<<" in element "<<el<<endl;
            exit(1);

	  }

	volume *= 0.25;

	// volume /= -F.det();

	rowAndColIdx(0) = F.n1();
	rowAndColIdx(1) = F.n2();
	rowAndColIdx(2) = F.n3();
	rowAndColIdx(3) = F.n4();


	resultElMat_Z(0) =  volume*F.dN1z();
	resultElMat_Z(5) =  volume*F.dN2z();
	resultElMat_Z(10) = volume*F.dN3z();
	resultElMat_Z(15) = volume*F.dN4z();

  	resultElMat_Z(1) = resultElMat_Z(2) = resultElMat_Z(3) = 0.0;
	resultElMat_Z(4) = resultElMat_Z(6) = resultElMat_Z(7) = 0.0;
	resultElMat_Z(8) = resultElMat_Z(9) = resultElMat_Z(11) = 0.0;
	resultElMat_Z(12) = resultElMat_Z(13) = resultElMat_Z(14) = 0.0;

	//cout<<"before matsetvalues"<<endl;

	ierr = MatSetValues(Gradient_z,nobasenodes,&rowAndColIdx(0),
			    nobasenodes,&rowAndColIdx(0),&resultElMat_Z(0),
			    ADD_VALUES);   
	CHKERRA(ierr);
	//	cout<<"after matsetvalues"<<endl;
      }
    }
  }

  ierr = MatAssemblyBegin(Gradient_z,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);
  ierr = MatAssemblyEnd(Gradient_z,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);

 // ierr=dirichletBCxxx(A); CHKERRA(ierr);
  //ierr=dirichletBCxxx(B); CHKERRA(ierr);
}



//=============  find gradient of function q for glob conv alg
//=================================================================================================

void WavesPoissonEqOp::GradFEM(Mat& Gradient, MV_Vector<double>& q, Vec& u_2 )
{
  
  int nno = gg.getNoNodes();
  PetscScalar zero = 0.0;
  Vec vec_q;


  VecCreate(PETSC_COMM_SELF, &vec_q);
  ierr = VecSetSizes(vec_q, PETSC_DECIDE, nno); CHKERRA(ierr); 
  ierr = VecSetType(vec_q, VECMPI);
  CHKERRA(ierr);

  VecCreate(PETSC_COMM_SELF, &u_2);
  ierr = VecSetSizes(u_2, PETSC_DECIDE, nno); CHKERRA(ierr); 
  ierr = VecSetType(u_2, VECMPI);
  CHKERRA(ierr);

  ierr = VecSet(vec_q, zero);
  CHKERRA(ierr);
  ierr = VecSet(u_2, zero);
 CHKERRA(ierr);

  PetscScalar *qvalues;
  ierr = VecGetArray(vec_q, &qvalues);
  CHKERRA(ierr);

  for (int i = 0; i < nno; i++) 
    qvalues[i] = q(i);
  
  //	cout<<"  bvalues("<<i<<")= "<<bvalues[i]<<" boundindex("<<i<<") = "<<boundindex[i]<<endl;
  ierr = VecRestoreArray(vec_q, &qvalues);
  CHKERRA(ierr);

//  cout << "before matmult" << endl;
  ierr = MatMult(Gradient, vec_q, u_2);   
  //  ierr = MatMult(A, vec_q, u_1);   
  // ierr = MatMultAdd(B, vec_q, u_1, u_2);
  CHKERRA(ierr);

}   
  
//====================================================================================
void WavesPoissonEqOp::Convert( Vec& u_2, double* solfdm)
{

  PetscScalar *u_2_values;
  ierr = VecGetArray(u_2, &u_2_values);

  CHKERRA(ierr);
  
  for (int i = 0; i < gg.getNoNodes(); i++) 
    {
        solfdm[i] = u_2_values[i];
//	cout<<"computed FEM-value"<<solfdm[i]<<endl;
    }
  ierr = VecRestoreArray(u_2, &u_2_values);
  CHKERRA(ierr);
  
}

//====================================================================================

void WavesPoissonEqOp::constructStiffnessMatrixOpt( Mat& A)
{
  int i,n,el;
  int ierr;
  nsd = gg.getNoSpaceDim();
  int nno = gg.getNoNodes();
  int nel = gg.getNoElms();

  
  PetscScalar zero=0.0;
 // cout<<"before vecset"<<endl;
  
  /* commented 18.10.2000 becouse of error in NeighborFE.init for
     case ellipse with hole
  */
  WavesNeighborFE& gridNeighbors = gg.getNeighbor ();
  gridNeighbors.init (gg,false,true,false);
  MV_Vector<int> number_of_nbs(nno);
  int count=0;
  for(i=0;i<nno;i++) // Node connecticivity
    number_of_nbs(count++) = (gridNeighbors.couplingsIrow (i+1) -
			      gridNeighbors.couplingsIrow (i));
  gridNeighbors.remove();
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nno,nno,0, &number_of_nbs(0), &A);
  CHKERRA(ierr);

  /*  commented 18.10.2000*/

  /*
  //added 18.10.2000
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nno,nno,0,PETSC_NULL, &A_);
  CHKERRA(ierr);
  */

 // cout << "Assembling matrix 01." << endl;
  int nobasenodes = nsd+1;// for triangles and tetras
  Array1dReal resultElMat(nobasenodes*nobasenodes); 
  Array1dInt rowAndColIdx(nobasenodes); 

//  cout<<"loop"<<endl;
  if(nsd==2){
    FET3n2D F(&gg);
    for( el = 0; el < nel; el++ ){
      if(gg.getElementType(el)==ELMTRI1) {
	F.refill(el);
	real volume = F.area();
	rowAndColIdx(0) = F.n1();
	rowAndColIdx(1) = F.n2();
	rowAndColIdx(2) = F.n3();
	resultElMat(0) = volume*F.a11();
	resultElMat(4) = volume*F.a22();
	resultElMat(8) = volume*F.a33();
	resultElMat(3) = resultElMat(1) = volume*F.a12();
	resultElMat(6) = resultElMat(2) = volume*F.a13();
	resultElMat(5) = resultElMat(7) = volume*F.a23();
	ierr = MatSetValues(A,nobasenodes,&rowAndColIdx(0),
			    nobasenodes,&rowAndColIdx(0),&resultElMat(0),
			    ADD_VALUES);  CHKERRA(ierr);

	
      }
    }
  } else if(nsd==3){
    FET4n3D F(&gg);
    for( el = 0; el < nel; el++ ){
      if(gg.getElementType(el)==ELMTET1) {
	F.refill(el);
	real volume = fabs(F.volume());

        if (volume < 0.0)
	  {
	    cout<<"negative volume !!! "<<volume<<" in element "<<el<<endl;
            exit(1);

	  }
	rowAndColIdx(0) = F.n1();
	rowAndColIdx(1) = F.n2();
	rowAndColIdx(2) = F.n3();
	rowAndColIdx(3) = F.n4();
	resultElMat(0) = volume*F.a11();
	resultElMat(5) = volume*F.a22();
	resultElMat(10) = volume*F.a33();
	resultElMat(15) = volume*F.a44();
	resultElMat(4) = resultElMat(1) = volume*F.a12();
	resultElMat(8) = resultElMat(2) = volume*F.a13();
	resultElMat(12) = resultElMat(3) = volume*F.a14();
	resultElMat(9) = resultElMat(6) = volume*F.a23();
	resultElMat(13) = resultElMat(7) = volume*F.a24();
	resultElMat(14) = resultElMat(11) = volume*F.a34();
	ierr = MatSetValues(A,nobasenodes,&rowAndColIdx(0),
			    nobasenodes,&rowAndColIdx(0),&resultElMat(0),
			    ADD_VALUES);  
	CHKERRA(ierr);

      }
    }
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);
  // cout<<"before dirichletBCxxx( A) "<<endl;

  //    ierr=dirichletBCxxx(A); CHKERRA(ierr);


  //    Scalar eps = 0.00001;
  //  ierr = MatShift(&eps,A); CHKERRA(ierr);
}
  



void WavesPoissonEqOp::constructStiffnessMatrixOpt( Mat& A,  MV_Vector<double>& grad_q_x,  
							  MV_Vector<double>& grad_q_y, 
							  MV_Vector<double>& grad_q_z, double eps_)
{
  int i,n,el;
  int ierr;
  nsd = gg.getNoSpaceDim();
  int nno = gg.getNoNodes();
  int nel = gg.getNoElms();

  
  PetscScalar zero=0.0;
 // cout<<"before vecset"<<endl;
  //=========================================================================
  WavesNeighborFE& gridNeighbors = gg.getNeighbor ();
  gridNeighbors.init (gg,false,true,false);
  MV_Vector<int> number_of_nbs(nno);
  int count=0;
  for(i=0;i<nno;i++) // Node connecticivity
    number_of_nbs(count++) = (gridNeighbors.couplingsIrow (i+1) -
			      gridNeighbors.couplingsIrow (i));
  gridNeighbors.remove();
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nno,nno,0, &number_of_nbs(0), &A);
  CHKERRA(ierr);

  /*
  // if number_of_nbs  does't work, this is more slower variant:
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nno,nno,0,PETSC_NULL, &A);
  CHKERRA(ierr);
  */
  //===================================================================
 // cout << "Assembling matrix 01." << endl;
  int nobasenodes = nsd+1;// for triangles and tetras
  Array1dReal resultElMat(nobasenodes*nobasenodes); 
  Array1dInt rowAndColIdx(nobasenodes); 

//  cout<<"loop"<<endl;
  if(nsd==2){
    FET3n2D F(&gg);
    for( el = 0; el < nel; el++ ){
      if(gg.getElementType(el)==ELMTRI1) {
	F.refill(el);
	real volume = F.area();
	rowAndColIdx(0) = F.n1();
	rowAndColIdx(1) = F.n2();
	rowAndColIdx(2) = F.n3();
	
	resultElMat(0) = volume*(F.a11() - 
				 (F.dN1x()*grad_q_x(F.n1()) 
				  + F.dN1y()*grad_q_y(F.n1())));
        resultElMat(1) = volume*(F.a12() -
				 (F.dN1x()*grad_q_x(F.n1()) 
				  + F.dN1y()*grad_q_y(F.n1())));
	resultElMat(2) = volume*(F.a13() - 
				 (F.dN1x()*grad_q_x(F.n1())
				  + F.dN1y()*grad_q_y(F.n1())));

	resultElMat(3) = volume*(F.a12()  -  
				 (F.dN2x()*grad_q_x(F.n2())
				  + F.dN2y()*grad_q_y(F.n2())));
	resultElMat(4) = volume*(F.a22()  -  
				 (F.dN2x()*grad_q_x(F.n2()) 
				  + F.dN2y()*grad_q_y(F.n2())));
	resultElMat(5) =  volume*(F.a23() -  
				  (F.dN2x()*grad_q_x(F.n2()) 
				   + F.dN2y()*grad_q_y(F.n2())));

	resultElMat(6) =  volume*(F.a13() - 
				  (F.dN3x()*grad_q_x(F.n3()) 
				   + F.dN3y()*grad_q_y(F.n3())));
	resultElMat(7) = volume*(F.a23() - 
				 (F.dN3x()*grad_q_x(F.n3()) 
				  + F.dN3y()*grad_q_y(F.n3())));
	resultElMat(8) = volume*(F.a33() - 
				 (F.dN3x()*grad_q_x(F.n3()) 
				  + F.dN3y()*grad_q_y(F.n3())));

	ierr = MatSetValues(A,nobasenodes,&rowAndColIdx(0),
			    nobasenodes,&rowAndColIdx(0),&resultElMat(0),
			    ADD_VALUES);  CHKERRA(ierr);

	
      }
    }
  } else if(nsd==3){
    FET4n3D F(&gg);
    for( el = 0; el < nel; el++ ){
      if(gg.getElementType(el)==ELMTET1) {
	F.refill(el);
	real volume = fabs(F.volume());

        if (volume < 0.0)
	  {
	    cout<<"negative volume !!! "<<volume<<" in element "<<el<<endl;
            exit(1);

	  }
	rowAndColIdx(0) = F.n1();
	rowAndColIdx(1) = F.n2();
	rowAndColIdx(2) = F.n3();
	rowAndColIdx(3) = F.n4();
	// assembling of matrice (grad_f_i,grad_f_j) - (grad f_i*grad_q_i, f_j)
	// note that f_j(n_i) = delta_ij, and actually previous formula reduses to
	//  (grad_f_i,grad_f_j) - (grad f_i*grad_q_i, 1)



	resultElMat(0) = volume*(F.a11() - (F.dN1x()*grad_q_x(F.n1()) +
					    F.dN1y()*grad_q_y(F.n1()) +
					    F.dN1z()*grad_q_z(F.n1())));

	
        resultElMat(1) = volume*(F.a12() - (F.dN1x()*grad_q_x(F.n1())  + 
					    F.dN1y()*grad_q_y(F.n1())  +
					    F.dN1z()*grad_q_z(F.n1())));
							      
				 
	resultElMat(2) = volume*(F.a13() - (F.dN1x()*grad_q_x(F.n1())  +
					    F.dN1y()*grad_q_y(F.n1()) +
					    F.dN1z()*grad_q_z(F.n1())));


       resultElMat(3) =  volume*(F.a14() - (F.dN1x()*grad_q_x(F.n1())  +
					    F.dN1y()*grad_q_y(F.n1()) +
					    F.dN1z()*grad_q_z(F.n1())));



	resultElMat(4) = volume*(F.a12() -  (F.dN2x()*grad_q_x(F.n2())  +
					      F.dN2y()*grad_q_y(F.n2()) +
					      F.dN2z()*grad_q_z(F.n2())));


	resultElMat(5) = volume*(F.a22() -  (F.dN2x()*grad_q_x(F.n2()) + 
					      F.dN2y()*grad_q_y(F.n2()) +
					      F.dN2z()*grad_q_z(F.n2())));


	resultElMat(6) =  volume*(F.a23() -  (F.dN2x()*grad_q_x(F.n2()) +
					      F.dN2y()*grad_q_y(F.n2())+
					      F.dN2z()*grad_q_z(F.n2())));



	resultElMat(7) = volume*(F.a24()  -  (F.dN2x()*grad_q_x(F.n2()) +
					      F.dN2y()*grad_q_y(F.n2())+
					      F.dN2z()*grad_q_z(F.n2()))); 




	resultElMat(8) =  volume*(F.a13() - (F.dN3x()*grad_q_x(F.n3()) +
					     F.dN3y()*grad_q_y(F.n3()) + 
					     F.dN3z()*grad_q_z(F.n3())));


	resultElMat(9) = volume*(F.a23() - (F.dN3x()*grad_q_x(F.n3()) +
					     F.dN3y()*grad_q_y(F.n3()) + 
					     F.dN3z()*grad_q_z(F.n3()))); 


	resultElMat(10) = volume*(F.a33() -  (F.dN3x()*grad_q_x(F.n3()) +
					     F.dN3y()*grad_q_y(F.n3()) + 
					     F.dN3z()*grad_q_z(F.n3())));


	resultElMat(11)= volume*(F.a34() - (F.dN3x()*grad_q_x(F.n3()) +
					     F.dN3y()*grad_q_y(F.n3()) + 
					     F.dN3z()*grad_q_z(F.n3())));



       resultElMat(12) = volume*(F.a14()  -  (F.dN4x()*grad_q_x(F.n4()) +
					      F.dN4y()*grad_q_y(F.n4()) + 
					      F.dN4z()*grad_q_z(F.n4())));

	resultElMat(13) = volume*(F.a24() -  (F.dN4x()*grad_q_x(F.n4()) +
					      F.dN4y()*grad_q_y(F.n4()) + 
					      F.dN4z()*grad_q_z(F.n4())));

	resultElMat(14) = volume*(F.a34() -  (F.dN4x()*grad_q_x(F.n4()) +
					      F.dN4y()*grad_q_y(F.n4()) + 
					      F.dN4z()*grad_q_z(F.n4())));

	resultElMat(15) = volume*(F.a44()  -  (F.dN4x()*grad_q_x(F.n4()) +
					      F.dN4y()*grad_q_y(F.n4()) + 
					      F.dN4z()*grad_q_z(F.n4())));

	
	ierr = MatSetValues(A,nobasenodes,&rowAndColIdx(0),
			    nobasenodes,&rowAndColIdx(0),&resultElMat(0),
			    ADD_VALUES);  
	CHKERRA(ierr);

      }
    }
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);
  // cout<<"before dirichletBCxxx( A) "<<endl;

  //    ierr=dirichletBCxxx(A); CHKERRA(ierr);

  //  regularization parameter
             PetscScalar eps = eps_;
	     
	     //old petsc
//   ierr = MatShift(&eps,A); CHKERRA(ierr);
   
ierr = MatShift(A,eps); CHKERRA(ierr);
}
  

void WavesPoissonEqOp::InitPoissonFEM()
{
  
  int ierr;
  // Mat Ahelp;
  
  int noeq = 1;
  int nno = noeq*gg.getNoNodes();
  // int ierr,
  int i;
  
  // iteration matrix for system
  
    constructStiffnessMatrixOpt(A);

  // ierr = dropZeroElements(Ahelp,A,1e-13);
 
  // ierr = MatDestroy(Ahelp); CHKERRA(ierr);

   Poisson_FEM_initialized = 1;
}

void WavesPoissonEqOp::InitGlobConvFEM(MV_Vector<double>& grad_q_x,
					     MV_Vector<double>& grad_q_y, 
					     MV_Vector<double>& grad_q_z, double eps)
{
  
  int ierr;
  // Mat Ahelp;
  
  // int nsd = 2;
   
  int nno = gg.getNoNodes();
  // int ierr,
  int i;
  
  // iteration matrix for system
 // cout<<"nsd "<<nsd<<endl;
  
      constructStiffnessMatrixOpt(A, grad_q_x,  grad_q_y, grad_q_z, eps);
 
  // ierr = dropZeroElements(Ahelp,A,1e-13);
 
  // ierr = MatDestroy(Ahelp); CHKERRA(ierr);

   Poisson_FEM_initialized = 1;
}





void WavesPoissonEqOp::ApplyDirichletFEM(Vec& exact)
{
  // create bndNodes  on the boundary 
  MV_Vector < int >markNodes;

 // cout<<" inside applydirichletfem "<<endl;
  opt.findBoundaryNodes(gg, markNodes);
  MV_Vector<int> bndNodes;
  opt.makeSortedNodesIndex(markNodes, bndNodes, 1);
  int nbn_ = bndNodes.size();
//  cout<<"number of boundary nodes in PoissonEqOp::ApplyDirichletFEM: nbn = "<<nbn_<<endl;
  nbn = nbn_;

  // boundindex = (int *) PetscMalloc(nbn * sizeof(int));
  // bvalues = (Scalar *) PetscMalloc(nbn * sizeof(Scalar));

    bvalues_initialised = 1;


  boundindex = new int[nbn_];
  bvalues = new double[nbn_];

  PetscScalar *tmpvec;
  ierr = VecGetArray(exact, &tmpvec);
  CHKERRA(ierr);

  for (int i = 0; i < nbn_; i++) {

    boundindex[i] = bndNodes(i);
    bvalues[i] = tmpvec[boundindex[i]];

    //    cout<<"  bvalues("<<i<<") "<< tmpvec[ boundindex[i]]<<"nbn "<<nbn_<<endl;
  
}
  //  cout<<"before applyDirichlet"<<endl;

  ierr = VecRestoreArray(exact, &tmpvec);
  CHKERRA(ierr);
}



void WavesPoissonEqOp::ApplyDirichletFEM(MV_Vector<double>& exact)
{
  // create bndNodes  on the boundary 
  MV_Vector < int >markNodes;

 // cout<<" inside applydirichletfem "<<endl;
  opt.findBoundaryNodes(gg, markNodes);
  MV_Vector<int> bndNodes;
  opt.makeSortedNodesIndex(markNodes, bndNodes, 1);
  int nbn_ = bndNodes.size();
//  cout<<"nbn = "<<nbn_<<endl;
  nbn = nbn_;

  // boundindex = (int *) PetscMalloc(nbn * sizeof(int));
  // bvalues = (Scalar *) PetscMalloc(nbn * sizeof(Scalar));

    bvalues_initialised = 1;


  boundindex = new int[nbn_];
  bvalues = new double[nbn_];



  for (int i = 0; i < nbn_; i++) {

    boundindex[i] = bndNodes(i);

    bvalues[i] = exact(boundindex[i]);
  
//      cout<<"  bvalues("<<i<<") ="<<exact(boundindex[i])<<"nbn "<<nbn_<<endl;
  
}
  //  cout<<"before applyDirichlet"<<endl;


}


void  WavesPoissonEqOp::DirichletFEM(MV_Vector<double>& exsol, real* function)
{
  
  //  exsol - values at the boundary, are assigned to bvalues in next function
  ApplyDirichletFEM(exsol);
  
  for (int i = 0; i < nbn; i++) {
    function[boundindex[i]] = bvalues[i];
    //	 boundvalues[boundindex[i]] = -1.0;
 //    cout<<"  boundary values "<<function[boundindex[i]] <<endl;
    
  }
  
}



void  WavesPoissonEqOp::FEM_Boundary(MV_Vector<double>& exsol,
				     MV_Vector<double>&  func_obs)
{
  
  //  exsol - values at the boundary, are assigned to bvalues in next function
  ApplyDirichletFEM(exsol);
  
  // boundindex - is array with boundary points
  // size of func_obs is size of nbn

    for (int i = 0; i < nbn; i++) {
      func_obs(boundindex[i]) = bvalues[i];
      //	 boundvalues[boundindex[i]] = -1.0;
      // cout<<"  boundary values "<<function[boundindex[i]] <<endl;
  
    }
 
  }
 

void WavesPoissonEqOp::evalRHS(Vec& Fn_, MV_Vector<double>& rhs_) 
{
  int n,e,q;
  //  int  ierr;
  nsd=gg.getNoSpaceDim();;
  int nel = gg.getNoElms();
  int nno = gg.getNoNodes();
  int noBasisFcn = nsd+1;

   noQpts = nsd+1;
  // change number of quadrature points 
   //noQpts = 36;
  //==================================
//  cout<<" noQpts in evalRHS = "<<noQpts<<" nsd = "<<nsd<<endl;

  MV_ColMat<real> points(noQpts,nsd+1);
  MV_Vector<real> weights(noQpts);

  if(nsd==2)
    getTriangleQuadrature(noQpts, 
			  points, 
			  weights);
  else  if(nsd==3)
    getTetraQuadrature2(noQpts, 
		       points, 
		       weights);
  PetscScalar zero = 0.0;
  ierr = VecSet(Fn_, zero); CHKERRA(ierr);
  PetscScalar  *tmpvec;
  ierr = VecGetArray( Fn_, &tmpvec); CHKERRA(ierr);
 

  if(nsd==2) {
  
  FET3n2D F(&gg);
  for(e=0;e<nel;e++){
    F.refill(e);
    real area = F.area();
    for(q=0;q<noQpts;q++){
       
      real areaWeightRhs = area*weights(q);
     
      tmpvec [F.n1()] += points(q,0)*areaWeightRhs*rhs_(F.n1());
      tmpvec [F.n2()] += points(q,1)*areaWeightRhs*rhs_(F.n2());
      tmpvec [F.n3()] += points(q,2)*areaWeightRhs*rhs_(F.n3());

     

    }
  }
  } 
  else if(nsd==3) {
    //cout<<"nsd ="<<nsd<<"time is "<<t<<endl;
    FET4n3D F(&gg);

      for(e=0;e<nel;e++){
	F.refill(e);
	real area = F.det();
         
        if ( area < 0.0)
	  {
	    cout<<" area of element is negative !!!"<<endl;
            exit(1);
	    
	  }

	for(q=0;q<noQpts;q++){
	  
      
	  real areaWeightRhs = area*weights(q);
 
      tmpvec [F.n1()] += points(q,0)*areaWeightRhs*rhs_(F.n1());
      tmpvec [F.n2()] += points(q,1)*areaWeightRhs*rhs_(F.n2());
      tmpvec [F.n3()] += points(q,2)*areaWeightRhs*rhs_(F.n3());
      tmpvec [F.n4()] += points(q,3)*areaWeightRhs*rhs_(F.n4());
    }}}  
  
  
  // for (int i=0; i < gg.getNoNodes(); i++)
  // cout<<"  rhs in evalRHS "<<tmpvec[i]<<endl;

  ierr = VecRestoreArray( Fn_, &tmpvec ); CHKERRA(ierr); 

} 

// Here we assign values from exsol to the values of exsol only at the boundary
// of the domain. At all inner nodes of gg initvalues=0.

void WavesPoissonEqOp::PoissonEqInitBoundaryValues(MV_Vector<double>& exsol, 
						   double* initbvalues)
{
  PetscScalar zero = 0.0;
  int nno = gg.getNoNodes();

  Vec u0;
  VecCreate(PETSC_COMM_SELF, &u0);
  ierr = VecSetSizes(u0, PETSC_DECIDE, nno); CHKERRA(ierr); 
  ierr = VecSetType(u0, VECMPI);
  CHKERRA(ierr);
  
  ierr = VecSet(u0, zero);
  CHKERRA(ierr);

  //  exsol has type MV_Vector<double> and is input parameter in this method.
  //  Dirichlet values at the boundary 
  //  are assigned to the bvalues in the  ApplyDirichletFEM(exsol)
  ApplyDirichletFEM(exsol);
  
  PetscScalar *boundvalues;
  ierr = VecGetArray(u0, &boundvalues);
  CHKERRA(ierr);
  
  // assign b.values at the inner boundary 
  
  for (int i = 0; i < nbn; i++) {
    boundvalues[boundindex[i]] = bvalues[i];
    //	 boundvalues[boundindex[i]] = -1.0;
    
    //  cout << "  boundvalues(" << i << ") " << boundvalues[boundindex[i]] <<
    //"in node " <<boundindex[i]  << endl;
  }

for (int i = 0; i < nno; i++)
  initbvalues[i] = boundvalues[i];

  ierr = VecRestoreArray(u0, &boundvalues);
    CHKERRA(ierr);

}
//===================================================================
// solution to the laplace u = 0 with non-zeros dirichlet conditions.
// In this method we solve equation A*puas_sol = 0 using PETSc  KSP-solver
// with non-zeros Dirichlet BC.
// Input parameter:
// exsol - is the initial guess from the homogeneous domain with b.c. from the
//         exact measured solution.  Should be prepared.


void WavesPoissonEqOp::PoissonEqSolverFEM(MV_Vector<double>& exsol)
{
 

  Vec b;
  Vec u0;
  Vec u0_left;
  Vec uhelp;
  Vec exact;
  
  PC pc;  /*preconditioner context */
  KSP ksp; /* linear solver context */
 
  int nsd =gg.getNoSpaceDim(); 
  int its;
  PetscScalar one = 1.0;
  PetscScalar two = 2.0;
  PetscScalar minusone = -1.0;
  PetscScalar zero = 0.0;
  
  int nno = gg.getNoNodes();


  VecCreate(PETSC_COMM_SELF, &b);
  ierr = VecSetSizes(b, PETSC_DECIDE, nno); CHKERRA(ierr); 
  ierr = VecSetType(b, VECMPI);
  CHKERRA(ierr);

  VecCreate(PETSC_COMM_SELF, &exact);
  ierr = VecSetSizes(exact, PETSC_DECIDE, nno); CHKERRA(ierr); 
  ierr = VecSetType(exact, VECMPI);
  CHKERRA(ierr); 

  VecCreate(PETSC_COMM_SELF, &puas_sol);
  ierr = VecSetSizes(puas_sol, PETSC_DECIDE, nno); CHKERRA(ierr); 
  ierr = VecSetType(puas_sol, VECMPI);
  CHKERRA(ierr);


  VecCreate(PETSC_COMM_SELF, &u0);
  ierr = VecSetSizes(u0, PETSC_DECIDE, nno); CHKERRA(ierr); 
  ierr = VecSetType(u0, VECMPI);
  CHKERRA(ierr);


  VecCreate(PETSC_COMM_SELF, &u0_left);
  ierr = VecSetSizes(u0_left, PETSC_DECIDE, nno); CHKERRA(ierr); 
  ierr = VecSetType(u0_left, VECMPI);
  CHKERRA(ierr);
  

  VecCreate(PETSC_COMM_SELF, &uhelp);
  ierr = VecSetSizes(uhelp, PETSC_DECIDE, nno); CHKERRA(ierr); 
  ierr = VecSetType(uhelp, VECMPI);
  CHKERRA(ierr);

  ierr = VecSet(u0, zero);
  CHKERRA(ierr);
  ierr = VecSet(u0_left, zero);
  CHKERRA(ierr);
  ierr = VecSet(exact, zero);
  CHKERRA(ierr);

 // cout << " after creating b and puas_sol" << endl;

  //*******************************************************
  // Create the linear solver 
  //********************************************************

  KSPCreate(PETSC_COMM_WORLD, &ksp); /* create linear solver context */
  
  // Set operators. Here the matrix A serves also as the preconditioning matrix

  KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN);

  // Set linear solver defaults for our problem
  KSPGetPC(ksp,&pc); 
  PCSetType(pc, PCJACOBI);
  KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
  KSPSetFromOptions(ksp);


    //  exsol has type MV_Vector<double> and is input parameter in this method.
    //  Dirichlet values at the boundary 
    //  are assigned to the bvalues in the  ApplyDirichletFEM(exsol)
    ApplyDirichletFEM(exsol);
    
    PetscScalar *boundvalues;
    ierr = VecGetArray(u0_left, &boundvalues);
    CHKERRA(ierr);
    
 // assign b.values at the inner boundary 

    for (int i = 0; i < nbn; i++) {
        boundvalues[boundindex[i]] = bvalues[i];
      //	 boundvalues[boundindex[i]] = -1.0;
     
	//  cout << "  boundvalues(" << i << ") " << boundvalues[boundindex[i]] <<
        //"in node " <<boundindex[i]  << endl;
    }


  // set initial guess for KSP solver  as homogeneous medium
  // we need arrange values of exsol for this medium in advance
  // example for the initial guess:  

    PetscScalar *initguess;
    ierr = VecGetArray(puas_sol, &initguess);
    CHKERRA(ierr);

    for (int i = 0; i < nno; i++)
      {
	initguess[i] = exsol(i);
	boundvalues[i] =  exsol(i);
      }


   ierr = VecRestoreArray(u0_left, &boundvalues);
    CHKERRA(ierr);
    ierr = VecRestoreArray(puas_sol, &initguess);
    CHKERRA(ierr);  


  
  //exsol=0.0;
    /* 
  PetscScalar* rhs;
  ierr =  VecGetArray(u0_left,&rhs);CHKERRA(ierr);
  
  for(int i=0; i < gg.getNoNodes(); i++) 
   exsol(i) = rhs[i];

  ierr =  VecRestoreArray(u0_left,&rhs);CHKERRA(ierr);
    */

 
  KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

  // this is not necessary:only left for some future programming with rhs
  if(USE_RHS)
    { 
      if (nsd == 2)
	{
	  
	  evalRHS(b, exsol);
	  
	  // ierr = writeInp("exact.inp", &gg, &exsol, 1); 
	  // CHKERRA(ierr);
	}
    }

//  cout << "before matmulttranspose" << endl;
   ierr= MatScale(A, minusone);

   ierr = MatConvert(A,MATSAME, MAT_INITIAL_MATRIX, &B);  

  ierr= MatScale(B, minusone);

  //ierr=rhsBCxxx(B); CHKERRA(ierr);
  //ierr = MatMultTranspose(A, u0_left, uhelp);
  ierr = MatMult(B, u0_left, uhelp);  //B*u0_left = uhelp
  CHKERRA(ierr);
  // computes uhelp = uhelp + 1*b
 
  ierr = VecAXPY(uhelp, one, b);
  CHKERRA(ierr);
 
//  cout << "before solving Ax=b" << endl;
  //  cout<<"before dirichletBCxxx(A), remove minor with boundary nodes "<<endl;

   ierr=dirichletBCxxx(A); CHKERRA(ierr);

     //********************************************************
// solve equation A*puas_sol = uhelp using PETSc KSP-solver
//***************************************************************
     ierr = KSPSolve(ksp, uhelp, puas_sol);
     CHKERRA(ierr);
//*************************************************************
     PetscScalar *bc;
     ierr = VecGetArray(puas_sol, &bc);
     
     CHKERRA(ierr);
     
    for (int i = 0; i < nbn; i++) {
        bc[boundindex[i]] = bvalues[i];
      //	 boundvalues[boundindex[i]] = -1.0;
     
	// cout << "  boundvalues(" << i << ") " << bc[boundindex[i]] <<
        //"in node " <<boundindex[i]  << endl;
    }
    ierr = VecRestoreArray(puas_sol, &bc);
    CHKERRA(ierr);
   
 /*
    
    PetscScalar* exact_sol;
    ierr =  VecGetArray(exact,&exact_sol);CHKERRA(ierr);
    
    for(int i=0; i < gg.getNoNodes(); i++) 
      exact_sol[i] = exsol(i);
    
    ierr =  VecRestoreArray(exact,&exact_sol);CHKERRA(ierr);
    */
    //     to check difference between computed solution and rhs
 
    // ierr = MatMult(A, puas_sol, u0);  
    //ierr = VecAXPY(uhelp, minusone, u0);
    //print_GID("poisson_diff_q.res",uhelp);


}   




// solve equation A*puas_sol = b using PETSc  KSP-solver
// Input parameters:
// exsol - is the initial guess from the homogeneous domain with b.c. from the
//         exact measured solution.  Should be prepared.
//  sdg - outer FDM geometry in hybrid method
//  grad_tail_square - square of the gradient of the computed tail on the previous pseudo-frequency
//                     iteration. Should be prepared. 


void WavesPoissonEqOp::PoissonEqSolverDirichletFEMnew(MV_Vector<double>& exsol, 
						      WavesSDGeometry& sdg, 
						      MV_Vector<double>& grad_tail_square)
{

  Vec b;
  Vec u0;
  Vec u0_left;
  Vec uhelp;
  Vec exact;
  
  PC pc;  /*preconditioner context */
  KSP ksp; /* linear solver context */
 
  int nsd =gg.getNoSpaceDim(); 
  int its;
  PetscScalar one = 1.0;
  PetscScalar two = 2.0;
  PetscScalar minusone = -1.0;
  PetscScalar zero = 0.0;
  
  int nno = gg.getNoNodes();


  VecCreate(PETSC_COMM_SELF, &b);
  ierr = VecSetSizes(b, PETSC_DECIDE, nno); CHKERRA(ierr); 
  ierr = VecSetType(b, VECMPI);
  CHKERRA(ierr);

  VecCreate(PETSC_COMM_SELF, &exact);
  ierr = VecSetSizes(exact, PETSC_DECIDE, nno); CHKERRA(ierr); 
  ierr = VecSetType(exact, VECMPI);
  CHKERRA(ierr); 

  VecCreate(PETSC_COMM_SELF, &puas_sol);
  ierr = VecSetSizes(puas_sol, PETSC_DECIDE, nno); CHKERRA(ierr); 
  ierr = VecSetType(puas_sol, VECMPI);
  CHKERRA(ierr);


  VecCreate(PETSC_COMM_SELF, &u0);
  ierr = VecSetSizes(u0, PETSC_DECIDE, nno); CHKERRA(ierr); 
  ierr = VecSetType(u0, VECMPI);
  CHKERRA(ierr);


  VecCreate(PETSC_COMM_SELF, &u0_left);
  ierr = VecSetSizes(u0_left, PETSC_DECIDE, nno); CHKERRA(ierr); 
  ierr = VecSetType(u0_left, VECMPI);
  CHKERRA(ierr);
  

  VecCreate(PETSC_COMM_SELF, &uhelp);
  ierr = VecSetSizes(uhelp, PETSC_DECIDE, nno); CHKERRA(ierr); 
  ierr = VecSetType(uhelp, VECMPI);
  CHKERRA(ierr);

  /* old petsc 2.1.0 format
     
  ierr = VecSet(&zero, u0);
  CHKERRA(ierr);
  ierr = VecSet(&zero, u0_left);
  CHKERRA(ierr);
  ierr = VecSet(&zero, exact);
 CHKERRA(ierr);
  
  */

  ierr = VecSet(u0, zero);
  CHKERRA(ierr);
  ierr = VecSet(u0_left, zero);
  CHKERRA(ierr);
  ierr = VecSet(exact, zero);
  CHKERRA(ierr);

//  cout << " after creating b and puas_sol" << endl;

  //*******************************************************
  // Create the linear solver 
  //********************************************************


  KSPCreate(PETSC_COMM_WORLD, &ksp); /* create linear solver context */
  
  // Set operators. Here the matrix A serves also as the preconditioning matrix

  KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN);

  // Set linear solver defaults for our problem
  KSPGetPC(ksp,&pc); 
  PCSetType(pc, PCJACOBI);
  KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
  KSPSetFromOptions(ksp);


  if (USE_DIRICHLET_FEM) {
    
    //  exsol has type MV_Vector<double> and is input parameter in this method.
    //  Dirichlet and Neumann values at the boundary 
    //  are assigned to the bvalues in the  ApplyDirichletFEM(exsol)
    ApplyDirichletFEM(exsol);
    
    PetscScalar *boundvalues;
    ierr = VecGetArray(u0_left, &boundvalues);
    CHKERRA(ierr);
    


    // assign b.values at the inner boundary 
    
    for (int i = 0; i < nbn; i++) {
      
      boundvalues[boundindex[i]] = bvalues[i];
      //  boundvalues[boundindex[i]] = bvalues[i] + grad_tail(boundindex[i]);
      //	 boundvalues[boundindex[i]] = -1.0;
      
//      cout << "  boundvalues(" << i << ") " <<boundvalues[boundindex[i]]<<"in global  node " <<boundindex[i]<< endl;
    }



  // set initial guess for KSP solver  as homogeneous medium
  // we need arrange values of p for this medium in advance
  // example for initial guess:  

    PetscScalar *initguess;
    ierr = VecGetArray(puas_sol, &initguess);
    CHKERRA(ierr);

    for (int i = 0; i < nno; i++)
      {
	initguess[i] = exsol(i);
	//	boundvalues[i] =  exsol(i);
      }


   ierr = VecRestoreArray(u0_left, &boundvalues);
    CHKERRA(ierr);
    ierr = VecRestoreArray(puas_sol, &initguess);
    CHKERRA(ierr);  

    
    KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

 /* // commented by Thanh
    GIDOutputMesh   Vmesh(&sdg,"u0.msh");
    bool  m = Vmesh.printMesh();
  

    if (nsd ==2)
      {
	GIDOutputOp   gid_out2(&sdg,"u0.res",  initguess,   boundvalues , initguess, 1);
	bool pr = gid_out2.printResults();
	
	
      }
    else if (nsd ==3)
      { 
	GIDOutputOp3D   gid1(&sdg,"u0.res",   initguess, initguess ,   boundvalues, initguess, 1);
	bool  pr1 = gid1.printResults();
      }
*/    

    /* 
       for(int i=0; i < gg.getNoNodes(); i++) 
       {
       grad_tail_square(i) = grad_tail_square(i) + boundvalues[i];
       //    grad_tail_square(i) = -grad_tail_square(i) ;
       
       }
       
    */
 
    
  }
  
  
  if(USE_RHS)
	   { 
            
		// cout<<"in evalRHS"<<endl;
		     evalRHS(b, grad_tail_square);
		    
		     // ierr = writeInp("exact.inp", &gg, &exsol, 1); 
		     // CHKERRA(ierr);
	   }

 // cout << "before matmulttranspose" << endl;
  Mat C;
  ierr = MatConvert(A,MATSAME, MAT_INITIAL_MATRIX, &C);    
  //old petsc 2.1.0 
  //ierr= MatScale(&minusone,A);
  
  ierr= MatScale(A, minusone);

  ierr = MatMult(C, u0_left, uhelp); // C*u0_left = uhelp
 
  // computes b= 1*uhelp + b
  ierr = VecAXPY(b,one, uhelp);
    
  // CHKERRA(ierr);
 
//  cout << "before solving  A*puas_sol = b" << endl;
  // cout<<"before dirichletBCxxx(A), remove minor with boundary nodes "<<endl;

   ierr=dirichletBCxxx(A); CHKERRA(ierr);
 
  // solve equation A*puas_sol = uhelp using PETSc KSP-solver
    //	ierr = KSPSolve(ksp, uhelp, puas_sol);
// solve equation A*puas_sol = b using PETSc KSP-solver
	 	ierr = KSPSolve(ksp, b, puas_sol);
	
	CHKERRA(ierr);

}   

/*
void WavesPoissonEqOp::ComputeL2NormPoisson()
{
  double l2norm;
  PetscScalar* poisson;
  ierr =  VecGetArray(puas_sol,&poisson);CHKERRA(ierr);
  
  l2norm = computeL2Norm(poisson); 

  cout<<" norm is "<<l2norm<<endl;


  // opt.ComputeL2Norm(poisson, l2norm);

 
  ierr =  VecRestoreArray(puas_sol,&poisson);CHKERRA(ierr);
  


}

*/




void WavesPoissonEqOp::writeSol(char* filename)
{
  ofstream outp;
  int i;
  
  outp.open(filename); 
  // if (!fp) SETERRA(1,1,"Cannot open matlab file\n");

  int nno = gg.getNoNodes();
  
  PetscScalar* puas;
  ierr =  VecGetArray(puas_sol,&puas);CHKERRA(ierr);
  
  for(i=0;i<nno;i++){ 
    outp<<puas[i]<<"\n";}
  outp.close();
  
  ierr =  VecRestoreArray(puas_sol,&puas);CHKERRA(ierr);
  
  //print("Solution at one point is written to file: %s\n", filename);
}

void WavesPoissonEqOp::writeSol(char* filename, Vec& sol)
{
  ofstream outp;
  int i;
  
  outp.open(filename); 
  // if (!fp) SETERRA(1,1,"Cannot open matlab file\n");

  int nno = gg.getNoNodes();
  
	 PetscScalar* tmp;
	 ierr =  VecGetArray(sol,&tmp);CHKERRA(ierr);
	  
	 for(i=0;i<nno;i++){ 
	   outp<<tmp[i]<<"\n";}
	 outp.close();
	 
	 ierr =  VecRestoreArray(sol,&tmp);CHKERRA(ierr);
	
  //printf("Solution at one point is written to file: %s\n", filename);
}


// save the solution to a given vector of type MV_Vector<double>. Added by Thanh
void WavesPoissonEqOp::extract_solution(MV_Vector<double>& Sol)
{
  int i;  
  int nno = gg.getNoNodes();
  
  PetscScalar* puas;
  ierr =  VecGetArray(puas_sol,&puas);CHKERRA(ierr);
  
  for(i=0;i<nno;i++)
	Sol(i) = puas[i];
  
  ierr =  VecRestoreArray(puas_sol,&puas);CHKERRA(ierr);
 // cout << " The solution has been saved to vector Sol: " << endl;

}

//added by Thanh
void  WavesPoissonEqOp::DirichletFEM(MV_Vector<double> exsol, MV_Vector<double>& function)
{
  
  //  exsol - values at the boundary, are assigned to bvalues in next function
  ApplyDirichletFEM(exsol);
  
  for (int i = 0; i < nbn; i++) {
    function(boundindex[i]) = bvalues[i];
    //	 boundvalues[boundindex[i]] = -1.0;
 //    cout<<"  boundary values "<<function[boundindex[i]] <<endl;
    
  }
  
}
