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
// Last changed: 2012-01-09 by Vladimir Timonov

#include "include/wavesOutputs.h"

using namespace std;

WavesOutputs::WavesOutputs()
{
}

WavesOutputs::WavesOutputs(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max) :
		x_min_(x_min), x_max_(x_max), y_min_(y_min), y_max_(y_max), z_min_(z_min), z_max_(z_max)
{
}

WavesOutputs::WavesOutputs(Grid& gg_, MV_Vector<double>& bb_)
{
	gg = gg_;
	bb = bb_;
}

WavesOutputs::WavesOutputs(Grid& gg_, MV_Vector<int>& b_)
{
	gg = gg_;
	b = b_;
}

WavesOutputs::WavesOutputs(Grid& gg_)
{
	gg = gg_;

}

WavesOutputs::WavesOutputs(WavesSDGeometry& sdg_)
{
	sdg = sdg_;

}

WavesOutputs::WavesOutputs(double norma_, int it_)
{
	norma = norma_;
	it = it_;
}

WavesOutputs::~WavesOutputs()

{
//	printf("WavesOutputs::~WavesOutputs: deleting object");
}

void WavesOutputs::CreateVecGlobNum(char *fname, MV_Vector<int>& Vec_Glob_Num)
{
	ifstream inp;
	inp.open(fname);
	int i, nonodes;
	//printf("Open the file: %s\n", fname);

	// nnodes to be supplied
	inp >> nonodes;

	Vec_Glob_Num.newsize(nonodes);

	for (i = 0; i < nonodes; i++)
		inp >> Vec_Glob_Num(i);

	inp.close();

}
void WavesOutputs::WriteMyPart(char* file, double func, int iter)
{

	int i, j;

	FILE *fp;
	fp = fopen(file, "a+");

	if (!fp)
		cout << endl << "WavesOutputs:in writemypart: cannot open file \"" << file << "\"\n";
	cout << endl << "opened file  " << file << endl;

	fseek(fp, offset(iter), SEEK_SET);

	if (iter == 1)
		fprintf(fp, "started new computations %f", func);
	else
		fprintf(fp, "%f ", func);

	cout << endl << "closed file" << file << endl;
	fclose(fp);
}

void WavesOutputs::WriteMyPart(char* file, int nr, double x, double y, double z, int iter, double norma)
{

	int i, j;

	FILE *fp;
	fp = fopen(file, "a+");

	if (!fp)
		cout << endl << "WavesOutputs:in writemypart: cannot open file \"" << file << "\"\n";
	cout << endl << "opened file  " << file << endl;

	fseek(fp, offset(iter), SEEK_SET);

	fprintf(fp, "  %i %f  %f %f %f \n", nr, x, y, z, norma);

	cout << endl << "closed file" << file << endl;
	fclose(fp);
}

int WavesOutputs::writeInp(char *file, Grid *grid, double* u, int nsys)
{
	FILE *fp;
	int i, j;
	int el;
	int NODE_DATA = 1;

	fp = fopen(file, "w");

	fprintf(fp, "%i %i %i 0 0\n", grid->getNoNodes(), grid->getNoElms(), nsys);

	for (i = 0; i < grid->getNoNodes(); i++)
	{
		fprintf(fp, "%i %f %f %f\n", i + 1, grid->getCoor(i, 0), grid->getCoor(i, 1), grid->getCoor(i, 2));
	}
	for (i = 0; i < grid->getNoElms(); i++)
	{
		fprintf(fp, "%i %i tet %i %i %i %i\n", i + 1, 1, grid->loc2glob(i, 0) + 1, grid->loc2glob(i, 1) + 1, grid->loc2glob(i, 2) + 1, grid->loc2glob(i, 3) + 1);
	}

	if (NODE_DATA)
	{
		fprintf(fp, "%i ", nsys);
		for (i = 0; i < nsys; i++)
			fprintf(fp, "1 ");
		fprintf(fp, "\n");
		for (i = 0; i < nsys; i++)
			fprintf(fp, "u %i, Unit\n", i + 1);
		for (i = 0; i < grid->getNoNodes(); i++)
		{
			fprintf(fp, "%i ", i + 1);
			for (j = 0; j < nsys; j++)
				fprintf(fp, "%f ", u[j + i * nsys]);
			fprintf(fp, "\n");
		}
	}
	else
	{
		fprintf(fp, "1 1");
		fprintf(fp, "\n");
		fprintf(fp, "u , Unit\n");
		for (i = 0; i < grid->getNoElms(); i++)
		{
			fprintf(fp, "%i ", i + 1);
			fprintf(fp, "%i \n", grid->getMaterialType(i));
		}
	}
	fclose(fp);

//	printf("Solution written to file: %s\n", file);

	return 0;
}

void WavesOutputs::WriteSolutionsToFile(char *fname, Mat_real& Ar1)
{
	ofstream outp;
	int i, j;

	outp.open(fname);

	// rows are vectors with values of the function
	// columns are timesteps

	for (j = 0; j < Ar1.size(1); j++)
	{
		for (i = 0; i < Ar1.size(0); i++)
		{

			outp << Ar1(i, j) << "   ";
		}
		outp << "\n";
	}

	outp.close();

	//print("Solution at one point is written to file: %s\n", fname);

}

void WavesOutputs::WriteToFile(char *fname, MV_Vector<double>& Ar1)
{
	ofstream outp;
	int i, j;

	outp.open(fname);

	for (i = 0; i < Ar1.size(); i++)
		outp << Ar1(i) << "   ";

	outp.close();

	//print("Solution at one point is written to file: %s\n", fname);
}

void WavesOutputs::WriteToFile(char *fname, MV_Vector<bool>& Ar1)
{
	ofstream outp;
	int i, j;

	outp.open(fname);

	for (i = 0; i < Ar1.size(); i++)
		outp << Ar1(i) << "   ";

	outp.close();

	//print("Solution at one point is written to file: %s\n", fname);
}

void WavesOutputs::WriteToFile(char *fname, MV_Vector<int>& Ar1)
{
	ofstream outp;
	int i, j;

	outp.open(fname);

	for (i = 0; i < Ar1.size(); i++)
		outp << Ar1(i) << "   ";

	outp.close();

	//print("Solution at one point is written to file: %s\n", fname);
}

void WavesOutputs::WriteToFile(char *fname, real* Ar1, int nno)
{
	ofstream outp;
	int i, j;

	outp.open(fname);

	for (i = 0; i < nno; i++)
		outp << Ar1[i] << "   ";

	outp.close();

	//print("Solution at one point is written to file: %s\n", fname);
}

//=============================================================
void WavesOutputs::ReadConst(char *fname, double& alfa)
{
	ifstream inp;
	inp.open(fname);
	int i;
	//printf("Open the file: %s\n", fname);

	inp >> alfa;

	inp.close();

}

//=========================================================================
void WavesOutputs::ReadSolfromFile(char *fname, Mat_real& Ar1)
{
	ifstream inp;
	inp.open(fname);
	int i, j;
	printf("Open the file: %s\n", fname);
	cout << endl << "number of columns" << Ar1.size(0) << "number of rows  " << Ar1.size(1) << endl;
	// rows are vectors with values of the function:
	//  Ar1.size(0) - number of nodes
	// columns are timesteps:  Ar1.size(1) - number of timesteps

	for (j = 0; j < Ar1.size(1); j++)
	{
		for (i = 0; i < Ar1.size(0); i++)
		{
			inp >> Ar1(i, j);

			if (Ar1(i, j) > 10000)
				cout << endl << " very big number: " << Ar1(i, j) << "in the point " << i << "and time moment" << j << endl;
		}
	}

	inp.close();
}

void WavesOutputs::ReadSolfromFile(char *fname, MV_Vector<double>& Ar1)
{
	ifstream inp;
	inp.open(fname);
	int i;
	//printf("Open the file: %s\n", fname);

	for (i = 0; i < Ar1.size(); i++)
	{
		inp >> Ar1(i);
		//cout << endl << " values of function" << Ar1(i) << endl;
	}

	inp.close();

}

void WavesOutputs::ReadSolfromFile(char *fname, MV_Vector<int>& Ar1)
{
	ifstream inp;
	inp.open(fname);
	int i;
	//printf("Open the file: %s\n", fname);

	for (i = 0; i < Ar1.size(); i++)
	{
		inp >> Ar1(i);
	}

	inp.close();

}

void WavesOutputs::WriteToGidFile1(char *fname, WavesGridB& gg, MV_Vector<int>& Vec_glob_Num)
{
	ofstream outp;
	outp.open(fname);

	int e, nel, sch;
	sch = 0;

	int n1, n2, n3, n4;

	nel = gg.getNoElms();
	int nno = gg.getNoNodes();

	Mat_int Gid_Elms(nel, 5);

	double x_max = 18;
	double y_max = 12.0;
	double z_max = 7.0;
	double x_min = 12.0;
	double y_min = 10.0;
	double z_min = 5.0;

	for (e = 0; e < nel; e++)
	{

		if (gg.getNoSpaceDim() == 3)
		{
			n1 = gg.loc2glob(e, 0);
			n2 = gg.loc2glob(e, 1);
			n3 = gg.loc2glob(e, 2);
			n4 = gg.loc2glob(e, 3);

			if (gg.getCoor(n1, 0) > x_min && gg.getCoor(n1, 0) < x_max && gg.getCoor(n2, 0) > x_min && gg.getCoor(n2, 0) < x_max && gg.getCoor(n3, 0) > x_min && gg.getCoor(n3, 0) < x_max && gg.getCoor(n4, 0) > x_min && gg.getCoor(n4, 0) < x_max
					&& gg.getCoor(n1, 1) > y_min && gg.getCoor(n1, 1) < y_max && gg.getCoor(n2, 1) > y_min && gg.getCoor(n2, 1) < y_max && gg.getCoor(n3, 1) > y_min && gg.getCoor(n3, 1) < y_max && gg.getCoor(n4, 1) > y_min && gg.getCoor(n4, 1) < y_max
					&& gg.getCoor(n1, 2) > z_min && gg.getCoor(n1, 2) < z_max && gg.getCoor(n2, 2) > z_min && gg.getCoor(n2, 2) < z_max && gg.getCoor(n3, 2) > z_min && gg.getCoor(n3, 2) < z_max && gg.getCoor(n4, 2) > z_min && gg.getCoor(n4, 2) < z_max)
			{
				Gid_Elms(sch, 0) = e;
				Gid_Elms(sch, 1) = n1;
				Gid_Elms(sch, 2) = n2;
				Gid_Elms(sch, 3) = n3;
				Gid_Elms(sch, 4) = n4;
				sch++;
				cout << endl << " nodes " << n1 << "  " << n2 << "  " << n3 << "   " << n4 << endl;

			}
		}
	}

	MV_Vector<int> Markers(nno);

	int i;
	int new_sch = 0;
	for (i = 0; i < sch; i++)
	{

		n1 = Gid_Elms(i, 1);
		n2 = Gid_Elms(i, 2);
		n3 = Gid_Elms(i, 3);
		n4 = Gid_Elms(i, 4);

		Markers(n1) = 1;
		Markers(n2) = 1;
		Markers(n3) = 1;
		Markers(n4) = 1;

	}

	for (i = 0; i < nno; i++)
		if ((Markers(i)) == 1)
			new_sch++;

	Vec_glob_Num.newsize(new_sch);

	outp << new_sch << endl;

	new_sch = 0;

	for (i = 0; i < nno; i++)
		if ((Markers(i)) == 1)
		{
			outp << i << endl;
			Vec_glob_Num(new_sch) = i;
			new_sch++;
		}
	outp.close();

	printf("Observation points are  written to file: %s\n", fname);
}

//==============================================================

void WavesOutputs::WriteToGidFile(char *fname, WavesGridB& gg, MV_Vector<int>& Vec_glob_Num, int code)
{
	ofstream outp;
	outp.open(fname);

	int e, nel, sch;
	sch = 0;

	int n1, n2, n3, n4;

	nel = gg.getNoElms();
	int nno = gg.getNoNodes();

	Mat_int Gid_Elms(nel, 5);

	// code == 1  obs. points located at the down boundary
	// code == 2  obs. points located at the right boundary
	// code == 3  obs. points located at the back boundary
	double x_max;
	double y_max;
	double z_max;
	double x_min;
	double y_min;
	double z_min;

	if (code == 1)
	{
		x_max = 18.0;
		y_max = 11.0;
		z_max = 7.0;
		x_min = 12.0;
		y_min = 10.0;
		z_min = 4.0;
	}
	else if (code == 2)
	{
		x_max = 20.0;
		y_max = 18.0;
		z_max = 7.0;
		x_min = 18.0;
		y_min = 12.0;
		z_min = 4.0;
	}
	else if (code == 3)
	{
		x_max = 18.0;
		y_max = 17.0;
		z_max = 10.0;
		x_min = 12.0;
		y_min = 15.0;
		z_min = 9.2;
	}

	for (e = 0; e < nel; e++)
	{

		if (gg.getNoSpaceDim() == 3)
		{
			n1 = gg.loc2glob(e, 0);
			n2 = gg.loc2glob(e, 1);
			n3 = gg.loc2glob(e, 2);
			n4 = gg.loc2glob(e, 3);

			if (gg.getCoor(n1, 0) > x_min && gg.getCoor(n1, 0) < x_max && gg.getCoor(n2, 0) > x_min && gg.getCoor(n2, 0) < x_max && gg.getCoor(n3, 0) > x_min && gg.getCoor(n3, 0) < x_max && gg.getCoor(n4, 0) > x_min && gg.getCoor(n4, 0) < x_max
					&& gg.getCoor(n1, 1) > y_min && gg.getCoor(n1, 1) < y_max && gg.getCoor(n2, 1) > y_min && gg.getCoor(n2, 1) < y_max && gg.getCoor(n3, 1) > y_min && gg.getCoor(n3, 1) < y_max && gg.getCoor(n4, 1) > y_min && gg.getCoor(n4, 1) < y_max
					&& gg.getCoor(n1, 2) > z_min && gg.getCoor(n1, 2) < z_max && gg.getCoor(n2, 2) > z_min && gg.getCoor(n2, 2) < z_max && gg.getCoor(n3, 2) > z_min && gg.getCoor(n3, 2) < z_max && gg.getCoor(n4, 2) > z_min && gg.getCoor(n4, 2) < z_max)
			{
				Gid_Elms(sch, 0) = e;
				Gid_Elms(sch, 1) = n1;
				Gid_Elms(sch, 2) = n2;
				Gid_Elms(sch, 3) = n3;
				Gid_Elms(sch, 4) = n4;
				sch++;
				cout << endl << " nodes " << n1 << "  " << n2 << "  " << n3 << "   " << n4 << endl;

			}
		}
	}

	MV_Vector<int> Markers(nno);

	int i;
	int new_sch = 0;
	for (i = 0; i < sch; i++)
	{

		n1 = Gid_Elms(i, 1);
		n2 = Gid_Elms(i, 2);
		n3 = Gid_Elms(i, 3);
		n4 = Gid_Elms(i, 4);

		Markers(n1) = 1;
		Markers(n2) = 1;
		Markers(n3) = 1;
		Markers(n4) = 1;

	}

	for (i = 0; i < nno; i++)
		if ((Markers(i)) == 1)
			new_sch++;

	Vec_glob_Num.newsize(new_sch);

	outp << new_sch << endl;

	new_sch = 0;

	for (i = 0; i < nno; i++)
		if ((Markers(i)) == 1)
		{
			outp << i << endl;
			Vec_glob_Num(new_sch) = i;
			new_sch++;
		}
	outp.close();

	printf("Observation points are  written to file: %s\n", fname);
}

//=================================================================

void WavesOutputs::WriteWavesElementArrays(WavesGridB& gg, Mat_int& Gid_Elms, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, int& sch)
{
	int e, n1, n2, n3, n4;
	int nel = gg.getNoElms();
	int nno = gg.getNoNodes();

	for (e = 0; e < nel; e++)
	{

		if (gg.getNoSpaceDim() == 3)
		{
			n1 = gg.loc2glob(e, 0);
			n2 = gg.loc2glob(e, 1);
			n3 = gg.loc2glob(e, 2);
			n4 = gg.loc2glob(e, 3);

			if (gg.getCoor(n1, 0) > x_min && gg.getCoor(n1, 0) < x_max && gg.getCoor(n2, 0) > x_min && gg.getCoor(n2, 0) < x_max && gg.getCoor(n3, 0) > x_min && gg.getCoor(n3, 0) < x_max && gg.getCoor(n4, 0) > x_min && gg.getCoor(n4, 0) < x_max
					&& gg.getCoor(n1, 1) > y_min && gg.getCoor(n1, 1) < y_max && gg.getCoor(n2, 1) > y_min && gg.getCoor(n2, 1) < y_max && gg.getCoor(n3, 1) > y_min && gg.getCoor(n3, 1) < y_max && gg.getCoor(n4, 1) > y_min && gg.getCoor(n4, 1) < y_max
					&& gg.getCoor(n1, 2) > z_min && gg.getCoor(n1, 2) < z_max && gg.getCoor(n2, 2) > z_min && gg.getCoor(n2, 2) < z_max && gg.getCoor(n3, 2) > z_min && gg.getCoor(n3, 2) < z_max && gg.getCoor(n4, 2) > z_min && gg.getCoor(n4, 2) < z_max)
			{
				Gid_Elms(sch, 0) = e;
				Gid_Elms(sch, 1) = n1;
				Gid_Elms(sch, 2) = n2;
				Gid_Elms(sch, 3) = n3;
				Gid_Elms(sch, 4) = n4;
				sch++;
				cout << endl << " nodes " << n1 << "  " << n2 << "  " << n3 << "   " << n4 << " sch " << sch << endl;

			}
		}
	}

}

void WavesOutputs::WriteToGidFile2(char *fname, WavesGridB& gg, MV_Vector<int>& Vec_glob_Num)
{
	ofstream outp;
	outp.open(fname);

	int e, nel, sch;
	sch = 0;

	int n1, n2, n3, n4;

	nel = gg.getNoElms();
	int nno = gg.getNoNodes();

	Mat_int Gid_Elms(nel, 5);

	// code == 1  obs. points located at the down boundary
	// code == 2  obs. points located at the right boundary
	// code == 3  obs. points located at the back boundary

	double x_max;
	double y_max;
	double z_max;
	double x_min;
	double y_min;
	double z_min;

	x_max = 16.0;
	y_max = 10.5;
	z_max = 7.0;
	x_min = 14.0;
	y_min = 10.0;
	z_min = 5.0;

	cout << endl << " write points at the down boundary " << endl;

	WriteWavesElementArrays(gg, Gid_Elms, x_min, x_max, y_min, y_max, z_min, z_max, sch);

	x_max = 20.0;
	y_max = 16.0;
	z_max = 7.0;
	x_min = 19.0;
	y_min = 14.0;
	z_min = 5.0;

	cout << endl << " write points at the right boundary " << endl;

	WriteWavesElementArrays(gg, Gid_Elms, x_min, x_max, y_min, y_max, z_min, z_max, sch);

	x_max = 11.0;
	y_max = 17.0;
	z_max = 7.0;
	x_min = 10.0;
	y_min = 15.0;
	z_min = 5.0;

	cout << endl << " write points at the left boundary " << endl;

	WriteWavesElementArrays(gg, Gid_Elms, x_min, x_max, y_min, y_max, z_min, z_max, sch);

	MV_Vector<int> Markers(nno);

	int i;
	int new_sch = 0;

	cout << endl << " sch " << sch << endl;

	for (i = 0; i < sch; i++)
	{

		n1 = Gid_Elms(i, 1);
		n2 = Gid_Elms(i, 2);
		n3 = Gid_Elms(i, 3);
		n4 = Gid_Elms(i, 4);

		Markers(n1) = 1;
		Markers(n2) = 1;
		Markers(n3) = 1;
		Markers(n4) = 1;

	}

	for (i = 0; i < nno; i++)
		if ((Markers(i)) == 1)
			new_sch++;

	Vec_glob_Num.newsize(new_sch);

	outp << new_sch << endl;

	new_sch = 0;

	for (i = 0; i < nno; i++)
		if ((Markers(i)) == 1)
		{
			outp << i << endl;
			Vec_glob_Num(new_sch) = i;
			new_sch++;
		}
	outp.close();

	printf("Observation points are  written to file: %s\n", fname);
}

//================================================================
void WavesOutputs::WriteToGidFile(char *fname, WavesSDGeometry& sdg, MV_Vector<int>& Vec_glob_Num)
{
	ofstream outp;
	outp.open(fname);

	int i, nel, sch;
	sch = 0;

	int n;

	nel = sdg.getNoElms();
	int nno = sdg.getNoNodes();

	MV_Vector<int> Gid_Elms(nno);

	// on the down boundary, mesh size h=0.2
	double x_max = 4.7;
	double y_max = 0.5;
	double z_max = 1.5;
	double x_min = 0.3;
	double y_min = 0.1;
	double z_min = 1.1;

	for (n = 0; n < nno; n++)
	{
		if (sdg.getNoSpaceDim() == 3)
		{

			if ((sdg.getCoor(n, 0) > x_min && sdg.getCoor(n, 0) < x_max) && (sdg.getCoor(n, 1) > y_min && sdg.getCoor(n, 1) < y_max) && (sdg.getCoor(n, 2) > z_min && sdg.getCoor(n, 2) < z_max))
			{

				Gid_Elms(sch) = n;

				sch++;
				cout << endl << " nodes " << n << endl;

			}
		}
	}

	Vec_glob_Num.newsize(sch);
	outp << sch << endl;

	for (i = 0; i < sch; i++)
	{

		outp << Gid_Elms(i) << endl;

		Vec_glob_Num(i) = Gid_Elms(i);

	}

	outp.close();
	printf("Observation points are  written to file: %s\n", fname);
}

void WavesOutputs::GidWrite(char *fname, WavesSDGeometry& sdg, MV_Vector<int>& Vec_glob_Num)
{

	FILE *fp;
	fp = fopen(fname, "a+");

	if (!fp)
	{
		cout << endl << " cannot open file \"" << fname << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	int i, nel, sch;
	sch = 0;

	int n;

	nel = sdg.getNoElms();
	int nno = sdg.getNoNodes();

	MV_Vector<int> Gid_Elms(nno);

	for (n = 0; n < nno; n++)
	{
		if (sdg.getNoSpaceDim() == 3)
		{

			if ((sdg.getCoor(n, 0) > x_min_ && sdg.getCoor(n, 0) < x_max_) && (sdg.getCoor(n, 1) > y_min_ && sdg.getCoor(n, 1) < y_max_) && (sdg.getCoor(n, 2) > z_min_ && sdg.getCoor(n, 2) < z_max_))
			{

				Gid_Elms(sch) = n;

				sch++;
				cout << endl << " nodes " << n << endl;

			}
		}
	}

	Vec_glob_Num.newsize(sch);

	fprintf(fp, " %i \n ", sch);

	for (i = 0; i < sch; i++)
	{
		fprintf(fp, " %i \n ", Gid_Elms(i));

		Vec_glob_Num(i) = Gid_Elms(i);

	}

	fclose(fp);
	printf("Observation points are  written to file: %s\n", fname);
}

//============== write results to gid file for FDM ================

int WavesOutputs::print_GID_FDM(char* filename, int sch, MV_Vector<double>& E1_new, MV_Vector<double>& E2_new)
{

	int i, j;

	int code, ierr;
	bool pr;
	int nno = sdg.getNoNodes();
	int nsd = sdg.getNoSpaceDim();

	real* E_common_FDM = new real[nno];
	real* E1 = new real[nno];
	real* E2 = new real[nno];

	// HRN
	ierr = 0;

	for (i = 0; i < nno; i++)
	{
		E_common_FDM[i] = sqrt(E1_new(i) * E1_new(i) + E2_new(i) * E2_new(i));
		E1[i] = E1_new(i);
		E2[i] = E2_new(i);
	}
	if (nsd == 2)
	{

		// here, we print out displacement of the elastic vector in  2D
		// Vector_results is results title
		// 2 - type of analysis (displacement analysis)
		// %i - time step
		// 2 - kind of results - vector
		// 1 - position of the data - in nodes
		// 0 - no description inside

		GIDOutputOp gid_out1(&sdg, filename, E1, E2, E_common_FDM, sch);

		pr = gid_out1.printResults();

	}

	delete[] E1;
	delete[] E2;
	delete[] E_common_FDM;

	return ierr;
}

//============= in 3D =============================================

int WavesOutputs::print_GID_FDM3D(char* filename, int sch, MV_Vector<double>& E1_new, MV_Vector<double>& E2_new, MV_Vector<double>& E3_new)
{
	int i, j;
	int code, ierr;
	bool pr;
	int nno = sdg.getNoNodes();
	int nsd = sdg.getNoSpaceDim();

	real* E_common_FDM = new real[nno];
	real* E1 = new real[nno];
	real* E2 = new real[nno];
	real* E3 = new real[nno];

	// HRN
	ierr = 0;

	for (i = 0; i < nno; i++)
	{
		E_common_FDM[i] = sqrt(E1_new(i) * E1_new(i) + E2_new(i) * E2_new(i) + E3_new(i) * E3_new(i));
		E1[i] = E1_new(i);
		E2[i] = E2_new(i);
		E3[i] = E3_new(i);
	}
	if (nsd == 3)
	{
		// here, we print out displacement of the elastic vector in  3D
		// Vector_results is results title
		// 2 - type of analysis (displacement analysis)
		// %i - time step
		// 2 - kind of results - vector
		// 1 - position of the data - in nodes
		// 0 - no description inside

		GIDOutputOp3D gid_out1(&sdg, filename, E1, E2, E3, E_common_FDM, sch);

		pr = gid_out1.printResults();

	}

	delete[] E1;
	delete[] E2;
	delete[] E3;
	delete[] E_common_FDM;
	return ierr;
}

// another version with Mat_real inputs

int WavesOutputs::print_GID_FDM(char* filename, int sch, Mat_real& E)
{

	int i, j;

	int code, ierr;
	bool pr;
	int nno = sdg.getNoNodes();
	int nsd = sdg.getNoSpaceDim();

	real* E_common_FDM = new real[nno];
	real* E1 = new real[nno];
	real* E2 = new real[nno];

	// HRN
	ierr = 0;

	for (i = 0; i < nno; i++)
	{
		E_common_FDM[i] = sqrt(E(i, 0) * E(i, 0) + E(i, 1) * E(i, 1));
		E1[i] = E(i, 0);
		E2[i] = E(i, 1);
	}
	if (nsd == 2)
	{

		// here, we print out displacement of the elastic vector in  2D
		// Vector_results is results title
		// 2 - type of analysis (displacement analysis)
		// %i - time step
		// 2 - kind of results - vector
		// 1 - position of the data - in nodes
		// 0 - no description inside

		GIDOutputOp gid_out1(&sdg, filename, E1, E2, E_common_FDM, sch);

		pr = gid_out1.printResults();

	}

	delete[] E1;
	delete[] E2;
	delete[] E_common_FDM;

	return ierr;
}

//====================================================================
// write results file for vizualization in Gid-postprocessor
// input parametrs: name of the file (example: file.flavia.res)
//                   results, number of iteration k
//  Values of the parameter "results" should be:
//                         1 - if "results" is  scalar results  ( in 
//                              output file  the resuts vector have
//                              following format:
//                              i   result(i) )
//                         2  - if "results" is vector, format:
//                         i   x(i) y(i) z(i), 
//                       where x(i),y(i),z(i) are results at the each
//                       component (x,y,z) at the node i.
//                         3 - if "results" is matrix , output formats
//                         i  M_xx(i) M_yy(i) M_zz(i) M_xy(i) M_yz(i) M_xz(i),
//                         where i is global node and M_... are values
// of the partial derivatives _xx,_yy, ... at this node. Used for stresses.
//=========================================================================== 
bool WavesOutputs::Inp_to_RES_Gid(char *inp_file, int results, int k)
{
	FILE *fp;

	int i, j;
	int el;

	int nsd = gg.getNoSpaceDim();

	Mat_real displ;

	cout << endl << " results " << results << "nsd " << nsd << "  k " << k << endl;
	if (results == 1)
		for (i = 0; i < gg.getNoNodes(); i++)
		{
			displ.newsize(gg.getNoNodes(), 1);

			displ(i, 0) = bb[i];

		}
	else if (results == 2)
	{
		displ.newsize(gg.getNoNodes(), 3);

		displ(i, 0) = bb[i];
		displ(i, 1) = bb[i];
		displ(i, 2) = bb[i];
	}
	else if (results == 3)
	{
		displ.newsize(gg.getNoNodes(), 6);

		displ(i, 0) = bb[i];
		displ(i, 1) = bb[i];
		displ(i, 2) = bb[i];
		displ(i, 3) = bb[i];
		displ(i, 4) = bb[i];
		displ(i, 5) = bb[i];
	}

	fp = fopen(inp_file, "a+");

	if (!fp)
	{
		cout << endl << " cannot open file \"" << inp_file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	// write 1 line to gid file


	if (results == 1)
	{
		fprintf(fp, " Scalar_result   1    %i       1   1   0 \n ", k);

		for (i = 0; i < gg.getNoNodes(); i++)
		{
			fprintf(fp, "%i %f \n", i + 1, displ(i, 0));
		}
	}
	else if (results == 2)
	{
		fprintf(fp, " Vector_result   1    %i       2   1   0 \n ", k);

		for (i = 0; i < gg.getNoNodes(); i++)
		{
			fprintf(fp, "%i %f %f %f\n", i + 1, displ(i, 0), displ(i, 1), displ(i, 2));
		}
	}
	else if (results == 3)
	{
		fprintf(fp, " Matrix_result   1    %i       3   1   0 \n ", k);

		for (i = 0; i < gg.getNoNodes(); i++)
		{
			fprintf(fp, "%i %f %f %f %f %f %f\n", i + 1, displ(i, 0), displ(i, 1), displ(i, 2), displ(i, 3), displ(i, 4), displ(i, 5));
		}
	}

	fprintf(fp, "\n");

	fclose(fp);

//	printf("Solution written to file: %s\n", inp_file);

	return true;

}
//==========================================================================

void WavesOutputs::WriteToFile(char *fname, WavesSDGeometry& sdg, WavesGridB& gg, int code)
{
	ofstream outp;
	int i, j;
	int counter = 0;

	int n_i = sdg.getN_i();
	int n_j = sdg.getN_j();
	outp.open(fname);
	int nsd = sdg.getNoSpaceDim();

	// code==4 if we don't have observation points: we write just 0
	// for example, we use 0 observation points from file Obs_points_fdm.m
	if (code == 4)
	{

		outp << 0 << endl;
	}

	// code==20 for some points in fem-domain
	// for example, we use 8 observation points from file Obs_points_left3d.m
	if (code == 20)
	{
		cout << endl << " works code==20" << endl;
		Mat_real CoordArray;

		CoordArray.newsize(38, 3);

		ReadSolfromFile((char *) "Obs_back.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

	// code==19 for some points in fem-domain
	// for example, we use 8 observation points from file Obs_points_left3d.m
	if (code == 19)
	{
		cout << endl << " works code==19" << endl;
		Mat_real CoordArray;

		CoordArray.newsize(38, 3);

		// ReadSolfromFile("Obs_plane_front3d.m",CoordArray);
		ReadSolfromFile((char *) "Obs_front.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

	// code==18 for some points in fem-domain
	// for example, we use 8 observation points from file Obs_points_left3d.m
	if (code == 18)
	{
		cout << endl << " works code==18" << endl;
		Mat_real CoordArray;
		CoordArray.newsize(18, 3);

		ReadSolfromFile((char *) "Obs_bot.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{

			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

	// code==17 for some points in fem-domain

	if (code == 17)
	{
		cout << endl << " works code==17" << endl;
		Mat_real CoordArray;
		CoordArray.newsize(18, 3);

		ReadSolfromFile((char *) "Obs_top.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

	// code==16 for some points in fem-domain
	// for example, we use 8 observation points from file Obs_points_left3d.m
	if (code == 16)
	{
		cout << endl << " works code==16" << endl;
		Mat_real CoordArray;

		CoordArray.newsize(38, 3);
		ReadSolfromFile((char *) "Obs_left.m", CoordArray);

		//for grid_ubot.dat

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

// code==15 for some points in fem-domain
	// for example, we use 15 observation points from file Obs_pointsl2.m
	if (code == 15)
	{
		int node;
		Mat_real CoordArray;

		//for grid_ubotref1.dat
		CoordArray.newsize(38, 3);
		ReadSolfromFile((char *) "Obs_right.m", CoordArray);
		//for grid_ubot.dat
		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << "coordinates:" << CoordArray(j, 0) << " " << CoordArray(j, 1) << " " << CoordArray(j, 2) << endl;
		}

	}

	outp.close();

	//print("Solution at one point is written to file: %s\n", fname);

}

void WavesOutputs::WriteToFile2(char *fname, WavesSDGeometry& sdg, WavesGridB& gg, int code)
{
	ofstream outp;
	int i, j;
	int counter = 0;

	int n_i = sdg.getN_i();
	int n_j = sdg.getN_j();
	outp.open(fname);
	int nsd = sdg.getNoSpaceDim();

	// code==4 if we don't have observation points: we write just 0
	// for example, we use 0 observation points from file Obs_points_fdm.m
	if (code == 4)
	{

		outp << 0 << endl;
	}

	// code==20 for some points in fem-domain
	// for example, we use 8 observation points from file Obs_points_left3d.m
	if (code == 20)
	{
		cout << endl << " works code==20" << endl;
		Mat_real CoordArray;
		CoordArray.newsize(38, 3);

		ReadSolfromFile((char *) "Obs_back.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

	// code==19 for some points in fem-domain
	// for example, we use 8 observation points from file Obs_points_left3d.m
	if (code == 19)
	{
		cout << endl << " works code==19" << endl;
		Mat_real CoordArray;

		CoordArray.newsize(38, 3);

		ReadSolfromFile((char *) "Obs_front.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

	// code==18 for some points in fem-domain
	// for example, we use 8 observation points from file Obs_points_left3d.m
	if (code == 18)
	{
		cout << endl << " works code==18" << endl;
		Mat_real CoordArray;
		CoordArray.newsize(10, 3);

		ReadSolfromFile((char *) "Obs_bot2.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

	// code==17 for some points in fem-domain

	if (code == 17)
	{
		cout << endl << " works code==17" << endl;
		Mat_real CoordArray;
		CoordArray.newsize(10, 3);

		ReadSolfromFile((char *) "Obs_top2.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

	// code==16 for some points in fem-domain
	// for example, we use 8 observation points from file Obs_points_left3d.m
	if (code == 16)
	{
		cout << endl << " works code==16" << endl;
		Mat_real CoordArray;

		CoordArray.newsize(20, 3);
		ReadSolfromFile((char *) "Obs_left2.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

// code==15 for some points in fem-domain
	// for example, we use 15 observation points from file Obs_pointsl2.m
	if (code == 15)
	{
		int node;
		Mat_real CoordArray;

		//for grid_ubotref1.dat
		CoordArray.newsize(20, 3);
		ReadSolfromFile((char *) "Obs_right2.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << "coordinates:" << CoordArray(j, 0) << " " << CoordArray(j, 1) << " " << CoordArray(j, 2) << endl;
		}

	}

	outp.close();

	//print("Solution at one point is written to file: %s\n", fname);

}

void WavesOutputs::WriteToFile002(char *fname, WavesSDGeometry& sdg, WavesGridB& gg, int code)
{
	ofstream outp;
	int i, j;
	int counter = 0;

	int n_i = sdg.getN_i();
	int n_j = sdg.getN_j();
	outp.open(fname);
	int nsd = sdg.getNoSpaceDim();

	// code==4 if we don't have observation points: we write just 0
	// for example, we use 0 observation points from file Obs_points_fdm.m
	if (code == 4)
	{

		outp << 0 << endl;
	}

	// code==20 for some points in fem-domain
	// for example, we use 8 observation points from file Obs_points_left3d.m
	if (code == 20)
	{
		cout << endl << " works code==20" << endl;
		Mat_real CoordArray;
		CoordArray.newsize(38, 3);

		ReadSolfromFile((char *) "Obs_back.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

	// code==19 for some points in fem-domain
	// for example, we use 8 observation points from file Obs_points_left3d.m
	if (code == 19)
	{
		cout << endl << " works code==19" << endl;
		Mat_real CoordArray;
		CoordArray.newsize(38, 3);

		ReadSolfromFile((char *) "Obs_front.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

	// code==18 for some points in fem-domain
	// for example, we use 8 observation points from file Obs_points_left3d.m
	if (code == 18)
	{
		cout << endl << " works code==18" << endl;
		Mat_real CoordArray;
		CoordArray.newsize(10, 3);

		ReadSolfromFile((char *) "Obs_bot2.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

	// code==17 for some points in fem-domain

	if (code == 17)
	{
		cout << endl << " works code==17" << endl;
		Mat_real CoordArray;
		CoordArray.newsize(10, 3);

		ReadSolfromFile((char *) "Obs_top2.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

	// code==16 for some points in fem-domain
	// for example, we use 8 observation points from file Obs_points_left3d.m
	if (code == 16)
	{
		cout << endl << " works code==16" << endl;
		Mat_real CoordArray;

		CoordArray.newsize(22, 3);
		ReadSolfromFile((char *) "Obs_left002.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

// code==15 for some points in fem-domain
	// for example, we use 15 observation points from file Obs_pointsl2.m
	if (code == 15)
	{
		int node;
		Mat_real CoordArray;

		//for grid_ubotref1.dat
		CoordArray.newsize(22, 3);
		ReadSolfromFile((char *) "Obs_right002.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << "coordinates:" << CoordArray(j, 0) << " " << CoordArray(j, 1) << " " << CoordArray(j, 2) << endl;
		}

	}

	outp.close();

	//print("Solution at one point is written to file: %s\n", fname);

}

void WavesOutputs::WriteToFile0025(char *fname, WavesSDGeometry& sdg, WavesGridB& gg, int code)
{
	ofstream outp;
	int i, j;
	int counter = 0;

	int n_i = sdg.getN_i();
	int n_j = sdg.getN_j();
	outp.open(fname);
	int nsd = sdg.getNoSpaceDim();

	// code==4 if we don't have observation points: we write just 0
	// for example, we use 0 observation points from file Obs_points_fdm.m
	if (code == 4)
	{

		outp << 0 << endl;
	}

	// code==20 for some points in fem-domain
	// for example, we use 8 observation points from file Obs_points_left3d.m
	if (code == 20)
	{
		cout << endl << " works code==20" << endl;
		Mat_real CoordArray;
		CoordArray.newsize(38, 3);

		ReadSolfromFile((char *) "Obs_back.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

	// code==19 for some points in fem-domain
	// for example, we use 8 observation points from file Obs_points_left3d.m
	if (code == 19)
	{
		cout << endl << " works code==19" << endl;
		Mat_real CoordArray;
		CoordArray.newsize(38, 3);

		ReadSolfromFile((char *) "Obs_front.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

	// code==18 for some points in fem-domain
	// for example, we use 8 observation points from file Obs_points_left3d.m
	if (code == 18)
	{
		cout << endl << " works code==18" << endl;
		Mat_real CoordArray;
		CoordArray.newsize(10, 3);

		ReadSolfromFile((char *) "Obs_bot2.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

	// code==17 for some points in fem-domain

	if (code == 17)
	{
		cout << endl << " works code==17" << endl;
		Mat_real CoordArray;
		CoordArray.newsize(10, 3);

		ReadSolfromFile((char *) "Obs_top2.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

	// code==16 for some points in fem-domain
	// for example, we use 8 observation points from file Obs_points_left3d.m
	if (code == 16)
	{
		cout << endl << " works code==16" << endl;
		Mat_real CoordArray;

		CoordArray.newsize(46, 3);
		ReadSolfromFile((char *) "Obs_left0025.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			int node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << endl;
		}

	}

// code==15 for some points in fem-domain
	// for example, we use 15 observation points from file Obs_pointsl2.m
	if (code == 15)
	{
		int node;
		Mat_real CoordArray;

		//for grid_ubotref1.dat
		CoordArray.newsize(46, 3);
		ReadSolfromFile((char *) "Obs_right0025.m", CoordArray);

		outp << CoordArray.size(0) << endl;

		cout << endl << CoordArray.size(0) << "  " << CoordArray.size(1) << endl;
		for (j = 0; j < CoordArray.size(0); j++)
		{
			node = sdg.coord2node(CoordArray(j, 0), CoordArray(j, 1), CoordArray(j, 2));
			outp << node << endl;
			cout << endl << "node " << node << "coordinates:" << CoordArray(j, 0) << " " << CoordArray(j, 1) << " " << CoordArray(j, 2) << endl;
		}

	}

	outp.close();

	//print("Solution at one point is written to file: %s\n", fname);

}

void WavesOutputs::ReadConst(char *fname, int& alfa)
{
	ifstream inp;
	inp.open(fname);
	int i;
	printf("Open the file: %s\n", fname);

	inp >> alfa;
	cout << endl << " const is " << alfa << endl;
	inp.close();

}

int WavesOutputs::offset(int j)
{
	return j * WIDTH;
}

int WavesOutputs::offset(int n_i, int n_j, int iter)
{
	return n_i * n_j * WIDTH * iter;
}

void WavesOutputs::Write_const(char *fname, double& alfa)
{
	ofstream outp;

	outp.open(fname);

	printf("Open the file: %s\n", fname);
	outp << alfa;

	outp.close();

}

void WavesOutputs::Write_const(char *fname, int alfa)
{
	ofstream outp;

	outp.open(fname);

	printf("Open the file: %s\n", fname);
	outp << alfa;

	outp.close();

}
void WavesOutputs::Write_my_part(char* file, double& value, int iter)
{

	int i, j;

	FILE *fp;
	fp = fopen(file, "a+");

	if (!fp)
		cout << endl << "WavesOutputs: cannot open file \"" << file << "\"\n";
	cout << endl << "opened file  " << file << endl;

	fseek(fp, offset(iter), SEEK_SET);

	if (iter == 1)
		fprintf(fp, "Started new computations:  %g ", value);
	else
		fprintf(fp, "%g ", value);

	cout << endl << "closed file" << file << endl;
	fclose(fp);
}

void WavesOutputs::Write_my_part(char* file)
{
	int i, j;

	FILE *fp;
	fp = fopen(file, "a+");

	if (!fp)
		cout << endl << "WavesOutputs: cannot open file \"" << file << "\"\n";
	cout << endl << "opened file  " << file << endl;

	fseek(fp, offset(it), SEEK_SET);

	fprintf(fp, "%g ", norma);

	cout << endl << "closed file" << file << endl;
	fclose(fp);
}

void WavesOutputs::Write_array_my_part(char* file, Mat_real& array_part, int iter)
{

	int i, j, n_i, n_j;
	n_i = array_part.size(0);
	n_j = array_part.size(1);

	FILE *fp;
	fp = fopen(file, "a+");

	if (!fp)
		cout << endl << "WavesOutputs: cannot open file \"" << file << "\"\n";
	//cout << endl << "opened file  " << file << " nnodes n_i = " << n_i << " nrtimesteps  n_j = " << n_j << endl;

	fseek(fp, offset(n_i, n_j, iter), SEEK_SET);

	for (j = 0; j < n_j; j++)
		for (i = 0; i < n_i; i++)
			fprintf(fp, "%15.8f ", array_part(i, j));

	//cout << endl << "closed file" << file << endl;
	fclose(fp);
}

void WavesOutputs::Write_vector_my_part(char* file, Mat_real& array_part, int iter)
{

	int i, j, n_i, n_j;
	n_i = array_part.size(0);

	FILE *fp;
	fp = fopen(file, "a+");

	if (!fp)
		cout << endl << "WavesOutputs: cannot open file \"" << file << "\"\n";
	//cout << endl << "opened file  " << file << " nnodes n_i = " << n_i << endl;

	fseek(fp, 1, SEEK_END);

	for (i = 0; i < n_i; i++)
		fprintf(fp, "%15.8f ", array_part(i, 1));

	fprintf(fp, "\n");

	//cout << endl << "closed file" << file << endl;
	fclose(fp);
}

void WavesOutputs::Write_array_my_part(char* file, MV_Vector<double>& array_part, int iter)
{

	int i, j, n_i, n_j;
	n_i = array_part.size();

	FILE *fp;
	fp = fopen(file, "a+");

	//cout << endl << "opened file  " << file << " nnodes n_i = " << n_i << endl;

	if (!fp)
		cout << endl << "WavesOutputs: cannot open file \"" << file << "\"\n";

	fseek(fp, 1, SEEK_END);

	for (i = 0; i < n_i; i++)
		fprintf(fp, "%20.15lf ", array_part(i));

	fprintf(fp, "\n");

	//cout << endl << "closed file" << file << endl;
	fclose(fp);
}




void WavesOutputs::Write_array_my_part(char* file, MV_Vector<int> array_part) // added by Thanh
{

	int i, j, n_i, n_j;
	n_i = array_part.size();

	FILE *fp;
	fp = fopen(file, "a+");

	//cout << endl << "opened file  " << file << " nnodes n_i = " << n_i << endl;

	if (!fp)
		cout << endl << "WavesOutputs: cannot open file \"" << file << "\"\n";

	fseek(fp, 1, SEEK_END);

	for (i = 0; i < n_i; i++)
		fprintf(fp, "%i  ", array_part(i));

	fprintf(fp, "\n");

	//cout << endl << "closed file" << file << endl;
	fclose(fp);
}

void WavesOutputs::Write_matrix_to_file(char* file, Mat_real& array_part) // write a matrix to a file. Thanh added
{

	int i, j, n_i, n_j;
	n_i = array_part.size(0);
	n_j = array_part.size(1);

	ofstream outp;
	outp.open(file);

	//cout << endl << "opened file  " << file << " No of row n_i = " << n_i << " No of column n_j = " << n_j << endl;

	for (i = 0; i < n_i; i++)
	  {
		for (j = 0; j < n_j; j++)
		  {
			outp << array_part(i,j) << "  ";
		  }
		outp << endl;
	  }

	//cout << endl << "closed file" << file << endl;
	outp.close();
}


double WavesOutputs::L2Norma()
{
	int i, n, el;
	int ierr;
	int nsd = gg.getNoSpaceDim();
	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();
	int n_1, n_2, n_3, n_4;

	cout << endl << "nsd = " << nsd << " nno " << nno << endl;

	double sol = 0.0;

	if (nsd == 2)
	{
		FET3n2D F(&gg);
		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTRI1)
			{
				F.refill(el);
				real volume = F.area();
				n_1 = F.n1();
				n_2 = F.n2();
				n_3 = F.n3();
				double val1 = bb(n_1) * bb(n_1);
				double val2 = bb(n_2) * bb(n_2);
				double val3 = bb(n_3) * bb(n_3);

				sol += ((val1 + val2 + val3) / 3.0) * volume;
				cout << endl << " element " << el << "bb(n1)*bb(n1) " << val1 << "bb(n2)*bb(n2) " << val2 << "bb(n3)*bb(n3) " << val3 << "  area " << volume << endl;

			}
		}
	}
	else if (nsd == 3)
	{
		FET4n3D F(&gg);
		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTET1)
			{
				F.refill(el);
				real volume = fabs(F.volume());
				n_1 = F.n1();
				n_2 = F.n2();
				n_3 = F.n3();
				n_4 = F.n4();

				volume *= 0.25;

				double val1 = bb(n_1) * bb(n_1);
				double val2 = bb(n_2) * bb(n_2);
				double val3 = bb(n_3) * bb(n_3);
				double val4 = bb(n_4) * bb(n_4);
				sol += (val1 + val2 + val3 + val4) * volume;

			}
		}
	}

	sol = sqrt(sol);
	cout << endl << " norma is  " << sol << endl;
	return (sol);

}

double WavesOutputs::L2Norma(MV_Vector<int>& glob_nodes, MV_Vector<double>& values_in_glob_nodes)
{
	int i, n, el;
	int ierr;
	int nsd = gg.getNoSpaceDim();
	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();
	int n_1, n_2, n_3, n_4;

	for (i = 0; i < glob_nodes.size(); i++)
	{
		n = glob_nodes(i);
		bb(n) = values_in_glob_nodes(i);
	}

	cout << endl << "nsd = " << nsd << " nno " << nno << endl;

	double sol = 0.0;

	if (nsd == 2)
	{
		FET3n2D F(&gg);
		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTRI1)
			{
				F.refill(el);
				real volume = F.area();
				n_1 = F.n1();
				n_2 = F.n2();
				n_3 = F.n3();
				double val1 = bb(n_1) * bb(n_1);
				double val2 = bb(n_2) * bb(n_2);
				double val3 = bb(n_3) * bb(n_3);
				sol += ((val1 + val2 + val3) / 3.0) * volume;

			}
		}
	}
	else if (nsd == 3)
	{
		FET4n3D F(&gg);
		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTET1)
			{
				F.refill(el);
				real volume = fabs(F.volume());
				n_1 = F.n1();
				n_2 = F.n2();
				n_3 = F.n3();
				n_4 = F.n4();

				volume *= 0.25;

				double val1 = bb(n_1) * bb(n_1);
				double val2 = bb(n_2) * bb(n_2);
				double val3 = bb(n_3) * bb(n_3);
				double val4 = bb(n_4) * bb(n_4);
				sol += (val1 + val2 + val3 + val4) * volume;

			}
		}
	}

	sol = sqrt(sol);
	cout << endl << " norma is  " << sol << endl;
	return (sol);

}

// configuration functions for loading input parameters: added by Thanh

// -- load the input parameters from the .dat file:
void WavesOutputs::Configurate(char *fname, Grid& grid, Grid& outer_gg, 
		 WavesSDGeometry& sdg, WavesSDGeometry& sdg_new, int& nsd, 
		 bool& USE_EXCHANGE, bool& USE_FEM, bool& USE_FDM, bool &USE_DIR_FEM, 
		 bool &USE_DIR_FDM, bool& USE_ABSORB, bool &USE_RHS,
		 bool &PRINT_FILES, int &NoTimeSteps, double &maxTime, double &rhs_code, 
		 int& NoObjects, int &TypeOfMat1, int &TypeOfMat2, int &TypeOfMat3, 
		 double& velocity1, double& velocity2, double& velocity3,
		 double& noiselevel, double& freq, MV_Vector<double>& Zposition)

{
  const int maxlen = 100;
  double dxFDM, dyFDM, dzFDM, x_minFDM, y_minFDM, z_minFDM;
  double dxFEM, dyFEM, dzFEM, x_minFEM, y_minFEM, z_minFEM;
  int NxFDM, NyFDM, NzFDM, NxFEM, NyFEM, NzFEM;
  char text[maxlen], 
       FDMGridFile[maxlen],	
       FEMGridFile[maxlen],
       BoundaryGridFile[maxlen];


  ifstream inpfile;
  inpfile.open(fname);
if (inpfile.is_open())
{
  inpfile >> text >> nsd; //number of space dimensions
  inpfile >> text >> FDMGridFile;		// FDM grid file with extension *.inp
  inpfile >> text >> FEMGridFile;		// FEM grid file with extension *.inp
  inpfile >> text >> BoundaryGridFile;		// file with extension *.inp for boundary nodes in the FEM domain
  inpfile >> text >> NxFDM >> NyFDM >> NzFDM >> dxFDM >> dyFDM >> dzFDM
	  >> x_minFDM >> y_minFDM >> z_minFDM; //FDM domain
  inpfile >> text >> NxFEM >> NyFEM >> NzFEM >> dxFEM >> dyFEM >> dzFEM 
	  >> x_minFEM >> y_minFEM >> z_minFEM; //FEM domain
  inpfile >> text >> USE_EXCHANGE; // use the hybrid FDM-FEM method  
  inpfile >> text >> USE_FEM; 
  inpfile >> text >> USE_FDM; 
  inpfile >> text >> USE_RHS; // use the right hand side (source term)
  inpfile >> text >> USE_DIR_FEM; 
  inpfile >> text >> USE_DIR_FDM; 
  inpfile >> text >> USE_ABSORB; // use the absorbing b.c.
  inpfile >> text >> PRINT_FILES;
  inpfile >> text >> NoTimeSteps;
  inpfile >> text >> maxTime; 
  inpfile >> text >> rhs_code; //code for different right hand side functions
  inpfile >> text >> NoObjects; //number of objects 

  double test; // unused variable!! 

  TypeOfMat1 = 0; TypeOfMat2 = 0; TypeOfMat3 = 0;
  velocity1 = 0; velocity2 = 0; velocity3 = 0;

  if (NoObjects==1)
    {
      inpfile >> text >> TypeOfMat1;
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> velocity1;
    }
else if (NoObjects==2)
    {
      inpfile >> text >> TypeOfMat1 >> TypeOfMat2;
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> velocity1  >> velocity2;
   
    }
 else if  (NoObjects==3)
   {
      inpfile >> text >> TypeOfMat1 >> TypeOfMat2 >> TypeOfMat3;
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> velocity1  >> velocity2 >> velocity3; 
  
    }
 inpfile >> text >> noiselevel; // noise level
 inpfile >> text >> freq; // frequency of the incident plane wave 
 inpfile >> text >> Zposition(0) >> Zposition(1) >> Zposition(2) >> Zposition(3); // values of z at which the data is saved

  inpfile.close();
  int dbg = 1;
	if (dbg)
		cout << endl << "Dimensions:" << nsd << ", max_time   " << maxTime << endl;

	grid.scan(FEMGridFile);
	outer_gg.scan(FDMGridFile);
	sdg.initialize(NxFDM, NyFDM, NzFDM, nsd, dxFDM, dyFDM, dzFDM, x_minFDM, y_minFDM, z_minFDM);
	sdg_new.initialize(NxFEM, NyFEM, NzFEM, nsd, dxFEM, dyFEM, dzFEM, x_minFEM, y_minFEM, z_minFEM); 

}
else 
  cout << endl << endl << "Error in WavesOutputs::Configurate: file cannot be opened" << endl << endl;

}


// -- load the input parameters from the .dat file:
void WavesOutputs::Configurate(char *fname, Grid& grid, Grid& outer_gg, 
		 WavesSDGeometry& sdg, WavesSDGeometry& sdg_new, int& nsd, 
		 bool& USE_EXCHANGE, bool& USE_FEM, bool& USE_FDM, bool &USE_DIR_FEM, 
		 bool &USE_DIR_FDM)

{
  const int maxlen = 100;
  double dxFDM, dyFDM, dzFDM, x_minFDM, y_minFDM, z_minFDM;
  double dxFEM, dyFEM, dzFEM, x_minFEM, y_minFEM, z_minFEM;
  int NxFDM, NyFDM, NzFDM, NxFEM, NyFEM, NzFEM;
  bool USE_RHS; 	
  char text[maxlen], 
       FDMGridFile[maxlen],	
       FEMGridFile[maxlen],
       BoundaryGridFile[maxlen];


  ifstream inpfile;
  inpfile.open(fname);
if (inpfile.is_open())
{
  inpfile >> text >> nsd; //number of space dimensions
  inpfile >> text >> FDMGridFile;		// FDM grid file with extension *.inp
  inpfile >> text >> FEMGridFile;		// FEM grid file with extension *.inp
  inpfile >> text >> BoundaryGridFile;		// file with extension *.inp for boundary nodes in the FEM domain
  inpfile >> text >> NxFDM >> NyFDM >> NzFDM >> dxFDM >> dyFDM >> dzFDM
	  >> x_minFDM >> y_minFDM >> z_minFDM; //FDM domain
  inpfile >> text >> NxFEM >> NyFEM >> NzFEM >> dxFEM >> dyFEM >> dzFEM 
	  >> x_minFEM >> y_minFEM >> z_minFEM; //FEM domain
  inpfile >> text >> USE_EXCHANGE; // use the hybrid FDM-FEM method  
  inpfile >> text >> USE_FEM; 
  inpfile >> text >> USE_FDM; 
  inpfile >> text >> USE_RHS; // use the right hand side (source term)
  inpfile >> text >> USE_DIR_FEM; 
  inpfile >> text >> USE_DIR_FDM; 
 
  inpfile.close();

	grid.scan(FEMGridFile);
	outer_gg.scan(FDMGridFile);
	sdg.initialize(NxFDM, NyFDM, NzFDM, nsd, dxFDM, dyFDM, dzFDM, x_minFDM, y_minFDM, z_minFDM);
	sdg_new.initialize(NxFEM, NyFEM, NzFEM, nsd, dxFEM, dyFEM, dzFEM, x_minFEM, y_minFEM, z_minFEM); 

}
else 
  cout << endl << endl << "Error in WavesOutputs::Configurate: file cannot be opened" << endl << endl;

}



// -- load the input parameters from the .dat file:
void WavesOutputs::Configurate(char *fname, Grid& grid, Grid& outer_gg, 
		 	 	WavesSDGeometry& sdg, WavesSDGeometry& sdg_new, int& nsd, 
			 	bool& USE_EXCHANGE, bool& USE_FEM, bool& USE_FDM, bool &USE_DIR_FEM, 
			 	bool &USE_DIR_FDM, bool& USE_ABSORB,
			 	int &NoTimeSteps, double &maxTime, double& freq)

{

  bool USE_RHS, PRINT_FILES;
  double rhs_code;
  int NoObjects, TypeOfMat1, TypeOfMat2, TypeOfMat3;
  double velocity1, velocity2, velocity3, noiselevel;
 
  const int maxlen = 100;
  double dxFDM, dyFDM, dzFDM, x_minFDM, y_minFDM, z_minFDM;
  double dxFEM, dyFEM, dzFEM, x_minFEM, y_minFEM, z_minFEM;
  int NxFDM, NyFDM, NzFDM, NxFEM, NyFEM, NzFEM;
  char text[maxlen], 
       FDMGridFile[maxlen],	
       FEMGridFile[maxlen],
       BoundaryGridFile[maxlen];


  ifstream inpfile;
  inpfile.open(fname);
if (inpfile.is_open())
{
  inpfile >> text >> nsd; //number of space dimensions
  inpfile >> text >> FDMGridFile;		// FDM grid file with extension *.inp
  inpfile >> text >> FEMGridFile;		// FEM grid file with extension *.inp
  inpfile >> text >> BoundaryGridFile;		// file with extension *.inp for boundary nodes in the FEM domain
  inpfile >> text >> NxFDM >> NyFDM >> NzFDM >> dxFDM >> dyFDM >> dzFDM
	  >> x_minFDM >> y_minFDM >> z_minFDM; //FDM domain
  inpfile >> text >> NxFEM >> NyFEM >> NzFEM >> dxFEM >> dyFEM >> dzFEM 
	  >> x_minFEM >> y_minFEM >> z_minFEM; //FEM domain
  inpfile >> text >> USE_EXCHANGE; // use the hybrid FDM-FEM method  
  inpfile >> text >> USE_FEM; 
  inpfile >> text >> USE_FDM; 
  inpfile >> text >> USE_RHS; // use the right hand side (source term)
  inpfile >> text >> USE_DIR_FEM; 
  inpfile >> text >> USE_DIR_FDM; 
  inpfile >> text >> USE_ABSORB; // use the absorbing b.c.
  inpfile >> text >> PRINT_FILES;
  inpfile >> text >> NoTimeSteps;
  inpfile >> text >> maxTime; 
  inpfile >> text >> rhs_code; //code for different right hand side functions
  inpfile >> text >> NoObjects; //number of objects 

  double test; // unused variable!! 

  TypeOfMat1 = 0; TypeOfMat2 = 0; TypeOfMat3 = 0;
  velocity1 = 0; velocity2 = 0; velocity3 = 0;

  if (NoObjects==1)
    {
      inpfile >> text >> TypeOfMat1;
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> velocity1;
    }
else if (NoObjects==2)
    {
      inpfile >> text >> TypeOfMat1 >> TypeOfMat2;
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> velocity1  >> velocity2;
   
    }
 else if  (NoObjects==3)
   {
      inpfile >> text >> TypeOfMat1 >> TypeOfMat2 >> TypeOfMat3;
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> velocity1  >> velocity2 >> velocity3; 
  
    }
 inpfile >> text >> noiselevel; // noise level
 inpfile >> text >> freq; // frequency of the incident plane wave 


int dbg = 1;
	if (dbg)
		cout << "Dimensions:" << nsd << ", max_time   " << maxTime << endl;

	grid.scan(FEMGridFile);
	outer_gg.scan(FDMGridFile);
	sdg.initialize(NxFDM, NyFDM, NzFDM, nsd, dxFDM, dyFDM, dzFDM, x_minFDM, y_minFDM, z_minFDM);
	sdg_new.initialize(NxFEM, NyFEM, NzFEM, nsd, dxFEM, dyFEM, dzFEM, x_minFEM, y_minFEM, z_minFEM); 

  inpfile.close();

}
else 
  cout << endl << "Error in WavesOutputs::Configurate: file cannot be opened" << endl << endl;

}





// load the boundary grid file:
string WavesOutputs::Configurate(char *fname)
{
  const int maxlen = 100;
	int nsd;  
  char text[maxlen], 
       FDMGridFile[maxlen],	
       FEMGridFile[maxlen];

  string BoundaryGridFile;
  
  ifstream inpfile;
  inpfile.open(fname);
if (inpfile.is_open())
{
  inpfile >> text >> nsd; //number of space dimensions
  inpfile >> text >> FDMGridFile;		// FDM grid file with extension *.inp
  inpfile >> text >> FEMGridFile;		// FEM grid file with extension *.inp
  inpfile >> text >> BoundaryGridFile;		// file with extension *.inp for boundary nodes in the FEM domain
 
  cout << endl << BoundaryGridFile << endl;
  //const char* filename = BoundaryGridFile.c_str();

}
else
{
  cout << endl << endl << "Error in WavesOutputs::Configurate: file cannot be opened" << endl << endl;
  BoundaryGridFile = " "; 
}	
  return BoundaryGridFile; 
}



//extract the data file names:
string WavesOutputs::Configurate(char *fname, int idx)
{
  const int maxlen = 100;
  double x;
  int Nx, NoObjects, i;
  char text[maxlen], FDMGridFile[maxlen], FEMGridFile[maxlen];
  
  string BoundaryGridFile, File1, File2, File3, File4; 

  ifstream inpfile;
  inpfile.open(fname);
if (inpfile.is_open())
{
  inpfile >> text >> Nx; //number of space dimensions
  inpfile >> text >> FDMGridFile;		// FDM grid file with extension *.inp
  inpfile >> text >> FEMGridFile;		// FEM grid file with extension *.inp
  inpfile >> text >> BoundaryGridFile;		// file with extension *.inp for boundary nodes in the FEM domain
  inpfile >> text >> Nx >> Nx >> Nx >> x >> x >> x >> x >> x >> x; //FDM domain
  inpfile >> text >> Nx >> Nx >> Nx >> x >> x >> x >> x >> x >> x; //FDM domain
  inpfile >> text >> Nx; // use the hybrid FDM-FEM method  
  inpfile >> text >> Nx; 
  inpfile >> text >> Nx; 
  inpfile >> text >> Nx; // use the right hand side (source term)
  inpfile >> text >> Nx; 
  inpfile >> text >> Nx; 
  inpfile >> text >> Nx; // use the absorbing b.c.
  inpfile >> text >> Nx;
  inpfile >> text >> Nx;
  inpfile >> text >> x; 
  inpfile >> text >> x; //code for different right hand side functions
  inpfile >> text >> NoObjects; //number of objects 

  double test; // unused variable!! 
  double TypeOfMat, velocity;

  if (NoObjects==1)
    {
      inpfile >> text >> TypeOfMat;
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> velocity;
    }
else if (NoObjects==2)
    {
      inpfile >> text >> TypeOfMat >> TypeOfMat;
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> velocity  >> velocity;
   
    }
 else if  (NoObjects==3)
   {
      inpfile >> text >> TypeOfMat >> TypeOfMat >> TypeOfMat;
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> test >> test >> test >> test >> test >> test; // coordinate of the subdomain, not used!
      inpfile >> text >> velocity  >> velocity >> velocity; 
  
    }
 inpfile >> text >> x; // noise level
 inpfile >> text >> Nx; // frequency of the incident plane wave 
 inpfile >> text >> x >> x >> x >> x; // values of z at which the data is saved
 inpfile >> text >> File1 >> File2 >> File3 >> File4; 

  if (idx == 1)
  {	
	return File1;
  }
  else if (idx == 2)
  {	
	return File2;
  }
  else if (idx == 3)
  {	
	return File3;
  }
  else 
  {	
	return File4;
  }
}
else 
{
  cout << endl << endl << "Error in WavesOutputs::Configurate: file cannot be opened" << endl << endl;
  return " ";
}

}


//===load the values of vector b for 1 inclusion

void WavesOutputs::VecGetTypeofMat(WavesGridB& gg, int nsd, int type_of_material, double velocity, MV_Vector<double>& b)
{
	// type of material for element is type of mat in file -1
	// example: in file difmat.inp ellipse have tofmat = 2
	// then gg.getMaterialType(ellipse) = 2 -1 = 1
	// and 1 we should put in pulse2D.dat file
	int nel = gg.getNoElms();
	int nno = gg.getNoNodes();

	int n1, n2, n3, n4, el, i;
	MV_Vector<int> Markers(nno);
	Markers = 0;
	double x, y;
	x = 0.0;
	y = 0.0;

	for (el = 0; el < nel; el++)
	{
		// for inner small ellipse and for outer ellipse
		if (gg.getMaterialType(el) == type_of_material)
		{
			if (nsd == 2)
			{
				n1 = gg.loc2glob(el, 0);
				n2 = gg.loc2glob(el, 1);
				n3 = gg.loc2glob(el, 2);

				cout << endl << "n1 " << n1 << " n2 " << n2 << " n3 " << n3 << endl;
				//to avoid repeating multiplications
				Markers(n1) = 1;
				Markers(n2) = 1;
				Markers(n3) = 1;
			}
			else if (nsd == 3)
			{
				n1 = gg.loc2glob(el, 0);
				n2 = gg.loc2glob(el, 1);
				n3 = gg.loc2glob(el, 2);
				n4 = gg.loc2glob(el, 3);

				//to avoid repeating multiplications
				Markers(n1) = 1;
				Markers(n2) = 1;
				Markers(n3) = 1;
				Markers(n4) = 1;

			}
		}
	}

	for (i = 0; i < gg.getNoNodes(); i++)
		if (Markers(i) == 1)
		{
			b(i) = velocity;
		}

}
//===== 2 types of materials:
void WavesOutputs::VecGetTypeofMat(WavesGridB& gg, int nsd,
		     int type_of_material1, int type_of_material2,
		     double velocity1, double velocity2, 
		     MV_Vector<double>& b)
{
	// type of material for element is type of mat in file -1
	// example: in file difmat.inp ellipse have tofmat = 2
	// then gg.getMaterialType(ellipse) = 2 -1 = 1
	// and 1 we should put in pulse2D.dat file
	int nel = gg.getNoElms();
	int nno = gg.getNoNodes();

	int n1, n2, n3, n4, el, i;
	MV_Vector<int> Markers(nno);
	Markers = 0;
	double x, y;
	x = 0.0;
	y = 0.0;

	for (el = 0; el < nel; el++)
	{
	
			if (nsd == 2)
			{
				n1 = gg.loc2glob(el, 0);
				n2 = gg.loc2glob(el, 1);
				n3 = gg.loc2glob(el, 2);

				cout << endl << "n1 " << n1 << " n2 " << n2 << " n3 " << n3 << endl;
				//to avoid repeating multiplications
				if (gg.getMaterialType(el) == type_of_material1)
				  {				  
				    Markers(n1) = 1;
				    Markers(n2) = 1;
				    Markers(n3) = 1;
				  }
				else if  (gg.getMaterialType(el) == type_of_material2)
				  {
				    Markers(n1) = 2;
				    Markers(n2) = 2;
				    Markers(n3) = 2;
				  }
			}
			else if (nsd == 3)
			  {
			    n1 = gg.loc2glob(el, 0);
			    n2 = gg.loc2glob(el, 1);
			    n3 = gg.loc2glob(el, 2);
			    n4 = gg.loc2glob(el, 3);
			    
			    //to avoid repeating multiplications
			    if (gg.getMaterialType(el) == type_of_material1)
			      {
				Markers(n1) = 1;
				Markers(n2) = 1;
				Markers(n3) = 1;
				Markers(n4) = 1;
			      }
			    else  if (gg.getMaterialType(el) == type_of_material2)
			      {
				Markers(n1) = 2;
				Markers(n2) = 2;
				Markers(n3) = 2;
				Markers(n4) = 2;
				
			      }
			  } // for nsd==3
	} // for el

	for (i = 0; i < gg.getNoNodes(); i++)
	  {
		if (Markers(i) == 1)
		  {
		    b(i) = velocity1;
		  }
		else if (Markers(i) == 2)
		  {
			b(i) = velocity2;
		}
	  }
}

//===== 3 types of materials:
void WavesOutputs::VecGetTypeofMat(WavesGridB& gg, int nsd,
		     int type_of_material1, int type_of_material2, int type_of_material3,
		     double velocity1, double velocity2, double velocity3,
		     MV_Vector<double>& b)
{
	// type of material for element is type of mat in file -1
	// example: in file difmat.inp ellipse have tofmat = 2
	// then gg.getMaterialType(ellipse) = 2 -1 = 1
	// and 1 we should put in pulse2D.dat file
	int nel = gg.getNoElms();
	int nno = gg.getNoNodes();

	int n1, n2, n3, n4, el, i;
	MV_Vector<int> Markers(nno);
	Markers = 0;
	double x, y;
	x = 0.0;
	y = 0.0;

	for (el = 0; el < nel; el++)
	{
	
			if (nsd == 2)
			{
				n1 = gg.loc2glob(el, 0);
				n2 = gg.loc2glob(el, 1);
				n3 = gg.loc2glob(el, 2);

				cout << endl << "n1 " << n1 << " n2 " << n2 << " n3 " << n3 << endl;
				//to avoid repeating multiplications
				if (gg.getMaterialType(el) == type_of_material1)
				  {				  
				    Markers(n1) = 1;
				    Markers(n2) = 1;
				    Markers(n3) = 1;
				  }
				else if  (gg.getMaterialType(el) == type_of_material2)
				  {
				    Markers(n1) = 2;
				    Markers(n2) = 2;
				    Markers(n3) = 2;
				  }
				else if (gg.getMaterialType(el) == type_of_material3)
				  {
				    Markers(n1) = 3;
				    Markers(n2) = 3;
				    Markers(n3) = 3;
				  }
			}
			else if (nsd == 3)
			  {
			    n1 = gg.loc2glob(el, 0);
			    n2 = gg.loc2glob(el, 1);
			    n3 = gg.loc2glob(el, 2);
			    n4 = gg.loc2glob(el, 3);
			    
			    //to avoid repeating multiplications
			    if (gg.getMaterialType(el) == type_of_material1)
			      {
				Markers(n1) = 1;
				Markers(n2) = 1;
				Markers(n3) = 1;
				Markers(n4) = 1;
			      }
			    else  if (gg.getMaterialType(el) == type_of_material2)
			      {
				Markers(n1) = 2;
				Markers(n2) = 2;
				Markers(n3) = 2;
				Markers(n4) = 2;
				
			      }
			    else  if (gg.getMaterialType(el) == type_of_material3)
			      {
				Markers(n1) = 3;
				Markers(n2) = 3;
				Markers(n3) = 3;
				Markers(n4) = 3;
				
			      }
			  } // for nsd==3
	} // for el

	for (i = 0; i < gg.getNoNodes(); i++)
	  {
		if (Markers(i) == 1)
		  {
		    b(i) = velocity1;
		  }
		else if (Markers(i) == 2)
		  {
			b(i) = velocity2;
		}
		else if (Markers(i) == 3)
		{
			b(i) = velocity3;
		}
	  }
}

// load the data on the boundary, added by Thanh. Each column is a time-domain waveform
void WavesOutputs::load_boundary_data(char *fname, Mat_real& Ar1, MV_Vector<int>& BoundNodes)
{
	
	int i, j, Nno, Nno2;
	Nno = Ar1.size(1); 
	Nno2 = BoundNodes.size();

	if (Nno != Nno2)
		cout << endl << " The input parameters are not consistent!" << "Test " << endl;
		
	else
	{
		ifstream inp;
		inp.open(fname);
	
		printf("Open the file: %s\n", fname);
		//  Ar1.size(1) - number of nodes
		//  Ar1.size(0) - number of timesteps
		
		// load the node numbers of the boundary nodes
		for (j = 0; j < Nno2; j++)
			inp >> BoundNodes(j);
		
		// load the data at the boundary nodes
		for (i = 0; i < Ar1.size(0); i++)
		{
			for (j = 0; j < Nno; j++)
			{
				inp >> Ar1(i, j);

				if (Ar1(i, j) > 1000)
					cout << endl << " very big number: " << Ar1(i, j) << "in the point " << j << "and time moment" << i << endl;
			}
		}
		inp.close();
	}	
}

// load the boundary grid:

int WavesOutputs::load_boundary_grid(char *bdgrid)
{
	int Nbnodes;
 	ifstream inp;
	printf("Open the grid file: %s\n", bdgrid);
 	inp.open(bdgrid);
	if (inp.is_open())
	{
 		inp >> Nbnodes; // number of boundary nodes
 		inp.close();
	}
	else
	{
		cout << endl << "Error in WavesOutputs::load_boundary_grid: file cannot be opened" << endl;
		Nbnodes = 0;
	}
	return Nbnodes;
}

void WavesOutputs::load_boundary_grid(char *bdgrid, MV_Vector<int>& BoundNodes)
{
 	ifstream inp;
 	inp.open(bdgrid);
	if (inp.is_open())
	{
		int Nbnodes, Nbnodes2, i; 
		Nbnodes2 = BoundNodes.size();

 		inp >> Nbnodes; // number of boundary nodes
		if (Nbnodes != Nbnodes2)
			cout << endl << "Error: the number of grid nodes in the grid file is not equal to the number of element in the input vector";
		else
		{
			double x;
 			for (i=0; i < Nbnodes; i++)
 				inp >> BoundNodes(i) >> x >> x >> x;		
		}
 		inp.close();
	}
	else 
  		cout << endl << "Error in WavesOutputs::load_boundary_grid: file cannot be opened" << endl << endl;

}

//load also the coordinates of the boundary grid nodes:
void WavesOutputs::load_boundary_grid(char *bdgrid, MV_Vector<int>& BoundNodes, MV_Vector<double>& X, MV_Vector<double>& Y, MV_Vector<double>& Z)
{
 	ifstream inp;
 	inp.open(bdgrid);
	if (inp.is_open())
	{
		int Nbnodes, Nbnodes2, i; 
		Nbnodes2 = BoundNodes.size();

 		inp >> Nbnodes; // number of boundary nodes
		if (Nbnodes != Nbnodes2)
			cout << endl << "Error: the number of grid nodes in the grid file is not equal to the number of element in the input vector";
		else
		{
			for (i=0; i < Nbnodes; i++)
 				inp >> BoundNodes(i) >> X(i) >> Y(i) >> Z(i);		
		}
 		inp.close();
	}
	else 
  		cout << endl << "Error in WavesOutputs::load_boundary_grid: file cannot be opened" << endl << endl;


}

// load inversion parameters
string WavesOutputs::load_inversion_parameters(char* fname, double& s_min, double& s_max, int& Ns)
{
  	char text[100]; 
  	string data_file; 

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> data_file;		// measured data file name with extension *.m
  	 	inpfile >> text >> text;		// file name with extension *.m  for true solution of the forward problem
  		inpfile >> text >> text;		// file name with extension *.m for incident wave 
 		inpfile >> text >> s_min >> s_max; 
		inpfile >> text >> Ns; 
		inpfile.close();
	}
  	else
	{
		cout << endl << "Error in WavesOutputs::load_inversion_parameters: file cannot be opened " << endl;
		data_file = " ";
	}
	//const char* datafile = data_file.c_str();
  	return data_file;
}


// load inversion parameters: load the file name with the incident wave
string WavesOutputs::load_inversion_parameters(char* fname)
{
  	char text[100]; 
  	string data_file; 

  	ifstream inpfile;
  	inpfile.open(fname);

 	if (inpfile.is_open())
	{
  		inpfile >> text >> text;		// measured data file name with extension *.m
  	 	inpfile >> text >> text;		// file name with extension *.m  for true solution of the forward problem
  		inpfile >> text >> data_file;		// file name with extension *.m for incident wave 
 		inpfile.close();
	}
  	else
	{
		cout << endl << "Error in WavesOutputs::load_inversion_parameters: file cannot be opened " << endl;
		data_file = " ";
	}
	//const char* datafile = data_file.c_str();
  	return data_file;

}

// load inversion parameters: load the file name with the incident wave
string WavesOutputs::load_inversion_parameters_fp(char* fname)
{
  	char text[100]; 
  	string data_file; 

  	ifstream inpfile;
  	inpfile.open(fname);

 	if (inpfile.is_open())
	{
  		inpfile >> text >> text;		// measured data file name with extension *.m
  	 	inpfile >> text >> data_file;		// file name with extension *.m  for true solution of the forward problem
 	
		inpfile.close();
	}
  	else
	{
		cout << endl << "Error in WavesOutputs::load_inversion_parameters: file cannot be opened " << endl;
		data_file = " ";
	}
	//const char* datafile = data_file.c_str();
  	return data_file;

}

// load inversion parameters: load the file name with the incident wave
string WavesOutputs::load_inversion_parameters_typebndmea(char* fname)
{
  	char text[100]; 
	double x; 
	int N, i, test; 
	string typeofboundarymeasurement = " "; 

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> text;		// measured data file name with extension *.m
 	 	inpfile >> text >> text;		// file name with extension *.m  for true solution of the forward problem
 		inpfile >> text >> text;		// file name with extension *.m for incident wave 
 		inpfile >> text >> x >> x; 
		inpfile >> text >> N; // number of pseudo-frequencies
	
		inpfile >> text; 
		for (i=0; i< N; i++)
			inpfile >> x;//load the Carlemann weighting factors
		inpfile >> text; 
		for (i=0; i< N; i++)
			inpfile >> x;//load the regularization parameters
		inpfile >> text; 
		for (i=0; i< N; i++)
			inpfile >> test;//load the maximum number of iterations
	
		inpfile >> text >> x >> x;
		inpfile >> text >> typeofboundarymeasurement; 

		inpfile.close();
	}
	else 
  		cout << endl << "Error in WavesOutputs::load_boundary_grid: file cannot be opened" << endl << endl;

	return typeofboundarymeasurement;
}

// load inversion parameters: load the text string for choice of first tail
string WavesOutputs::load_inversion_parameters_firsttail(char* fname)
{
  	char text[100]; 
	double x; 
	int N, i, test; 
	string firsttail = " "; 
	

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> text;		// measured data file name with extension *.m
 	 	inpfile >> text >> text;		// file name with extension *.m  for true solution of the forward problem
 		inpfile >> text >> text;		// file name with extension *.m for incident wave 
 		inpfile >> text >> x >> x; 
		inpfile >> text >> N; // number of pseudo-frequencies
	
		inpfile >> text; 
		for (i=0; i< N; i++)
			inpfile >> x;//load the Carlemann weighting factors
		inpfile >> text; 
		for (i=0; i< N; i++)
			inpfile >> x;//load the regularization parameters
		inpfile >> text; 
		for (i=0; i< N; i++)
			inpfile >> test;//load the maximum number of iterations
	
		inpfile >> text >> x >> x;
		inpfile >> text >> text; 
		inpfile >> text >> x >> x >> x >> x >> x >> x; 
		inpfile >> text >> firsttail;
		inpfile.close();

	}
	else 
  		cout << endl << "Error in WavesOutputs::load_boundary_grid: file cannot be opened" << endl << endl;

	return firsttail;
}

void WavesOutputs::load_inversion_parameters(char* fname, MV_Vector<double>& Lambda_CWF, MV_Vector<double>& RegPar, 
					     MV_Vector<int>& MaxIter, double& C_lb, double& C_ub)
{
  	char text[100]; 
	double x; 
	int N, i; 

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> text;		// measured data file name with extension *.m
 	 	inpfile >> text >> text;		// file name with extension *.m  for true solution of the forward problem
 		inpfile >> text >> text;		// file name with extension *.m for incident wave 
 		inpfile >> text >> x >> x; 
		inpfile >> text >> N; // number of pseudo-frequencies
	
		if (N != Lambda_CWF.size())
			cout << endl << "Error: The input parameters are not consistent" << endl;
		else	
		{
			inpfile >> text;
			for (i=0; i< N; i++)
				inpfile >> Lambda_CWF(i);//load the Carlemann weighting factors

			inpfile >> text;
			for (i=0; i< N; i++)
				inpfile >> RegPar(i);//load the regularization parameters

			inpfile >> text;
			for (i=0; i< N; i++)
				inpfile >> MaxIter(i);//load the maximum number of iterations
	
			inpfile >> text >> C_lb >> C_ub; 

		} 
		inpfile.close();
	}
	else 
  		cout << endl << "Error in WavesOutputs::load_boundary_grid: file cannot be opened" << endl << endl;

}


void WavesOutputs::load_inversion_parameters(char* fname, MV_Vector<double>& Lambda_CWF, MV_Vector<double>& RegPar, 
					     MV_Vector<int>& MaxIter, double& C_lb, double& C_ub, MV_Vector<double>& Subdomain)
{
  	char text[100]; 
	double x; 
	int N, i; 

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> text;		// measured data file name with extension *.m
 	 	inpfile >> text >> text;		// file name with extension *.m  for true solution of the forward problem
 		inpfile >> text >> text;		// file name with extension *.m for incident wave 
 		inpfile >> text >> x >> x; 
		inpfile >> text >> N; // number of pseudo-frequencies

		if (N != Lambda_CWF.size())
			cout << endl << "Error in WavesOutputs::load_inversion_parameters: The input parameters are not consistent" << endl;
		else	
		{
			inpfile >> text;
			for (i=0; i< N; i++)
				inpfile >> Lambda_CWF(i);//load the Carlemann weighting factors

			inpfile >> text;
			for (i=0; i< N; i++)
				inpfile >> RegPar(i);//load the regularization parameters

			inpfile >> text;
			for (i=0; i< N; i++)
				inpfile >> MaxIter(i);//load the maximum number of iterations
	
			inpfile >> text >> C_lb >> C_ub; 
			
			inpfile >> text >> text; 
			inpfile >> text >> Subdomain(0) >> Subdomain(1) >> Subdomain(2) >> Subdomain(3) >> Subdomain(4) >> Subdomain(5); 

		} 
		inpfile.close();
	}
	else 
  		cout << endl << "Error in WavesOutputs::load_boundary_grid: file cannot be opened" << endl << endl;

}




// new input parameter file type: added by Thanh.

// -- load the input parameters from the .dat file:
void WavesOutputs::Configure(char *fname, Grid& grid, Grid& outer_gg, 
		 WavesSDGeometry& sdg, WavesSDGeometry& sdg_new, int& nsd, 
		 bool& USE_EXCHANGE, bool& USE_FEM, bool& USE_FDM, bool &USE_DIR_FEM, 
		 bool &USE_DIR_FDM, bool& USE_ABSORB, bool &USE_RHS,
		 bool &PRINT_FILES, int &NoTimeSteps, double &maxTime, double &rhs_code, 
		 double& noiselevel, double& freq, int& NoObjects, MV_Vector<double>& Zposition)

{
  	const int maxlen = 100;
  	double dxFDM, dyFDM, dzFDM, x_minFDM, y_minFDM, z_minFDM;
  	double dxFEM, dyFEM, dzFEM, x_minFEM, y_minFEM, z_minFEM;
 	int NxFDM, NyFDM, NzFDM, NxFEM, NyFEM, NzFEM;
	
	char text[maxlen], FDMGridFile[maxlen],	FEMGridFile[maxlen]; 
       	string	BoundaryGridFile, TypeOfBoundaryCondition, File1, File2, File3, File4;


  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> nsd; //number of space dimensions
 	 	inpfile >> text >> FDMGridFile;		// FDM grid file with extension *.inp
  		inpfile >> text >> FEMGridFile;		// FEM grid file with extension *.inp
	  	inpfile >> text >> BoundaryGridFile;		// file with extension *.inp for boundary nodes in the FEM domain
  		inpfile >> text >> NxFDM >> NyFDM >> NzFDM >> dxFDM >> dyFDM >> dzFDM
	  		>> x_minFDM >> y_minFDM >> z_minFDM; //FDM domain
  		inpfile >> text >> NxFEM >> NyFEM >> NzFEM >> dxFEM >> dyFEM >> dzFEM 
		  	>> x_minFEM >> y_minFEM >> z_minFEM; //FEM domain

  		inpfile >> text >> USE_EXCHANGE; // use the hybrid FDM-FEM method  
  		inpfile >> text >> USE_FEM; 
  		inpfile >> text >> USE_FDM; 
  		inpfile >> text >> USE_RHS; // use the right hand side (source term)
  		inpfile >> text >> USE_DIR_FEM; 
  		inpfile >> text >> USE_DIR_FDM; 
  		inpfile >> text >> USE_ABSORB; // use the absorbing b.c.
  		inpfile >> text >> PRINT_FILES;

  		inpfile >> text >> NoTimeSteps;
  		inpfile >> text >> maxTime; 
  		inpfile >> text >> rhs_code; //code for different right hand side functions

 		inpfile >> text >> TypeOfBoundaryCondition; 

  		inpfile >> text >> noiselevel; // noise level
 		inpfile >> text >> freq; // frequency of the incident plane wave 
  		inpfile >> text >> Zposition(0) >> Zposition(1) >> Zposition(2) >> Zposition(3); // values of z at which the data is saved
		inpfile >> text >> File1 >> File2 >> File3 >> File4; 		
		inpfile >> text >> NoObjects; 
		inpfile.close();

		grid.scan(FEMGridFile);
		outer_gg.scan(FDMGridFile);
		sdg.initialize(NxFDM, NyFDM, NzFDM, nsd, dxFDM, dyFDM, dzFDM, x_minFDM, y_minFDM, z_minFDM);
		sdg_new.initialize(NxFEM, NyFEM, NzFEM, nsd, dxFEM, dyFEM, dzFEM, x_minFEM, y_minFEM, z_minFEM); 

	}
	else 
	  	cout << endl << endl << "Error in WavesOutputs::Configure: file cannot be opened" << endl << endl;

}

// -- load the input parameters from the .dat file:
void WavesOutputs::Configure(char *fname, Grid& grid, Grid& outer_gg, WavesSDGeometry& sdg, WavesSDGeometry& sdg_new, int& nsd, 
		 bool& USE_EXCHANGE, bool& USE_FEM, bool& USE_FDM, bool &USE_DIR_FEM, bool &USE_DIR_FDM, bool& USE_ABSORB, 
		 int &NoTimeSteps, double &maxTime, double& freq)

{
  	const int maxlen = 100;
  	double dxFDM, dyFDM, dzFDM, x_minFDM, y_minFDM, z_minFDM;
  	double dxFEM, dyFEM, dzFEM, x_minFEM, y_minFEM, z_minFEM;
 	int NxFDM, NyFDM, NzFDM, NxFEM, NyFEM, NzFEM;
	
	char text[maxlen], FDMGridFile[maxlen],	FEMGridFile[maxlen]; 
       	string	BoundaryGridFile, TypeOfBoundaryCondition, File1, File2, File3, File4;

	int int_v; 	
	double doub_v;	

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> nsd; //number of space dimensions
 	 	inpfile >> text >> FDMGridFile;		// FDM grid file with extension *.inp
  		inpfile >> text >> FEMGridFile;		// FEM grid file with extension *.inp
	  	inpfile >> text >> BoundaryGridFile;		// file with extension *.inp for boundary nodes in the FEM domain
  		inpfile >> text >> NxFDM >> NyFDM >> NzFDM >> dxFDM >> dyFDM >> dzFDM
	  		>> x_minFDM >> y_minFDM >> z_minFDM; //FDM domain
  		inpfile >> text >> NxFEM >> NyFEM >> NzFEM >> dxFEM >> dyFEM >> dzFEM 
		  	>> x_minFEM >> y_minFEM >> z_minFEM; //FEM domain

  		inpfile >> text >> USE_EXCHANGE; // use the hybrid FDM-FEM method  
  		inpfile >> text >> USE_FEM; 
  		inpfile >> text >> USE_FDM; 
  		inpfile >> text >> int_v; // use the right hand side (source term)
  		inpfile >> text >> USE_DIR_FEM; 
  		inpfile >> text >> USE_DIR_FDM; 
  		inpfile >> text >> USE_ABSORB; // use the absorbing b.c.
  		inpfile >> text >> int_v;

  		inpfile >> text >> NoTimeSteps;
  		inpfile >> text >> maxTime; 
  		inpfile >> text >> doub_v; //code for different right hand side functions

 		inpfile >> text >> TypeOfBoundaryCondition; 

  		inpfile >> text >> doub_v; // noise level
 		inpfile >> text >> freq; // frequency of the incident plane wave 
 		inpfile.close();

		grid.scan(FEMGridFile);
		outer_gg.scan(FDMGridFile);
		sdg.initialize(NxFDM, NyFDM, NzFDM, nsd, dxFDM, dyFDM, dzFDM, x_minFDM, y_minFDM, z_minFDM);
		sdg_new.initialize(NxFEM, NyFEM, NzFEM, nsd, dxFEM, dyFEM, dzFEM, x_minFEM, y_minFEM, z_minFEM); 

	}
	else 
	  	cout << endl << endl << "Error in WavesOutputs::Configure: file cannot be opened" << endl << endl;

}



// -- load the input parameters from the .dat file:
void WavesOutputs::Configure(char *fname, double& x_minFEM, double& y_minFEM, double& z_minFEM, 
			    double& x_maxFEM, double& y_maxFEM, double& z_maxFEM, double& maxTime, int& NoTimeSteps)

{
  	const int maxlen = 100;
  	double dxFDM, dyFDM, dzFDM, x_minFDM, y_minFDM, z_minFDM;
  	double dxFEM, dyFEM, dzFEM;
 	int NxFDM, NyFDM, NzFDM, NxFEM, NyFEM, NzFEM;
	
	char text[maxlen]; 
  

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> NxFDM; //number of space dimensions
 	 	inpfile >> text >> text;		// FDM grid file with extension *.inp
  		inpfile >> text >> text;		// FEM grid file with extension *.inp
	  	inpfile >> text >> text;		// file with extension *.inp for boundary nodes in the FEM domain
  		inpfile >> text >> NxFDM >> NyFDM >> NzFDM >> dxFDM >> dyFDM >> dzFDM
	  		>> x_minFDM >> y_minFDM >> z_minFDM; //FDM domain
  		inpfile >> text >> NxFEM >> NyFEM >> NzFEM >> dxFEM >> dyFEM >> dzFEM 
		  	>> x_minFEM >> y_minFEM >> z_minFEM; //FEM domain

		for (int i = 0; i < 8; i++)
	  		inpfile >> text >> NxFDM;

  		inpfile >> text >> NoTimeSteps;
  		inpfile >> text >> maxTime; 
 		inpfile.close();
		
		x_maxFEM = x_minFEM + (NxFEM-1)*dxFEM;
		y_maxFEM = y_minFEM + (NyFEM-1)*dyFEM;
		z_maxFEM = z_minFEM + (NzFEM-1)*dzFEM;

	}
	else 
	  	cout << endl << endl << "Error in WavesOutputs::Configure: file cannot be opened" << endl << endl;

}



// load the boundary grid file:
string WavesOutputs::boundary_grid_file_name(char *fname)
{
  	const int maxlen = 100;
	int nsd;  
  	char text[maxlen], FDMGridFile[maxlen],	FEMGridFile[maxlen];

  	string BoundaryGridFile;
  
  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> nsd; //number of space dimensions
  		inpfile >> text >> FDMGridFile;		// FDM grid file with extension *.inp
  		inpfile >> text >> FEMGridFile;		// FEM grid file with extension *.inp
  		inpfile >> text >> BoundaryGridFile;		// file with extension *.inp for boundary nodes in the FEM domain
 
//  		cout << endl << BoundaryGridFile << endl;
 
	}
	else
	{
  		cout << endl << endl << "Error in WavesOutputs::boundary_grid_file_name: file cannot be opened" << endl << endl;
  		BoundaryGridFile = " "; 
	}	
  	return BoundaryGridFile; 
}

// load the boundary grid file:
string WavesOutputs::FEM_grid_file_name(char *fname)
{
  	const int maxlen = 100;
	int nsd;  
  	char text[maxlen], FDMGridFile[maxlen];

  	string FEMGridFile;
  
  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> nsd; //number of space dimensions
  		inpfile >> text >> FDMGridFile;		// FDM grid file with extension *.inp
  		inpfile >> text >> FEMGridFile;		// FEM grid file with extension *.inp
 
//  		cout << endl << FEMGridFile << endl;
 
	}
	else
	{
  		cout << endl << endl << "Error in WavesOutputs::FEM_grid_file_name: file cannot be opened" << endl << endl;
  		FEMGridFile = " "; 
	}	
  	return FEMGridFile; 
}


//extract the data file names for saving the simulated solution:
string WavesOutputs::Configure(char *fname, int idx)
{
  	const int maxlen = 100;
  	double dxFDM, dyFDM, dzFDM, x_minFDM, y_minFDM, z_minFDM;
  	double dxFEM, dyFEM, dzFEM, x_minFEM, y_minFEM, z_minFEM;
 	int NxFDM, NyFDM, NzFDM, NxFEM, NyFEM, NzFEM, nsd, i;
  	
	char text[maxlen], FDMGridFile[maxlen],	FEMGridFile[maxlen]; 
       	string	BoundaryGridFile, TypeOfBoundaryCondition, File1, File2, File3, File4;


  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> nsd; //number of space dimensions
 	 	inpfile >> text >> FDMGridFile;		// FDM grid file with extension *.inp
  		inpfile >> text >> FEMGridFile;		// FEM grid file with extension *.inp
	  	inpfile >> text >> BoundaryGridFile;		// file with extension *.inp for boundary nodes in the FEM domain
  		inpfile >> text >> NxFDM >> NyFDM >> NzFDM >> dxFDM >> dyFDM >> dzFDM
	  		>> x_minFDM >> y_minFDM >> z_minFDM; //FDM domain
  		inpfile >> text >> NxFEM >> NyFEM >> NzFEM >> dxFEM >> dyFEM >> dzFEM 
		  	>> x_minFEM >> y_minFEM >> z_minFEM; //FEM domain

		for (i = 0; i < 8; i++)
	  		inpfile >> text >> nsd;
		for (i = 0; i < 3; i++)
	  		inpfile >> text >> dxFDM;

 		inpfile >> text >> TypeOfBoundaryCondition; 

  		inpfile >> text >> dxFDM; // noise level
 		inpfile >> text >> dxFDM; // frequency of the incident plane wave 
  		inpfile >> text >> dxFDM >> dxFDM >> dxFDM >> dxFDM; // values of z at which the data is saved
		inpfile >> text >> File1 >> File2 >> File3 >> File4; 	
		inpfile.close();
  		if (idx == 1)
  			return File1;
  		else if (idx == 2)
  			return File2;
  		else if (idx == 3)
  			return File3;
  		else 
  			return File4;
  	}
	else 
	{
  		cout << endl << endl << "Error in WavesOutputs::Configure: file cannot be opened" << endl << endl;
  		return " ";
	}

}

// load the boundary condition type:
string WavesOutputs::type_of_boundary_condition(char *fname)
{
  	const int maxlen = 100;
  	double dxFDM, dyFDM, dzFDM, x_minFDM, y_minFDM, z_minFDM;
  	double dxFEM, dyFEM, dzFEM, x_minFEM, y_minFEM, z_minFEM;
 	int NxFDM, NyFDM, NzFDM, NxFEM, NyFEM, NzFEM, nsd, i;
  	
	char text[maxlen], FDMGridFile[maxlen],	FEMGridFile[maxlen]; 
       	string	BoundaryGridFile, TypeOfBoundaryCondition;


  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> nsd; //number of space dimensions
 	 	inpfile >> text >> FDMGridFile;		// FDM grid file with extension *.inp
  		inpfile >> text >> FEMGridFile;		// FEM grid file with extension *.inp
	  	inpfile >> text >> BoundaryGridFile;		// file with extension *.inp for boundary nodes in the FEM domain
  		inpfile >> text >> NxFDM >> NyFDM >> NzFDM >> dxFDM >> dyFDM >> dzFDM
	  		>> x_minFDM >> y_minFDM >> z_minFDM; //FDM domain
  		inpfile >> text >> NxFEM >> NyFEM >> NzFEM >> dxFEM >> dyFEM >> dzFEM 
		  	>> x_minFEM >> y_minFEM >> z_minFEM; //FEM domain

		for (i = 0; i < 8; i++)
	  		inpfile >> text >> nsd; // use the hybrid FDM-FEM method  
		for (i = 0; i < 3; i++)
	  		inpfile >> text >> dxFDM;

 		inpfile >> text >> TypeOfBoundaryCondition; 

 		inpfile.close();
		return TypeOfBoundaryCondition; 
  	}
	else 
	{
  		cout << endl << endl << "Error in WavesOutputs::Configure: file cannot be opened" << endl << endl;
  		return " ";
	}

}

// -- load the input parameters from the .dat file:
void WavesOutputs::load_FEM_parameters(char *fname, int& nsd, int& NxFEM, int& NyFEM, int& NzFEM, 
			double& dxFEM, double& dyFEM, double& dzFEM, 
			double& x_minFEM, double& y_minFEM, double& z_minFEM, int& NoObjects)

{
  	const int maxlen = 100;
  	double dxFDM, dyFDM, dzFDM, x_minFDM, y_minFDM, z_minFDM;
  	int NxFDM, NyFDM, NzFDM, i;
	
	char text[maxlen], FDMGridFile[maxlen]; 

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> nsd; //number of space dimensions
		for (i = 0; i < 3; i++)
	 	 	inpfile >> text >> FDMGridFile;		// FDM grid file with extension *.inp
  		
		inpfile >> text >> NxFDM >> NyFDM >> NzFDM >> dxFDM >> dyFDM >> dzFDM
	  		>> x_minFDM >> y_minFDM >> z_minFDM; //FDM domain
  		inpfile >> text >> NxFEM >> NyFEM >> NzFEM >> dxFEM >> dyFEM >> dzFEM 
		  	>> x_minFEM >> y_minFEM >> z_minFEM; //FEM domain

		for (i=0;i <8; i++)
			inpfile >> text >> NxFDM; 
		for (i=0;i <3; i++)
			inpfile >> text >> dxFDM;  
  		
 		inpfile >> text >> text; 

  		inpfile >> text >> dxFDM; // noise level
 		inpfile >> text >> dxFDM; // frequency of the incident plane wave 
  		inpfile >> text >> dxFDM >> dxFDM >> dxFDM >> dxFDM; // values of z at which the data is saved
		inpfile >> text >> text >> text >> text >> text;	
		inpfile >> text >> NoObjects; 
		inpfile.close();

	}
	else 
	  	cout << endl << endl << "Error in WavesOutputs::load_FEM_parameters: file cannot be opened" << endl << endl;

}

// -- load the object parameters from the .dat file:
string WavesOutputs::load_object_type(char *fname, int objidx)
{
 	const int maxlen = 100;
  	double dxFDM, dyFDM, dzFDM, x_minFDM, y_minFDM, z_minFDM;
  	int NxFDM, NyFDM, NzFDM, i;
	
	char text[maxlen], FDMGridFile[maxlen]; 
	string ObjectType;

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> NxFDM; //number of space dimensions
		for (i = 0; i < 3; i++)
	 	 	inpfile >> text >> FDMGridFile;		// FDM grid file with extension *.inp
  		
		inpfile >> text >> NxFDM >> NyFDM >> NzFDM >> dxFDM >> dyFDM >> dzFDM
	  		>> x_minFDM >> y_minFDM >> z_minFDM; //FDM domain
 		inpfile >> text >> NxFDM >> NyFDM >> NzFDM >> dxFDM >> dyFDM >> dzFDM
	  		>> x_minFDM >> y_minFDM >> z_minFDM; //FDM domain
 
		for (i=0;i <8; i++)
			inpfile >> text >> NxFDM; 
		for (i=0;i <3; i++)
			inpfile >> text >> dxFDM;  
  		
 		inpfile >> text >> text; 

  		inpfile >> text >> dxFDM; // noise level
 		inpfile >> text >> dxFDM; // frequency of the incident plane wave 
  		inpfile >> text >> dxFDM >> dxFDM >> dxFDM >> dxFDM; // values of z at which the data is saved
		inpfile >> text >> text >> text >> text >> text;	
		inpfile >> text >> NxFDM; 

		// avoid loading previous objects:
		if (objidx > 0)
		{	for (i = 0; i < objidx; i++)
			{
				inpfile >> text >> text; // object type
				inpfile >> text >> NxFDM;
				inpfile >> text >> dxFDM;
				inpfile >> text >> dxFDM >> dxFDM >> dxFDM >> dxFDM >> dxFDM >> dxFDM;
			}
		}
		inpfile >> text >> ObjectType; 	
	
		inpfile.close();
	}
	else 
	{
	  	cout << endl << endl << "Error in WavesOutputs::load_object_type: file cannot be opened" << endl << endl;
		ObjectType = " ";
	}

	return ObjectType; 

}

// -- load the object parameters from the .dat file:
void WavesOutputs::load_object_parameters(char *fname, int objidx, int& TypeOfMat, double& velocity, MV_Vector<double>& Coord)
{
	const int maxlen = 100;
  	double dxFDM, dyFDM, dzFDM, x_minFDM, y_minFDM, z_minFDM;
  	int NxFDM, NyFDM, NzFDM, i;
	
	char text[maxlen], FDMGridFile[maxlen]; 
	string ObjectType;

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> NxFDM; //number of space dimensions
		for (i = 0; i < 3; i++)
	 	 	inpfile >> text >> FDMGridFile;		// FDM grid file with extension *.inp
  		
		inpfile >> text >> NxFDM >> NyFDM >> NzFDM >> dxFDM >> dyFDM >> dzFDM
	  		>> x_minFDM >> y_minFDM >> z_minFDM; //FDM domain
 		inpfile >> text >> NxFDM >> NyFDM >> NzFDM >> dxFDM >> dyFDM >> dzFDM
	  		>> x_minFDM >> y_minFDM >> z_minFDM; //FDM domain
 
		for (i=0;i <8; i++)
			inpfile >> text >> NxFDM; 
		for (i=0;i <3; i++)
			inpfile >> text >> dxFDM;  
  		
 		inpfile >> text >> text; 

  		inpfile >> text >> dxFDM; // noise level
 		inpfile >> text >> dxFDM; // frequency of the incident plane wave 
  		inpfile >> text >> dxFDM >> dxFDM >> dxFDM >> dxFDM; // values of z at which the data is saved
		inpfile >> text >> text >> text >> text >> text;	
		inpfile >> text >> NxFDM; 

		// avoid loading previous objects:
		if (objidx > 0)
		{	for (i = 0; i < objidx; i++)
			{
				inpfile >> text >> text; // object type
				inpfile >> text >> NxFDM;
				inpfile >> text >> dxFDM;
				inpfile >> text >> dxFDM >> dxFDM >> dxFDM >> dxFDM >> dxFDM >> dxFDM;
			}
		}
		inpfile >> text >> ObjectType; 		
		inpfile >> text >> TypeOfMat;  
		inpfile >> text >> velocity; 
		inpfile >> text >> Coord(0) >> Coord(1) >> Coord(2) >> Coord(3) >> Coord(4) >> Coord(5); 
		
		inpfile.close();
	}
	else 
	  	cout << endl << endl << "Error in WavesOutputs::load_object_parameters: file cannot be opened" << endl << endl;

}

// -- load the object parameters from the .dat file:
void WavesOutputs::load_object_parameters(char *fname, int objidx, int& TypeOfMat, double& velocity)
{
	const int maxlen = 100;
  	double dxFDM, dyFDM, dzFDM, x_minFDM, y_minFDM, z_minFDM;
  	int NxFDM, NyFDM, NzFDM, i;
	
	char text[maxlen], FDMGridFile[maxlen]; 
	string ObjectType;

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> NxFDM; //number of space dimensions
		for (i = 0; i < 3; i++)
	 	 	inpfile >> text >> FDMGridFile;		// FDM grid file with extension *.inp
  		
		inpfile >> text >> NxFDM >> NyFDM >> NzFDM >> dxFDM >> dyFDM >> dzFDM
	  		>> x_minFDM >> y_minFDM >> z_minFDM; //FDM domain
 		inpfile >> text >> NxFDM >> NyFDM >> NzFDM >> dxFDM >> dyFDM >> dzFDM
	  		>> x_minFDM >> y_minFDM >> z_minFDM; //FDM domain
 
		for (i=0;i <8; i++)
			inpfile >> text >> NxFDM; 
		for (i=0;i <3; i++)
			inpfile >> text >> dxFDM;  
  		
 		inpfile >> text >> text; 

  		inpfile >> text >> dxFDM; // noise level
 		inpfile >> text >> dxFDM; // frequency of the incident plane wave 
  		inpfile >> text >> dxFDM >> dxFDM >> dxFDM >> dxFDM; // values of z at which the data is saved
		inpfile >> text >> text >> text >> text >> text;	
		inpfile >> text >> NxFDM; 

		// avoid loading previous objects:
		if (objidx > 0)
		{	for (i = 0; i < objidx; i++)
			{
				inpfile >> text >> text; // object type
				inpfile >> text >> NxFDM;
				inpfile >> text >> dxFDM;
				inpfile >> text >> dxFDM >> dxFDM >> dxFDM >> dxFDM >> dxFDM >> dxFDM;
			}
		}
		inpfile >> text >> ObjectType; 		
		inpfile >> text >> TypeOfMat;  
		inpfile >> text >> velocity; 
			
		inpfile.close();
	}
	else 
	  	cout << endl << endl << "Error in WavesOutputs::load_object_parameters: file cannot be opened" << endl << endl;

}


void WavesOutputs::load_boundary_data_inversion(char* fname, MV_ColMat<double>& psi, MV_Vector<int>& BndNodes, MV_Vector<double>& bdvalue_tail)
{
	int NoBndNodes = BndNodes.size(); 
	int Ns = psi.size(0);

 	ifstream inpfile;

  	inpfile.open(fname);
	if (inpfile.is_open())
	{
		int i,j; 
		for (j=0;j < NoBndNodes; j++)
		{	inpfile >> BndNodes(j); }
		for (i = 0; i < Ns; i++)
		{	for (j=0;j < NoBndNodes; j++)
			{	inpfile >> psi(i,j); 
			}
		}
		for (j=0;j < NoBndNodes; j++)
		{	inpfile >> bdvalue_tail(BndNodes(j)-1); }
		inpfile.close();

	}
	else
		cout << endl << endl << "Error in WavesOutputs::load_boundary_data: file cannot be opened" << endl << endl;
}

// load the information of a grid file *.inp 
void WavesOutputs::load_grid_file(char *fname, MV_ColMat<double>& Nodes, MV_ColMat<int>& Elements)
{
	ifstream inp;
	inp.open(fname);
	int i, j, Nno, Nel;
	char test[10];
	//printf("Open the file: %s\n", fname);
	
	// load the number of nodes, number of elements:
	inp >> Nno >> Nel >> i >> i >> i; // the variable i is not used. 

	//test resizing the Nodes and Elements:
	Nodes.newsize(Nno,4); 
	Elements.newsize(Nel,4);

	//load the nodes' coordinates:
	for (i = 0; i < Nno; i++)
	{
		inp >> j >> Nodes(i, 0) >> Nodes(i,1) >> Nodes(i,2); // coordinate of the nodes
	}

	//load the element's nodes:
	for (i = 0; i < Nel; i++)
	{
		inp >> j >> j >> test >> Elements(i,0) >> Elements(i,1) >> Elements(i,2) >> Elements(i,3); 
	}
	
	// load the material properties of the nodes:
	inp >> j >> j; //ignore these two data lines
	inp >> test >> j >> test >> test; 

	for (i = 0; i < Nno; i++)
	{
		inp >> j >> Nodes(i,3);
	}
	inp.close();
}


void WavesOutputs::write2tecplot(char* file, MV_ColMat<double> Nodes, MV_ColMat<int> Elements) // write the data to the tecplot FE data format
{

	int i, j, Nno, Nel;
	Nno = Nodes.size(0); //number of nodes
	Nel = Elements.size(0); // number of elements

	ofstream outp;
	outp.open(file);

	cout << "opened file  " << file << endl;

	// write 3 lines of the header: 
	outp << "TITLE = \" The density\" " << endl;
	outp << "VARIABLES = \"X\", \"Y\", \"Z\", \"Density\" " << endl;
	outp << "ZONE N = " << Nno << ", E = " << Nel << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON " << endl;

	// write the nodes: 
	outp << "# Nodes: " << endl; 
	for (i = 0; i < Nno; i++)
	  {
		outp << Nodes(i,0) << " " << Nodes(i,1) << " " << Nodes(i,2) << " " << Nodes(i,3) << " " << endl;
	  }
	// write the elements:	
	outp << "# Elements: " << endl;
	for (i = 0; i < Nel; i++)
	  {
		outp << Elements(i,0) << " " << Elements(i,1) << " " << Elements(i,2) << " " << Elements(i,3) << " " << endl;
	  }

	cout << "closed file " << file << endl;
	outp.close();
}

