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

#include "include/wavesGridB.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

using namespace std;


WavesGridB::WavesGridB()

{
	// initialize the WavesGridB data
	nsd = nno = nel = maxnne = nbind = 0;
	onemat = false;
	//  surface_triangulation = false;
	elm_type_for_grid = NO_ELEMENT;
//	printf("WavesGridB::WavesGridB creating object %#lx\n", (long) this);
}


WavesGridB::~WavesGridB()

{
//	printf("WavesGridB::~WavesGridB: deleting object %#lx, %d nodes and %d elements\n", (long) this, nno, nel);
}


bool WavesGridB::ok() const

{
// initialize the WavesGridB data
	if (nno > 0)
		return true;
	else
		return false;
}


void WavesGridB::operator =(const WavesGridB& grid)

{
	if (!grid.ok())
		printf(" fatal error WavesGridB:: operator = Empty grid is given\n");

	this->redim(grid.nsd, grid.nno, grid.nel, grid.maxnne, grid.nbind, grid.elm_type_for_grid);
	this->elm_type_for_grid = grid.elm_type_for_grid;
	this->element_type = grid.element_type;
	this->elid = grid.elid;
	this->nodpek = grid.nodpek;
	this->coord = grid.coord;
//	int aa = grid.nbind;
	if (grid.bind.size(1))
		this->bind = grid.bind;
//	printf("WavesGridB:: operator = \n");
}


bool WavesGridB::redim(const int nsd_, const int nno_, int const nel_, const int maxnne_, const int nbind_, const ElementType e_type)

{
	//printf("WavesGridB:: redim, nsd=%i,nno=%i,nel=%i,maxnne=%i,nbind=%i\n", nsd_, nno_, nel_, maxnne_, nbind_);

	nsd = nsd_;
	nno = nno_;
	nel = nel_;
	maxnne = maxnne_;
	nbind = nbind_;
	elid.newsize(nel);
	elid = 0; // initialize to material 0
	nodpek.newsize(nel, maxnne);
	coord.newsize(nno, nsd);
	coord = 0.0;
	bind.newsize(nno, nbind);
	bind = 0;
	if (e_type == NO_ELEMENT)
	{
		elm_type_for_grid = NO_ELEMENT;
		element_type.newsize(nel);
	}
	else
	{
		elm_type_for_grid = e_type;
		element_type.newsize(0);
	}
	return true;
}

int str_is_in(const char *str, char *tag)
{
	while (*str != '\0')
	{
		if (*str++ == *tag)
		{
			tag++;
			while (*tag != '\0')
				if (*str++ != *tag++)
					return 0;
			return 1;
		}
	}
	return 0;
}

void WavesGridB::scan(const char *file)

{
	const char *fileptr = file;
	if (str_is_in(fileptr, (char *) ".inp"))
	{
		scan_avs(file);
	}
	else
	{
		printf("WavesGridB:: scan, failed to find scan routine for %s\n", file);
		exit(-1);
	}
	return;
}

ElementType fromAVStext2eltype(char *c)
{
	ElementType etype = NO_ELEMENT;
	if (str_is_in(c, (char *) "tet"))
	{
		etype = ELMTET1;
	}
	else if (str_is_in(c, (char *) "hex"))
	{
		etype = ELMHEX1;
	}
	else if (str_is_in(c, (char *) "pyr"))
	{
		etype = ELMPYR1;
	}
	else if (str_is_in(c, (char *) "prism"))
	{
		etype = ELMPRISM1;
	}
	else if (str_is_in(c, (char *) "tri"))
	{
		etype = ELMTRI1;
	}
	else if (str_is_in(c, (char *) "quad"))
	{
		etype = ELMQUAD1;
	}
	return etype;
}


void read_max_elm_nodes(const char *file, MV_Vector<ElementType>& etypes, int& maximum_node_in_elm, int& nospacedim)

{

	FILE *fp;
	int i, j;
	int el;
	int d, id;
	double fdummy;
	char c[5];
	fp = fopen(file, "r");
	//printf("Reading from %s\n", file);

	ifstream isnode(file, ios::in);
	if (!isnode)
	{
		printf("error: read_avs The node file was not found\n");
		return;
	}
	int nnode, elnum;
	fscanf(fp, "%i %i %i %i %i", &nnode, &elnum, &i, &i, &i);
	isnode >> nnode >> elnum >> i >> i >> i;
	//printf(" There are %d nodes and %d elements\n", nnode, elnum);

	maximum_node_in_elm = 0; // just an initialization
	int locnsd = 0;
	nospacedim = 0; //initialization
	etypes.newsize(elnum);
	for (i = 0; i < nnode; i++)
	{
		fscanf(fp, "%i %lf %lf %lf\n", &j, &fdummy, &fdummy, &fdummy);
	}
	int nbf = 0;
	for (el = 0; el < elnum; el++)
	{
		fscanf(fp, "%i %i %s\n", &d, &id, c);
		ElementType etype = fromAVStext2eltype(c);
		etypes(el) = etype;
		switch (etype)
		{
			case ELMTET1:
				nbf = 4;
				locnsd = 3;
				break;
			case ELMPYR1:
				nbf = 5;
				locnsd = 3;
				break;
			case ELMHEX1:
				nbf = 8;
				locnsd = 3;
				break;
			case ELMPRISM1:
				nbf = 6;
				locnsd = 3;
				break;
			case ELMTRI1:
				nbf = 3;
				locnsd = 2;
				break;
			case ELMQUAD1:
				nbf = 4;
				locnsd = 2;
				break;
			case NO_ELEMENT:
				printf("error read_max_elm_nodes wrong type of element read\n");
				break;
		}
		for (i = 0; i < nbf; i++)
			fscanf(fp, "%lf", &fdummy);
		if (maximum_node_in_elm < nbf)
			maximum_node_in_elm = nbf;
		nospacedim = locnsd;
	}
	fclose(fp);
}


void WavesGridB::scan_avs(const char *file)

{

	FILE *fp;
	int i, j;
	int el;
	int d, id;
	double dd;
	char c[5];
	fp = fopen(file, "r");
	//printf("Reading from %s\n", file);

	ifstream isnode(file, ios::in);
	if (!isnode)
	{
		printf("error: read_avs The node file was not found\n");
		return;
	}
	int nnode, elnum;
	fscanf(fp, "%i %i %i %i %i", &nnode, &elnum, &i, &i, &i);
	isnode >> nnode >> elnum >> i >> i >> i;
	//printf(" There are %d nodes and %d elements\n", nnode, elnum);

	nno = nnode;
	nel = elnum;
	MV_Vector<ElementType> etypes;
	read_max_elm_nodes(file, etypes, maxnne, nsd);
	nbind = 0;

	redim(nsd, nno, nel, maxnne, nbind, NO_ELEMENT);
	element_type = etypes;

	if (nsd == 3)
		for (i = 0; i < nnode; i++)
		{
			fscanf(fp, "%i %lf %lf %lf\n", &j, &coord(i, 0), &coord(i, 1), &coord(i, 2));
		}
	else if (nsd == 2)
		for (i = 0; i < nnode; i++)
		{
			fscanf(fp, "%i %lf %lf %lf\n", &j, &coord(i, 0), &coord(i, 1), &dd);
		}

	int nbf = 0;
	for (el = 0; el < elnum; el++)
	{
		fscanf(fp, "%i %i %s\n", &d, &id, c);
		ElementType etype = element_type(el);
		switch (etype)
		{
			case ELMTET1:
				nbf = 4;
				fscanf(fp, "%i %i %i %i\n", &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3), &nodpek(el, 0));
				break;
			case ELMPYR1:
				nbf = 5;
				fscanf(fp, "%i %i %i %i %i\n", &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3), &nodpek(el, 4), &nodpek(el, 0));
				break;
			case ELMHEX1:
				nbf = 8;
				fscanf(fp, "%i %i %i %i %i %i %i %i\n", &nodpek(el, 4), &nodpek(el, 5), &nodpek(el, 7), &nodpek(el, 6), &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 3), &nodpek(el, 2));
				break;
			case ELMPRISM1:
				nbf = 6;
				fscanf(fp, "%i %i %i %i %i %i\n", &nodpek(el, 3), &nodpek(el, 4), &nodpek(el, 5), &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2));
				break;
			case ELMTRI1:
				nbf = 3;
				fscanf(fp, "%i %i %i\n", &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2));
				break;
			case ELMQUAD1:
				nbf = 4;
				fscanf(fp, "%i %i %i %i\n", &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3));
				break;
			case NO_ELEMENT:
				printf("error scan_avs wrong type of element read\n");
				break;
		}
		for (i = 0; i < nbf; i++)
			nodpek(el, i)--;elid
		(el) = id - 1;
	}

	fclose(fp);

	checkIfOnlyOneElmTypeInGrid();

}


void WavesGridB::scan_avs_amira(const char *file)

{

	FILE *fp;
	int i, j;
	int el;
	int d, id;
	double dd;
	char c[5];
	fp = fopen(file, "r");
	//printf("Reading from %s\n", file);

	ifstream isnode(file, ios::in);
	if (!isnode)
	{
		printf("error: read_avs The node file was not found\n");
		return;
	}
	int nnode, elnum;
	fscanf(fp, "%i %i %i %i %i", &nnode, &elnum, &i, &i, &i);
	isnode >> nnode >> elnum >> i >> i >> i;
	//printf(" There are %d nodes and %d elements\n", nnode, elnum);

	nno = nnode;
	nel = elnum;
	MV_Vector<ElementType> etypes;
	read_max_elm_nodes(file, etypes, maxnne, nsd);
	nbind = 0;

	redim(nsd, nno, nel, maxnne, nbind, NO_ELEMENT);
	element_type = etypes;

	if (nsd == 3)
		for (i = 0; i < nnode; i++)
		{
			fscanf(fp, "%i %lf %lf %lf\n", &j, &coord(i, 0), &coord(i, 1), &coord(i, 2));
		}
	else if (nsd == 2)
		for (i = 0; i < nnode; i++)
		{
			fscanf(fp, "%i %lf %lf %lf\n", &j, &coord(i, 0), &coord(i, 1), &dd);
		}

	int nbf = 0;
	for (el = 0; el < elnum; el++)
	{
		fscanf(fp, "%i %i %s\n", &d, &id, c);
		ElementType etype = element_type(el);
		switch (etype)
		{
			case ELMTET1:
				nbf = 4;
				fscanf(fp, "%i %i %i %i\n", &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3));
				break;
			case ELMPYR1:
				nbf = 5;
				fscanf(fp, "%i %i %i %i %i\n", &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3), &nodpek(el, 4), &nodpek(el, 0));
				break;
			case ELMHEX1:
				nbf = 8;
				fscanf(fp, "%i %i %i %i %i %i %i %i\n", &nodpek(el, 4), &nodpek(el, 5), &nodpek(el, 6), &nodpek(el, 7), &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3));
				break;
			case ELMPRISM1:
				nbf = 6;
				fscanf(fp, "%i %i %i %i %i %i\n", &nodpek(el, 3), &nodpek(el, 4), &nodpek(el, 5), &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2));
				break;
			case ELMTRI1:
				nbf = 3;
				fscanf(fp, "%i %i %i\n", &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2));
				break;
			case ELMQUAD1:
				nbf = 4;
				fscanf(fp, "%i %i %i %i\n", &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3));
				break;
			case NO_ELEMENT:
				printf("error scan_avs wrong type of element read\n");
				break;
		}
		for (i = 0; i < nbf; i++)
			elid(el) = id;
	}

	fclose(fp);

	checkIfOnlyOneElmTypeInGrid();

}


void WavesGridB::scan_avs_vel(const char *file, MV_Vector<double>& velocity)

{

	FILE *fp;
	int i, j;
	int el;
	int d, id;
	double dd;
	char c[5];
	fp = fopen(file, "r");
	//printf("Reading from %s\n", file);

	ifstream isnode(file, ios::in);
	if (!isnode)
	{
		printf("error: read_avs The node file was not found\n");
		return;
	}
	int nnode, elnum;
	fscanf(fp, "%i %i %i %i %i", &nnode, &elnum, &i, &i, &i);
	isnode >> nnode >> elnum >> i >> i >> i;
	//printf(" There are %d nodes and %d elements\n", nnode, elnum);

	nno = nnode;
	nel = elnum;
	MV_Vector<ElementType> etypes;
	read_max_elm_nodes(file, etypes, maxnne, nsd);
	nbind = 0;

	redim(nsd, nno, nel, maxnne, nbind, NO_ELEMENT);
	element_type = etypes;

	if (nsd == 3)
		for (i = 0; i < nnode; i++)
		{
			fscanf(fp, "%i %lf %lf %lf\n", &j, &coord(i, 0), &coord(i, 1), &coord(i, 2));
		}
	else if (nsd == 2)
		for (i = 0; i < nnode; i++)
		{
			fscanf(fp, "%i %lf %lf %lf\n", &j, &coord(i, 0), &coord(i, 1), &dd);
		}

	int nbf = 0;
	for (el = 0; el < elnum; el++)
	{
		fscanf(fp, "%i %i %s\n", &d, &id, c);
		ElementType etype = element_type(el);
		switch (etype)
		{
			case ELMTET1:
				nbf = 4;
				fscanf(fp, "%i %i %i %i\n", &nodpek(el, 1), &nodpek(el, 3), &nodpek(el, 2), &nodpek(el, 0));
				break;
			case ELMPYR1:
				nbf = 5;
				fscanf(fp, "%i %i %i %i %i\n", &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3), &nodpek(el, 4), &nodpek(el, 0));
				break;
			case ELMHEX1:
				nbf = 8;
				fscanf(fp, "%i %i %i %i %i %i %i %i\n", &nodpek(el, 4), &nodpek(el, 5), &nodpek(el, 6), &nodpek(el, 7), &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3));
				break;
			case ELMPRISM1:
				nbf = 6;
				fscanf(fp, "%i %i %i %i %i %i\n", &nodpek(el, 3), &nodpek(el, 4), &nodpek(el, 5), &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2));
				break;
			case ELMTRI1:
				nbf = 3;
				fscanf(fp, "%i %i %i\n", &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2));
				break;
			case ELMQUAD1:
				nbf = 4;
				fscanf(fp, "%i %i %i %i\n", &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3));
				break;
			case NO_ELEMENT:
				printf("error scan_avs wrong type of element read\n");
				break;
		}
		for (i = 0; i < nbf; i++)
			nodpek(el, i)--;elid
		(el) = id - 1;
	}
	char ee[5];
	int a, b;

	fscanf(fp, " %i %i\n", &a, &b);

	fscanf(fp, "%s\n", ee);

	cout << ee << endl;

	velocity.newsize(nnode);
	int count;
	for (i = 0; i < nnode; i++)
	{
		fscanf(fp, "%i  %lf\n", &count, &velocity(i));
		cout << "velocity(" << count << ") = " << velocity(i) << endl;
	}

	fclose(fp);

	checkIfOnlyOneElmTypeInGrid();

}


void WavesGridB::scan_gid(const char *file, int elnum1)

{

	FILE *fp;
	int i, j;
	int el;
	int d, id;
	double dd;
	fp = fopen(file, "r");
	//printf("Reading from %s\n", file);

	ifstream isnode(file, ios::in);
	if (!isnode)
	{
		printf("error: read_gid The node file was not found\n");
		return;
	}
	int nnode, elnum;
	fscanf(fp, "%i %i %i %i %i", &nnode, &elnum, &i, &i, &i);
	isnode >> nnode >> elnum >> i >> i >> i;
	//printf(" There are %d nodes and %d elements\n", nnode, elnum);
	maxnne = 4;
	nsd = 3;
	nno = nnode;
	nel = elnum;
	MV_Vector<ElementType> etypes(elnum);
	nbind = 0;

	redim(nsd, nno, nel, maxnne, nbind, NO_ELEMENT);

	for (el = 0; el < elnum; el++)
	{
		etypes(el) = ELMTET1;
		cout << "etypes" << etypes(el) << endl;
	}

	if (nsd == 3)
		for (i = 0; i < nnode; i++)
		{
			fscanf(fp, "%i %lf %lf %lf\n", &j, &coord(i, 0), &coord(i, 1), &coord(i, 2));
		}
	else if (nsd == 2)
		for (i = 0; i < nnode; i++)
		{
			fscanf(fp, "%i %lf %lf %lf\n", &j, &coord(i, 0), &coord(i, 1), &dd);
		}

	cout << "before for elnum" << elnum << endl;

	int nbf = 0;

	for (el = 0; el < elnum1; el++)
	{

		// read number of element
		ElementType etype = etypes(el);
		switch (etype)
		{
			case ELMTET1:
				nbf = 4;
				fscanf(fp, "%i %i %i %i %i %i \n", &d, &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3), &nodpek(el, 0), &id);
				cout << nodpek(el, 0) << "  " << nodpek(el, 1) << "  " << nodpek(el, 2) << "  " << nodpek(el, 3) << endl;
				break;
			case ELMPYR1:
				nbf = 5;
				fscanf(fp, "%i %i %i %i %i\n", &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3), &nodpek(el, 4), &nodpek(el, 0));
				break;
			case ELMHEX1:
				nbf = 8;
				fscanf(fp, "%i %i %i %i %i %i %i %i\n", &nodpek(el, 4), &nodpek(el, 5), &nodpek(el, 6), &nodpek(el, 7), &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3));
				break;
			case ELMPRISM1:
				nbf = 6;
				fscanf(fp, "%i %i %i %i %i %i\n", &nodpek(el, 3), &nodpek(el, 4), &nodpek(el, 5), &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2));
				break;
			case ELMTRI1:
				nbf = 3;
				fscanf(fp, "%i %i %i\n", &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2));
				break;
			case ELMQUAD1:
				nbf = 4;
				fscanf(fp, "%i %i %i %i\n", &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3));
				break;
			case NO_ELEMENT:
				printf("error scan_avs wrong type of element read\n");
				break;
		}
		for (i = 0; i < nbf; i++)
			nodpek(el, i)--;elid
		(el) = id;

	}

	for (el = elnum1; el < elnum; el++)
	{

		// read number of element
		ElementType etype = etypes(el);

		switch (etype)
		{
			case ELMTET1:
				nbf = 4;
				fscanf(fp, "%i %i %i %i %i  \n", &d, &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3), &nodpek(el, 0));
				cout << nodpek(el, 0) << "  " << nodpek(el, 1) << "  " << nodpek(el, 2) << "  " << nodpek(el, 3) << endl;
				break;
			case ELMPYR1:
				nbf = 5;
				fscanf(fp, "%i %i %i %i %i\n", &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3), &nodpek(el, 4), &nodpek(el, 0));
				break;
			case ELMHEX1:
				nbf = 8;
				fscanf(fp, "%i %i %i %i %i %i %i %i\n", &nodpek(el, 4), &nodpek(el, 5), &nodpek(el, 6), &nodpek(el, 7), &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3));
				break;
			case NO_ELEMENT:
				printf("error scan_avs wrong type of element read\n");
				break;
		}
		id = 2;

		for (i = 0; i < nbf; i++)
			nodpek(el, i)--;elid
		(el) = id;
	}
	elm_type_for_grid = ELMTET1;

	cout << "at the end " << endl;

	fclose(fp);

	checkIfOnlyOneElmTypeInGrid();

}

//===========================================================================


void WavesGridB::scan_2Dgid(const char *file, int elnum1)

{

	FILE *fp;
	int i, j;
	int el;
	int d, id;
	double dd;
	fp = fopen(file, "r");
	//printf("Reading from %s\n", file);

	ifstream isnode(file, ios::in);
	if (!isnode)
	{
		printf("error: read_gid The node file was not found\n");
		return;
	}
	int nnode, elnum;
	fscanf(fp, "%i %i %i %i %i", &nnode, &elnum, &i, &i, &i);
	isnode >> nnode >> elnum >> i >> i >> i;
	//printf(" There are %d nodes and %d elements\n", nnode, elnum);
	maxnne = 3;
	nsd = 2;
	nno = nnode;
	nel = elnum;
	MV_Vector<ElementType> etypes(elnum);
	nbind = 0;

	redim(nsd, nno, nel, maxnne, nbind, NO_ELEMENT);

	if (nsd == 2)
		for (el = 0; el < elnum; el++)
		{
			etypes(el) = ELMTRI1;
			cout << "etypes" << etypes(el) << endl;
		}

	if (nsd == 3)
		for (i = 0; i < nnode; i++)
		{
			fscanf(fp, "%i %lf %lf %lf\n", &j, &coord(i, 0), &coord(i, 1), &coord(i, 2));
		}
	else if (nsd == 2)
		for (i = 0; i < nnode; i++)
		{
			fscanf(fp, "%i %lf %lf %lf\n", &j, &coord(i, 0), &coord(i, 1), &dd);
		}

	cout << "before for elnum" << elnum << endl;

	int nbf = 0;

	for (el = 0; el < elnum1; el++)
	{

		// read number of element
		ElementType etype = etypes(el);
		switch (etype)
		{
			case ELMTET1:
				nbf = 4;
				fscanf(fp, "%i %i %i %i %i %i \n", &d, &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3), &nodpek(el, 0), &id);
				cout << nodpek(el, 0) << "  " << nodpek(el, 1) << "  " << nodpek(el, 2) << "  " << nodpek(el, 3) << endl;
				break;
			case ELMPYR1:
				nbf = 5;
				fscanf(fp, "%i %i %i %i %i\n", &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3), &nodpek(el, 4), &nodpek(el, 0));
				break;
			case ELMHEX1:
				nbf = 8;
				fscanf(fp, "%i %i %i %i %i %i %i %i\n", &nodpek(el, 4), &nodpek(el, 5), &nodpek(el, 6), &nodpek(el, 7), &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3));
				break;
			case ELMPRISM1:
				nbf = 6;
				fscanf(fp, "%i %i %i %i %i %i\n", &nodpek(el, 3), &nodpek(el, 4), &nodpek(el, 5), &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2));
				break;
			case ELMTRI1:
				nbf = 3;
				fscanf(fp, "%i %i %i\n", &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2));
				break;
			case ELMQUAD1:
				nbf = 4;
				fscanf(fp, "%i %i %i %i\n", &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3));
				break;
			case NO_ELEMENT:
				printf("error scan_avs wrong type of element read\n");
				break;
		}
		for (i = 0; i < nbf; i++)
			nodpek(el, i)--;elid
		(el) = id;

	}

	for (el = elnum1; el < elnum; el++)
	{

		// read number of element
		ElementType etype = etypes(el);

		switch (etype)
		{
			case ELMTET1:
				nbf = 4;
				fscanf(fp, "%i %i %i %i %i  \n", &d, &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3), &nodpek(el, 0));
				cout << nodpek(el, 0) << "  " << nodpek(el, 1) << "  " << nodpek(el, 2) << "  " << nodpek(el, 3) << endl;
				break;
			case ELMPYR1:
				nbf = 5;
				fscanf(fp, "%i %i %i %i %i\n", &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3), &nodpek(el, 4), &nodpek(el, 0));
				break;
			case ELMHEX1:
				nbf = 8;
				fscanf(fp, "%i %i %i %i %i %i %i %i\n", &nodpek(el, 4), &nodpek(el, 5), &nodpek(el, 6), &nodpek(el, 7), &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3));
				break;
			case NO_ELEMENT:
				printf("error scan_avs wrong type of element read\n");
				break;
		}
		id = 2;

		for (i = 0; i < nbf; i++)
			nodpek(el, i)--;elid
		(el) = id;
	}
	elm_type_for_grid = ELMTRI1;

	cout << "at the end " << endl;

	fclose(fp);

	checkIfOnlyOneElmTypeInGrid();

}


void WavesGridB::scan_gid_difmat(const char *file)

{

	FILE *fp;
	int i, j;
	int el;
	int d;
	double dd;
	fp = fopen(file, "r");
	//printf("Reading from %s\n", file);

	ifstream isnode(file, ios::in);
	if (!isnode)
	{
		printf("error: read_gid The node file was not found\n");
		return;
	}
	int nnode, elnum;
	fscanf(fp, "%i %i %i %i %i", &nnode, &elnum, &i, &i, &i);
	isnode >> nnode >> elnum >> i >> i >> i;
	//printf(" There are %d nodes and %d elements\n", nnode, elnum);
	maxnne = 4;
	nsd = 3;
	nno = nnode;
	nel = elnum;
	MV_Vector<ElementType> etypes(elnum);
	MV_Vector<int> id(elnum);

	nbind = 0;

	redim(nsd, nno, nel, maxnne, nbind, NO_ELEMENT);

	for (el = 0; el < elnum; el++)
	{
		etypes(el) = ELMTET1;
		cout << "etypes" << etypes(el) << endl;
	}

	if (nsd == 3)
		for (i = 0; i < nnode; i++)
		{
			fscanf(fp, "%i %lf %lf %lf\n", &j, &coord(i, 0), &coord(i, 1), &coord(i, 2));
		}
	else if (nsd == 2)
		for (i = 0; i < nnode; i++)
		{
			fscanf(fp, "%i %lf %lf %lf\n", &j, &coord(i, 0), &coord(i, 1), &dd);
		}

	cout << "before for elnum" << elnum << endl;

	int nbf = 0;

	for (el = 0; el < elnum; el++)
	{

		// read number of element
		ElementType etype = etypes(el);
		switch (etype)
		{
			case ELMTET1:
				nbf = 4;
				fscanf(fp, "%i %i %i %i %i %i \n", &d, &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3), &nodpek(el, 0), &id(el));
				cout << nodpek(el, 0) << "  " << nodpek(el, 1) << "  " << nodpek(el, 2) << "  " << nodpek(el, 3) << endl;
				break;
			case ELMPYR1:
				nbf = 5;
				fscanf(fp, "%i %i %i %i %i\n", &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3), &nodpek(el, 4), &nodpek(el, 0));
				break;
			case ELMHEX1:
				nbf = 8;
				fscanf(fp, "%i %i %i %i %i %i %i %i\n", &nodpek(el, 4), &nodpek(el, 5), &nodpek(el, 6), &nodpek(el, 7), &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3));
				break;
			case ELMPRISM1:
				nbf = 6;
				fscanf(fp, "%i %i %i %i %i %i\n", &nodpek(el, 3), &nodpek(el, 4), &nodpek(el, 5), &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2));
				break;
			case ELMTRI1:
				nbf = 3;
				fscanf(fp, "%i %i %i\n", &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2));
				break;
			case ELMQUAD1:
				nbf = 4;
				fscanf(fp, "%i %i %i %i\n", &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3));
				break;
			case NO_ELEMENT:
				printf("error scan_avs wrong type of element read\n");
				break;
		}
		for (i = 0; i < nbf; i++)
			nodpek(el, i)--;elid
		(el) = id(el);

	}

	elm_type_for_grid = ELMTET1;

	cout << "at the end " << endl;

	fclose(fp);

	checkIfOnlyOneElmTypeInGrid();

}

//=================================================

void WavesGridB::scan_gid(const char *file, bool type_of_mat)

{

	FILE *fp;
	int i, j;
	int el;
	int d, id;
	double dd;
	fp = fopen(file, "r");
	//printf("Reading from %s\n", file);

	ifstream isnode(file, ios::in);
	if (!isnode)
	{
		printf("error: read_gid The node file was not found\n");
		return;
	}
	int nnode, elnum;
	fscanf(fp, "%i %i %i %i %i", &nnode, &elnum, &i, &i, &i);
	isnode >> nnode >> elnum >> i >> i >> i;
	//printf(" There are %d nodes and %d elements\n", nnode, elnum);
	maxnne = 4;
	nsd = 3;
	nno = nnode;
	nel = elnum;
	MV_Vector<ElementType> etypes(elnum);
	nbind = 0;

	redim(nsd, nno, nel, maxnne, nbind, NO_ELEMENT);

	for (el = 0; el < elnum; el++)
	{
		etypes(el) = ELMTET1;
		cout << "etypes" << etypes(el) << endl;
	}

	if (nsd == 3)
		for (i = 0; i < nnode; i++)
		{
			fscanf(fp, "%i %lf %lf %lf\n", &j, &coord(i, 0), &coord(i, 1), &coord(i, 2));
		}
	else if (nsd == 2)
		for (i = 0; i < nnode; i++)
		{
			fscanf(fp, "%i %lf %lf %lf\n", &j, &coord(i, 0), &coord(i, 1), &dd);
		}

	cout << "before for elnum" << elnum << endl;

	int nbf = 0;

	for (el = 0; el < elnum; el++)
	{

		// read number of element
		ElementType etype = etypes(el);
		switch (etype)
		{
			case ELMTET1:
				nbf = 4;
				if (type_of_mat == 1)
					fscanf(fp, "%i %i %i %i %i %i \n", &d, &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3), &nodpek(el, 0), &id);
				cout << nodpek(el, 0) << "  " << nodpek(el, 1) << "  " << nodpek(el, 2) << "  " << nodpek(el, 3) << endl;
				if (type_of_mat == 0)
				{
					fscanf(fp, "%i %i %i %i %i  \n", &d, &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3), &nodpek(el, 0));
					cout << nodpek(el, 0) << "  " << nodpek(el, 1) << "  " << nodpek(el, 2) << "  " << nodpek(el, 3) << endl;
					id = 1;
				}
				break;
			case ELMPYR1:
				nbf = 5;
				fscanf(fp, "%i %i %i %i %i\n", &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3), &nodpek(el, 4), &nodpek(el, 0));
				break;
			case ELMHEX1:
				nbf = 8;
				fscanf(fp, "%i %i %i %i %i %i %i %i\n", &nodpek(el, 4), &nodpek(el, 5), &nodpek(el, 6), &nodpek(el, 7), &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3));
				break;
			case ELMPRISM1:
				nbf = 6;
				fscanf(fp, "%i %i %i %i %i %i\n", &nodpek(el, 3), &nodpek(el, 4), &nodpek(el, 5), &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2));
				break;
			case ELMTRI1:
				nbf = 3;
				fscanf(fp, "%i %i %i\n", &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2));
				break;
			case ELMQUAD1:
				nbf = 4;
				fscanf(fp, "%i %i %i %i\n", &nodpek(el, 0), &nodpek(el, 1), &nodpek(el, 2), &nodpek(el, 3));
				break;
			case NO_ELEMENT:
				printf("error scan_avs wrong type of element read\n");
				break;
		}
		for (i = 0; i < nbf; i++)
			nodpek(el, i)--;elid
		(el) = id;

	}

	elm_type_for_grid = ELMTET1;

	cout << "at the end " << endl;

	fclose(fp);

	checkIfOnlyOneElmTypeInGrid();

}


bool WavesGridB::checkIfOnlyOneElmTypeInGrid()

{
	// check if there is only one element type in the grid.
	// if so it is unnnecessary to store the elementtype in a vector

	if (oneElementTypeInGrid())
		return true;
	ElementType etype = getElementType(0); // get the candidate for single elmtype
	int e;
	for (e = 0; e < nel; e++)
	{
		if (element_type(e) != etype)
			return false;
	}
	// if not yet finished there is only one element type namely etype
	// set the 
	elm_type_for_grid = etype;
	element_type.newsize(0);
	return true;
}

void WavesGridB::print(const char *file) const

{
	FILE *fp;
	int i;
	int el;
	fp = fopen(file, "r");
	//printf("Writing data to: %s\n", file);
	//printf("Nodes: %i, WavesElements: %i\n", nno, nel);
	int nnode = nno;
	int elnum = nel;
	fp = fopen(file, "w");
	fprintf(fp, "%i %i %i 0 0\n", nnode, elnum, 1);

	if (nsd == 3)
	{
		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n ", i + 1, getCoor(i, 0), getCoor(i, 1), getCoor(i, 2));
	}
	else if (nsd == 2)
	{
		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n ", i + 1, getCoor(i, 0), getCoor(i, 1), 0.0);
	}

	for (el = 0; el < elnum; el++)
	{
		switch (getElementType(el))
		{
			case ELMTET1:
				fprintf(fp, "%i %i  tet  %i %i %i %i\n", el + 1, elid(el) + 1, loc2glob(el, 1) + 1, loc2glob(el, 2) + 1, loc2glob(el, 3) + 1, loc2glob(el, 0) + 1);
				break;
			case ELMPYR1:
				fprintf(fp, "%i %i  pyr  %i %i %i %i %i\n", el + 1, elid(el) + 1, loc2glob(el, 1) + 1, loc2glob(el, 2) + 1, loc2glob(el, 3) + 1, loc2glob(el, 4) + 1, loc2glob(el, 0) + 1);
				break;
			case ELMHEX1:
				fprintf(fp, "%i %i  hex  %i %i %i %i %i %i %i %i\n", el + 1, elid(el) + 1, loc2glob(el, 4) + 1, loc2glob(el, 5) + 1, loc2glob(el, 6) + 1, loc2glob(el, 7) + 1, loc2glob(el, 0) + 1, loc2glob(el, 1) + 1, loc2glob(el, 2) + 1, loc2glob(el, 3) + 1);
				break;
			case ELMPRISM1:
				fprintf(fp, "%i %i  prism  %i %i %i %i %i %i\n", el + 1, elid(el) + 1, loc2glob(el, 3) + 1, loc2glob(el, 4) + 1, loc2glob(el, 5) + 1, loc2glob(el, 0) + 1, loc2glob(el, 1) + 1, loc2glob(el, 2) + 1);
				break;
			case ELMTRI1:
				fprintf(fp, "%i %i  tri  %i %i %i\n", el + 1, elid(el) + 1, loc2glob(el, 0) + 1, loc2glob(el, 1) + 1, loc2glob(el, 2) + 1);
				break;
			case ELMQUAD1:
				fprintf(fp, "%i %i  quad  %i %i %i %i\n", el + 1, elid(el) + 1, loc2glob(el, 0) + 1, loc2glob(el, 1) + 1, loc2glob(el, 2) + 1, loc2glob(el, 3) + 1);
				break;
			case NO_ELEMENT:
				printf("error print wrong type of element\n");
				break;
		}
	}
	fprintf(fp, "1 1\n");
	fprintf(fp, " u0, none\n");

	for (i = 0; i < nnode; i++)
	{
		fprintf(fp, "%i %f\n ", i + 1, 1.0);

	}
	fclose(fp);
}


void WavesGridB::print_inp_amira(const char *file) const

{
	FILE *fp;
	int i;
	int el;
	fp = fopen(file, "r");
	//printf("Writing data to: %s\n", file);
	//printf("Nodes: %i, WavesElements: %i\n", nno, nel);
	int nnode = nno;
	int elnum = nel;
	fp = fopen(file, "w");
	fprintf(fp, "%i %i %i 0 0\n", nnode, elnum, 1);

	if (nsd == 3)
	{
		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n ", i + 1, getCoor(i, 0), getCoor(i, 1), getCoor(i, 2));
	}
	else if (nsd == 2)
	{
		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n ", i + 1, getCoor(i, 0), getCoor(i, 1), 0.0);
	}

	for (el = 0; el < elnum; el++)
	{
		switch (getElementType(el))
		{
			case ELMTET1:
				fprintf(fp, "%i %i  tet  %i %i %i %i\n", el + 1, elid(el) + 2, loc2glob(el, 1) + 2, loc2glob(el, 2) + 2, loc2glob(el, 3) + 2, loc2glob(el, 0) + 2);
				break;
			case ELMPYR1:
				fprintf(fp, "%i %i  pyr  %i %i %i %i %i\n", el + 1, elid(el) + 1, loc2glob(el, 1) + 1, loc2glob(el, 2) + 1, loc2glob(el, 3) + 1, loc2glob(el, 4) + 1, loc2glob(el, 0) + 1);
				break;
			case ELMHEX1:
				fprintf(fp, "%i %i  hex  %i %i %i %i %i %i %i %i\n", el + 1, elid(el) + 1, loc2glob(el, 4) + 1, loc2glob(el, 5) + 1, loc2glob(el, 6) + 1, loc2glob(el, 7) + 1, loc2glob(el, 0) + 1, loc2glob(el, 1) + 1, loc2glob(el, 2) + 1, loc2glob(el, 3) + 1);
				break;
			case ELMPRISM1:
				fprintf(fp, "%i %i  prism  %i %i %i %i %i %i\n", el + 1, elid(el) + 1, loc2glob(el, 3) + 1, loc2glob(el, 4) + 1, loc2glob(el, 5) + 1, loc2glob(el, 0) + 1, loc2glob(el, 1) + 1, loc2glob(el, 2) + 1);
				break;
			case ELMTRI1:
				fprintf(fp, "%i %i  tri  %i %i %i\n", el + 1, elid(el) + 1, loc2glob(el, 0) + 1, loc2glob(el, 1) + 1, loc2glob(el, 2) + 1);
				break;
			case ELMQUAD1:
				fprintf(fp, "%i %i  quad  %i %i %i %i\n", el + 1, elid(el) + 1, loc2glob(el, 0) + 1, loc2glob(el, 1) + 1, loc2glob(el, 2) + 1, loc2glob(el, 3) + 1);
				break;
			case NO_ELEMENT:
				printf("error print wrong type of element\n");
				break;
		}
	}
	fprintf(fp, "1 1\n");
	fprintf(fp, " u0, none\n");

	for (i = 0; i < nnode; i++)
	{
		fprintf(fp, "%i %f\n ", i + 1, 1.0);

	}
	fclose(fp);
}


void WavesGridB::print_vel(const char *file, MV_Vector<int>& velocity) const

{
	FILE *fp;
	int i;
	int el;
	fp = fopen(file, "r");
	//printf("Writing data to: %s\n", file);
	//printf("Nodes: %i, WavesElements: %i\n", nno, nel);
	int nnode = nno;
	int elnum = nel;
	fp = fopen(file, "w");
	fprintf(fp, "%i %i %i 0 0\n", nnode, elnum, 1);

	if (nsd == 3)
	{
		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n ", i + 1, getCoor(i, 0), getCoor(i, 1), getCoor(i, 2));
	}
	else if (nsd == 2)
	{
		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n ", i + 1, getCoor(i, 0), getCoor(i, 1), 0.0);
	}

	for (el = 0; el < elnum; el++)
	{
		switch (getElementType(el))
		{
			case ELMTET1:
				fprintf(fp, "%i %i  tet  %i %i %i %i\n", el + 1, elid(el) + 1, loc2glob(el, 1) + 1, loc2glob(el, 0) + 1, loc2glob(el, 3) + 1, loc2glob(el, 2) + 1);

				break;
			case ELMPYR1:
				fprintf(fp, "%i %i  pyr  %i %i %i %i %i\n", el + 1, elid(el) + 1, loc2glob(el, 1) + 1, loc2glob(el, 2) + 1, loc2glob(el, 3) + 1, loc2glob(el, 4) + 1, loc2glob(el, 0) + 1);
				break;
			case ELMHEX1:
				fprintf(fp, "%i %i  hex  %i %i %i %i %i %i %i %i\n", el + 1, elid(el) + 1, loc2glob(el, 4) + 1, loc2glob(el, 5) + 1, loc2glob(el, 6) + 1, loc2glob(el, 7) + 1, loc2glob(el, 0) + 1, loc2glob(el, 1) + 1, loc2glob(el, 2) + 1, loc2glob(el, 3) + 1);
				break;
			case ELMPRISM1:
				fprintf(fp, "%i %i  prism  %i %i %i %i %i %i\n", el + 1, elid(el) + 1, loc2glob(el, 3) + 1, loc2glob(el, 4) + 1, loc2glob(el, 5) + 1, loc2glob(el, 0) + 1, loc2glob(el, 1) + 1, loc2glob(el, 2) + 1);
				break;
			case ELMTRI1:
				fprintf(fp, "%i %i  tri  %i %i %i\n", el + 1, elid(el) + 1, loc2glob(el, 0) + 1, loc2glob(el, 1) + 1, loc2glob(el, 2) + 1);
				break;
			case ELMQUAD1:
				fprintf(fp, "%i %i  quad  %i %i %i %i\n", el + 1, elid(el) + 1, loc2glob(el, 0) + 1, loc2glob(el, 1) + 1, loc2glob(el, 2) + 1, loc2glob(el, 3) + 1);
				break;
			case NO_ELEMENT:
				printf("error print wrong type of element\n");
				break;
		}
	}
	fprintf(fp, "1 1");
	fprintf(fp, "\n");
	fprintf(fp, "u , Unit\n");
	for (i = 0; i < getNoNodes(); i++)
	{
		fprintf(fp, "%i ", i + 1);
		fprintf(fp, "%i \n", velocity(i));
	}
	fclose(fp);
}

int WavesGridB::getNoNodesForElmType(const ElementType etype) const
{
	switch (etype)
	{
		case NO_ELEMENT:
			return 0;
		case ELMTRI1:
			return 3;
		case ELMQUAD1:
			return 4;
		case ELMTET1:
			return 4;
		case ELMPYR1:
			return 5;
		case ELMPRISM1:
			return 6;
		case ELMHEX1:
			return 8;
	}
	return 0;
}

int WavesGridB::Coord2Node(double x, double y, double z) const
{
	int i, node;
	double eps = 0.0000001;
	node = 0;

	if (nsd == 2)
	{
		for (i = 0; i < getNoNodes(); i++)
		{
			double x_i = getCoor(i, 0);
			double y_i = getCoor(i, 1);

			if (fabs(x - x_i) < eps && fabs(y - y_i) < eps)
			{
				node = i;
				break;
			}
		}
	}
	else if (nsd == 3)
	{
		for (i = 0; i < getNoNodes(); i++)
		{
			double x_i = getCoor(i, 0);
			double y_i = getCoor(i, 1);
			double z_i = getCoor(i, 2);

			if (fabs(x - x_i) < eps && fabs(y - y_i) < eps && fabs(z - z_i) < eps)
			{
				node = i;
				break;
			}
		}
	}
	return node;

}

