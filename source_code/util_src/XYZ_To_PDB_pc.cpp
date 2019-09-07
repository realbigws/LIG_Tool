#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;


//-------- utility ------//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}


//--------- Parse_Str_Str --------//
int Parse_Str_Str(string &in,vector <string> &out, char separator)
{
	istringstream www(in);
	out.clear();
	int count=0;
	for(;;)
	{
		string buf;
		if(!getline(www,buf,separator))break;
		out.push_back(buf);
		count++;
	}
	return count;
}
//--------- Parse_Str_Str (automatic) --------//
int Parse_Str_Str(string &in,vector <string> &out)
{
	istringstream www(in);
	out.clear();
	int count=0;
	for(;;)
	{
		string buf;
		if(! (www>>buf) )break;
		out.push_back(buf);
		count++;
	}
	return count;
}

//----- check insert_code ------//
int Check_Ins(string &in)
{
	int i=(int)in.length()-1;
	if(in[i]>='0'&&in[i]<='9')return 0;
	else return 1;
}


//========= XYZ format for point-cloud ==========//
/*
    6 DOg  -22.389   29.406   56.176 0  58
    6 DOh  -23.720   31.160   56.319 0  58
    7 KNa  -22.679   27.076   60.272 0  55
    7 KCb  -21.836   26.578   61.340 0  55
    7 KCc  -20.356   26.635   60.950 0  55
    7 KOd  -19.746   25.631   60.582 0  55
    7 KCe  -22.307   25.189   61.783 0  55
    7 KCf  -23.751   25.237   62.298 0  55
    7 KCg  -24.043   24.207   63.371 0  55
    7 KCh  -25.357   24.538   64.095 0  55
    7 KNi  -25.610   23.600   65.237 0  55
    8 QNa  -19.785   27.829   61.073 0  24
    8 QCb  -18.384   28.110   60.732 0  24
    8 QCc  -17.293   27.383   61.502 0  24
    8 QOd  -16.145   27.388   61.063 0  24
...
*/

//--------- load XYZ ----------//
int Load_XYZ(string &fn,
	vector <vector <vector <double> > > &xyz,
	vector <vector <string> > &str, vector <vector <int> > &lab, 
	vector <string> &resi)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(fn.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list %s not found!!\n",fn.c_str());
		exit(-1);
	}
	resi.clear();
	xyz.clear();
	str.clear();
	lab.clear();
	vector <double> point(3);
	vector <vector <double> > xyz_tmp;
	vector <string> str_tmp;
	vector <int> lab_tmp;
	string prev="";
	string str_rec;
	int first=1;
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		//-> get input
		vector <string> tmp_rec;
		int retv=Parse_Str_Str(buf,tmp_rec);
		if(retv<5)
		{
			fprintf(stderr,"xyz %s format bad!!\n",fn.c_str());
			exit(-1);
		}
		temp=tmp_rec[0];
		//check ins_code
		if(Check_Ins(temp)!=1)temp.push_back(' ');
		//record first
		if(first==1)
		{
			first=0;
			prev=temp;
		}
		if(temp!=prev)
		{
			resi.push_back(prev);
			xyz.push_back(xyz_tmp);
			str.push_back(str_tmp);
			lab.push_back(lab_tmp);
			xyz_tmp.clear();
			str_tmp.clear();
			lab_tmp.clear();
			prev=temp;
			count++;
		}
		temp=tmp_rec[1];
		str_tmp.push_back(temp);
		point[0]=atof(tmp_rec[2].c_str());
		point[1]=atof(tmp_rec[3].c_str());
		point[2]=atof(tmp_rec[4].c_str());
		xyz_tmp.push_back(point);
		//-> addi
		int label=0;
		if(retv>5)label=atoi(tmp_rec[5].c_str());
		lab_tmp.push_back(label);
	}
	//termi
	if(first==0)
	{
		resi.push_back(prev);
		xyz.push_back(xyz_tmp);
		str.push_back(str_tmp);
		lab.push_back(lab_tmp);
		count++;
	}
	//return
	return count;
}


//=========== PDB relevant atom labels ===========//
const char *WWW_backbone_atom_name_decode(int pos)
{
	switch (pos)
	{
		case 0:return "N  ";
		case 1:return "CA ";
		case 2:return "C  ";
		case 3:return "O  ";
		default:return 0;
	}
}
const char *WWW_sidechain_atom_name_decode(int pos, char amino)
{
	switch(amino)
	{
		case 'A':
		{
			switch(pos)
			{
				case 0:return "CB ";
				default:return 0;
			}
		}
		case 'R':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD ";
				case 3:return "NE ";
				case 4:return "CZ ";
				case 5:return "NH1";
				case 6:return "NH2";
				default:return 0;
			}
		}
		case 'N':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "OD1";
				case 3:return "ND2";
				default:return 0;
			}
		}
		case 'D':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "OD1";
				case 3:return "OD2";
				default:return 0;
			}
		}
		case 'C':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "SG ";
				default:return 0;
			}
		}
		case 'Q':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD ";
				case 3:return "OE1";
				case 4:return "NE2";
				default:return 0;
			}
		}
		case 'E':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD ";
				case 3:return "OE1";
				case 4:return "OE2";
				default:return 0;
			}
		}
		case 'G':
		{
			switch(pos)   //note: pseudo CB (=CA)
			{
				case 0:return "CB ";
				default:return 0;
			}
		}
		case 'H':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "ND1";
				case 3:return "CD2";
				case 4:return "CE1";
				case 5:return "NE2";
				default:return 0;
			}
		}
		case 'I':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG1";
				case 2:return "CG2";
				case 3:return "CD1";
				default:return 0;
			}
		}
		case 'L':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD1";
				case 3:return "CD2";
				default:return 0;
			}
		}
		case 'K':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD ";
				case 3:return "CE ";
				case 4:return "NZ ";
				default:return 0;
			}
		}
		case 'M':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "SD ";
				case 3:return "CE ";
				default:return 0;
			}
		}
		case 'F':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD1";
				case 3:return "CD2";
				case 4:return "CE1";
				case 5:return "CE2";
				case 6:return "CZ ";
				default:return 0;
			}
		}
		case 'P':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD ";
				default:return 0;
			}
		}
		case 'S':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "OG ";
				default:return 0;
			}
		}
		case 'T':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "OG1";
				case 2:return "CG2";
				default:return 0;
			}
		}
		case 'W':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD1";
				case 3:return "CD2";
				case 4:return "NE1";
				case 5:return "CE2";
				case 6:return "CE3";
				case 7:return "CZ2";
				case 8:return "CZ3";
				case 9:return "CH2";
				default:return 0;
			}
		}
		case 'Y':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG ";
				case 2:return "CD1";
				case 3:return "CD2";
				case 4:return "CE1";
				case 5:return "CE2";
				case 6:return "CZ ";
				case 7:return "OH ";
				default:return 0;
			}
		}
		case 'V':
		{
			switch(pos)
			{
				case 0:return "CB ";
				case 1:return "CG1";
				case 2:return "CG2";
				default:return 0;
			}
		}
		case 'X':  //only consider CB!! (alanine model)
		{
			switch(pos)
			{
				case 0:return "CB ";
				default:return 0;
			}
		}
		default:return 0;
	}
}

//-------- one2three decode ------------//
const char* WWW_One2Three_III(char c)
{
	//switch
	switch(c)
	{
		case 'A':return "ALA";
		case 'R':return "ARG";
		case 'N':return "ASN";
		case 'D':return "ASP";
		case 'C':return "CYS";
		case 'Q':return "GLN";
		case 'E':return "GLU";
		case 'G':return "GLY";
		case 'H':return "HIS";
		case 'I':return "ILE";
		case 'L':return "LEU";
		case 'K':return "LYS";
		case 'M':return "MET";
		case 'F':return "PHE";
		case 'P':return "PRO";
		case 'S':return "SER";
		case 'T':return "THR";
		case 'W':return "TRP";
		case 'Y':return "TYR";
		case 'V':return "VAL";
		default:return "UNK";
	}
}

//------------------ from XYZ to PDB ----------------//
//-> PDB format:
/*
ATOM    496  N   MET A   1     -30.766  33.288  73.186  1.00 19.31      A    N
ATOM    497  CA  MET A   1     -29.721  33.119  72.137  1.00 18.15      A    C
ATOM    498  C   MET A   1     -30.322  33.150  70.741  1.00 19.44      A    C
ATOM    499  O   MET A   1     -31.535  33.139  70.580  1.00 20.69      A    O
ATOM    500  CB  MET A   1     -28.881  31.843  72.350  1.00 16.43      A    C
ATOM    501  CG  MET A   1     -29.636  30.528  72.448  1.00 17.15      A    C
...
*/
void XYZ_To_PDB_pointcloud(string &xyz_file, string &pdb_file, int chain)
{
	//--- load XYZ ---//
	vector <vector <vector <double> > > xyz;
	vector <vector <string> > str;
	vector <vector <int> > lab;
	vector <string> resi;
	int moln=Load_XYZ(xyz_file,xyz,str,lab,resi);
	//--- output PDB ----//
	double rfactor=1;
	int atom_num=1;
	char rel_chain;
	if(chain<=0)rel_chain=' ';
	else rel_chain=chain;
	FILE *fp=fopen(pdb_file.c_str(),"wb");
	for(int i=0;i<moln;i++)
	{
		const char* resi_name;
		const char* atom_name;
		for(int j=0;j<(int)xyz[i].size();j++)
		{
			resi_name=WWW_One2Three_III(str[i][j][0]);
			int atom_lab=str[i][j][2]-'a';
			if(atom_lab<4) atom_name=WWW_backbone_atom_name_decode(atom_lab);
			else atom_name=WWW_sidechain_atom_name_decode(atom_lab-4, str[i][j][0]);
			fprintf(fp,"ATOM  %5d %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f      %c    %c\n",
				atom_num,atom_name,resi_name,rel_chain,resi[i].c_str(),
				xyz[i][j][0],xyz[i][j][1],xyz[i][j][2],rfactor,
				10.0*lab[i][j],rel_chain,str[i][j][1]);
			atom_num++;
		}
	}
	fclose(fp);
}


//---------- main ----------//
int main(int argc,char **argv)
{
	//---- XYZ_To_PDB_pc ----//
	{
		if(argc<4)
		{
			fprintf(stderr,"XYZ_To_PDB_pc <xyz_file> <pdb_file> <chain> \n");
			fprintf(stderr,"[note]: set chain -1 to use 'blank' chain identifier \n");
			exit(-1);
		}
		string xyz_file=argv[1];
		string pdb_file=argv[2];
		int chain=argv[3][0];
		//proc
		XYZ_To_PDB_pointcloud(xyz_file,pdb_file,chain);
		//exit
		exit(0);
	}
}

