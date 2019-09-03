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


//----------- AA three to one ---------//
char WWW_Three2One_III(const char *input)
{
	int i;
	int len;
	int result;
	//encoding
	len=(int)strlen(input);
	if(len!=3)return 'X';
	result=0;
	for(i=0;i<len;i++)result+=(input[i]-'A')*(int)pow(1.0*26,1.0*i);
	//switch
	switch(result)
	{
		case 286:return 'A';
		case 4498:return 'R';
		case 9256:return 'N';
		case 10608:return 'D';
		case 12794:return 'C';
		case 9080:return 'Q';
		case 13812:return 'E';
		case 16516:return 'G';
		case 12383:return 'H';
		case 2998:return 'I';
		case 13635:return 'L';
		case 12803:return 'K';
		case 12960:return 'M';
		case 2901:return 'F';
		case 9921:return 'P';
		case 11614:return 'S';
		case 11693:return 'T';
		case 10601:return 'W';
		case 12135:return 'Y';
		case 7457:return 'V';
		default:return 'X';
	}
}

//========= XYZ format for point-cloud ==========//
/*
    1 MN  -30.766   33.288   73.186 0  43
    1 MC  -29.721   33.119   72.137 0  43
    1 MC  -30.322   33.150   70.741 0  43
    1 MO  -31.535   33.139   70.580 0  43
    1 MC  -28.881   31.843   72.350 0  43
    1 MC  -29.636   30.528   72.448 0  43
    1 MS  -28.543   29.109   72.790 0  43
    1 MC  -28.005   28.655   71.145 0  43
    2 IN  -29.448   33.231   69.748 0  38
    2 IC  -29.800   33.305   68.344 0  38
*/

//--------- PDB_To_XYZ_fix ----------//
void PDB_To_XYZ_fix(string &pdb,FILE *fp,int resi_type)
{
	ifstream fin;
	string buf,temp,name;
	//read
	fin.open(pdb.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list %s not found!!\n",pdb.c_str());
		exit(-1);
	}
	int len;
	int resi;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len<3)continue;
		//check TER
		temp=buf.substr(0,3);
		if(temp=="TER"||temp=="END")break;
		//check ATOM
		if(len<4)continue;
		temp=buf.substr(0,4);
		if(temp!="ATOM"&&temp!="HETA")continue;
		//ami
		temp=buf.substr(17,3);
		char c=WWW_Three2One_III(temp.c_str());
		char h=buf[77];
		//resi
		if(resi_type<0)
		{
			temp=buf.substr(22,4);
			resi=atoi(temp.c_str());
		}
		else
		{
			resi=resi_type;
		}
		//mol
		temp=buf.substr(30,8);
		double x=atof(temp.c_str());
		temp=buf.substr(38,8);
		double y=atof(temp.c_str());
		temp=buf.substr(46,8);
		double z=atof(temp.c_str());
		//output
		fprintf(fp,"%5d %c%c %8.3f %8.3f %8.3f \n",
			resi,c,h,x,y,z);
	}
}

//-------- main ---------//
int main(int argc,char **argv)
{
	//------ PDB_To_XYZ_fix -------//
	{
		if(argc<4)
		{
			fprintf(stderr,"PDB_To_XYZ_fix <pdb_file> <xyz_file> <fix_resi> \n");
			fprintf(stderr,"[note]: set fix_resi -1 to use original residue numbering \n");
			exit(-1);
		}
		string pdb_file=argv[1];
		string xyz_file=argv[2];
		int fix_resi=atoi(argv[3]);
		//proc
		FILE *fp=fopen(xyz_file.c_str(),"wb");
		PDB_To_XYZ_fix(pdb_file,fp,fix_resi);
		fclose(fp);
		//exit
		exit(0);
	}
}

