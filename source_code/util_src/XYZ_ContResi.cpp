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



//========= XYZ format for point-cloud ==========//
/*
    1 MH  -30.302   30.594   73.150
    1 MH  -30.113   30.371   71.618
    1 MH  -27.405   27.894   71.199
    1 MH  -28.776   28.421   70.605
    1 MH  -27.541   29.403   70.738
    2 IN  -29.448   33.231   69.748
    2 IC  -29.800   33.305   68.344
    2 IC  -29.435   31.994   67.649
    2 IO  -28.554   31.259   68.126
    2 IC  -29.064   34.526   67.741
*/

//--------- load XYZ ----------//
void Load_XYZ(string &fn,vector <vector <vector <double> > > &xyz,
	vector <vector <string> > &str, vector <vector <string> > &remain,
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
	xyz.clear();
	str.clear();
	resi.clear();
	remain.clear();
	vector <double> point(3);
	vector <vector <double> > xyz_tmp;
	vector <string> str_tmp;
	vector <string> remain_tmp;
	string prev="";
	string str_rec;
	int first=1;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		www>>temp;
		if(first==1)
		{
			first=0;
			prev=temp;
		}
		if(temp!=prev)
		{
			xyz.push_back(xyz_tmp);
			str.push_back(str_tmp);
			remain.push_back(remain_tmp);
			resi.push_back(prev);
			xyz_tmp.clear();
			str_tmp.clear();
			remain_tmp.clear();
			prev=temp;
		}
		www>>temp;
		str_tmp.push_back(temp);
		www>>point[0]>>point[1]>>point[2];
		xyz_tmp.push_back(point);
		//remain
		string remain_rec="";
		for(;;)
		{
			if(! (www>>temp) )break;
			remain_rec=remain_rec+temp+" ";
		}
		remain_tmp.push_back(remain_rec);
	}
	//termi
	if(first==0)
	{
		xyz.push_back(xyz_tmp);
		str.push_back(str_tmp);
		remain.push_back(remain_tmp);
		resi.push_back(prev);
	}
}

//------ get distance -----//
double Distance2(vector <double> &in1,vector <double> &in2)
{
	double dist2=0;
	for(int i=0;i<3;i++)dist2+=(in1[i]-in2[i])*(in1[i]-in2[i]);
	return dist2;
}
double Distance(vector <double> &in1,vector <double> &in2)
{
	double dist=0;
	for(int i=0;i<3;i++)dist+=(in1[i]-in2[i])*(in1[i]-in2[i]);
	return sqrt(1.0*dist);
}
//------- check distance between two moleculae -----//
double TwoMol_Distance2(double thres,
	vector <vector <double> >  &in1,
	vector <vector <double> >  &in2,
	vector <int> &atom1,
	vector <int> &atom2)
{
	//init
	atom1.clear();
	atom2.clear();
	//proc
	double mindist2=9999999.0*9999999.0;
	double thres2=thres*thres;
	int found=0;
	vector <int> atom1_;
	vector <int> atom2_;
	for(int i=0;i<(int)in1.size();i++)
	{
		for(int j=0;j<(int)in2.size();j++)
		{
			double dist2=Distance2(in1[i],in2[j]);
			if(dist2<mindist2)mindist2=dist2;
			if(dist2<thres2)
			{
				atom1_.push_back(i);
				atom2_.push_back(j);
				found=1;
			}
		}
	}
	//final
	if(found==1)
	{
		atom1.resize((int)in1.size(),-1);
		for(int k=0;k<(int)atom1_.size();k++)atom1[atom1_[k]]=atom2_[k];
		atom2.resize((int)in2.size(),-1);
		for(int k=0;k<(int)atom2_.size();k++)atom2[atom2_[k]]=atom1_[k];
	}
	//return
	return mindist2;
}
double TwoMol_Distance(double thres,
	vector <vector <double> >  &in1,
	vector <vector <double> >  &in2,
	vector <int> &atom1,
	vector <int> &atom2)
{
	//init
	atom1.clear();
	atom2.clear();
	//proc
	double mindist=9999999.0;
	int found=0;
	vector <int> atom1_;
	vector <int> atom2_;
	for(int i=0;i<(int)in1.size();i++)
	{
		for(int j=0;j<(int)in2.size();j++)
		{
			double dist=Distance(in1[i],in2[j]);
			if(dist<mindist)mindist=dist;
			if(dist<thres)
			{
				atom1_.push_back(i);
				atom2_.push_back(j);
				found=1;
			}
		}
	}
	//final
	if(found==1)
	{
		atom1.resize((int)in1.size(),-1);
		for(int k=0;k<(int)atom1_.size();k++)atom1[atom1_[k]]=atom2_[k];
		atom2.resize((int)in2.size(),-1);
		for(int k=0;k<(int)atom2_.size();k++)atom2[atom2_[k]]=atom1_[k];
	}
	//return
	return mindist;
}

//------ search binding residues ----------//
void Search_Binding_Residue(double thres,
	vector <vector <vector <double> > > &xyz1,
	vector <vector <vector <double> > > &xyz2,
	vector <pair<int,int> > &binding_resi, 
	vector <int> &binding_resi1,
	vector <int> &binding_resi2,
	vector <vector <int> > &binding_lig1,
	vector <vector <int> > &binding_lig2,
	vector <vector <vector <int> > > &binding_rec1,
	vector <vector <vector <int> > > &binding_rec2)
{
	binding_resi.clear();
	binding_resi1.resize((int)xyz1.size(),0);
	binding_resi2.resize((int)xyz2.size(),0);
	binding_lig1.resize((int)xyz1.size());
	binding_lig2.resize((int)xyz2.size());
	binding_rec1.resize((int)xyz1.size());
	binding_rec2.resize((int)xyz2.size());
	double thres2=thres*thres;
	int count=0;
	for(int i=0;i<(int)xyz1.size();i++)
		for(int j=0;j<(int)xyz2.size();j++)
		{
			vector <int> atom1;
			vector <int> atom2;
			double retdist2=TwoMol_Distance2(thres,xyz1[i],xyz2[j],atom1,atom2);
			if(retdist2<thres2)
			{
				binding_resi.push_back(pair<int,int>(i,j));
				binding_rec1[i].push_back(atom1);
				binding_rec2[j].push_back(atom2);
				binding_lig1[i].push_back(j);
				binding_lig2[j].push_back(i);
				binding_resi1[i]=1;
				binding_resi2[j]=1;
			}
		}
}

//--------- output to PointCloud XYZ format --------//
//-> example
/*
   24 FC  -16.868   31.976   78.613 P   6
   24 FC  -15.683   33.130   76.402 0   6
   24 FC  -16.188   33.184   78.720 P   6
   24 FC  -15.594   33.762   77.616 P   6
   25 RN  -19.819   27.596   76.450 P   7
   25 RC  -20.075   26.333   75.734 0   7
   25 RC  -21.358   26.330   74.879 0   7
   25 RO  -21.401   25.695   73.836 0   7
   25 RC  -20.113   25.128   76.685 0   7
   25 RC  -19.951   23.799   75.949 0   7

*/

//--- atom_level ----//
void Output_PointCloud_AtomLevel(FILE *fp, 	
	vector <vector <vector <double> > > &xyz,
	vector <vector <string> > &str,	
	vector <vector <string> > &remain,
	vector <string> &resi,
	vector <vector <int> > &binding_lig,
	vector <vector <vector <int> > > &binding_rec)
{
	for(int i=0;i<(int)xyz.size();i++)
	{
		//check binding
		int size=(int)binding_rec[i].size();
		int len=(int)xyz[i].size();
		vector <int> atom_rec(len,0);
		if(size>0)  //-> contain binding atoms for this residue
		{
			for(int k=0;k<size;k++)
			{
				int label=binding_lig[i][k];
				for(int j=0;j<len;j++)
				{
					if(binding_rec[i][k][j]>=0)atom_rec[j]=label+1;
				}
			}
		}
		//output to file
		for(int j=0;j<len;j++)
		{
			fprintf(fp,"%5s %3s %8.3f %8.3f %8.3f %1d %s\n",
				resi[i].c_str(),str[i][j].c_str(),
				xyz[i][j][0],xyz[i][j][1],xyz[i][j][2],
				atom_rec[j],remain[i][j].c_str());
		}
	}
}

//--- resi_level ----//
void Output_PointCloud_ResiLevel(FILE *fp, 	
	vector <vector <vector <double> > > &xyz,
	vector <vector <string> > &str,	
	vector <vector <string> > &remain,
	vector <string> &resi,
	vector <int> &binding_rec)
{
	for(int i=0;i<(int)xyz.size();i++)
	{
		//output to file
		for(int j=0;j<(int)xyz[i].size();j++)
		{
			fprintf(fp,"%5s %3s %8.3f %8.3f %8.3f %1d %s\n",
				resi[i].c_str(),str[i][j].c_str(),
				xyz[i][j][0],xyz[i][j][1],xyz[i][j][2],
				binding_rec[i],remain[i][j].c_str());
		}
	}
}

//========== XYZ_ContResi ===========//
//-> [example]: thres could be set to 6.5A according to scPDB
//-> [note]: we ALWAYS fix the first input as protein
void XYZ_ContResi(string &protein_xyz, string &ligand_xyz, double thres,
	string &protein_out, int RESI_or_ATOM)
{
	//---- load protein -----//
	vector <vector <vector <double> > > xyz1;
	vector <vector <string> > str1;
	vector <vector <string> > remain1;
	vector <string> resi1;
	Load_XYZ(protein_xyz,xyz1,str1,remain1,resi1);
	//---- load ligand -----//
	vector <vector <vector <double> > > xyz2;
	vector <vector <string> > str2;
	vector <vector <string> > remain2;
	vector <string> resi2;
	Load_XYZ(ligand_xyz,xyz2,str2,remain2,resi2);
	//---- calculate binding residues ----//
	vector <pair<int,int> > binding_resi;
	vector <int> binding_resi1;
	vector <int> binding_resi2;
	vector <vector <int> > binding_lig1;
	vector <vector <int> > binding_lig2;
	vector <vector <vector <int> > > binding_rec1;
	vector <vector <vector <int> > > binding_rec2;
	Search_Binding_Residue(thres,xyz1,xyz2,binding_resi,
		binding_resi1,binding_resi2,
		binding_lig1,binding_lig2,
		binding_rec1,binding_rec2);
	//---- output to file -----//
	FILE *fp=fopen(protein_out.c_str(),"wb");
	if(RESI_or_ATOM==0) //-> output at resi-level
	{
		Output_PointCloud_ResiLevel(fp,xyz1,str1,remain1,resi1,binding_resi1);
	}
	else                //-> output at atom-level
	{
		Output_PointCloud_AtomLevel(fp,xyz1,str1,remain1,resi1,binding_lig1,binding_rec1);
	}
	//---- screenout binding residues -----//
	for(int i=0;i<(int)binding_resi1.size();i++)
	{
		if(binding_resi1[i]==1)printf("%c|%s\n",str1[i][0][0],resi1[i].c_str());
	}
}



//-------- main ---------//
int main(int argc,char **argv)
{
	//------ XYZ_ContResi -------//
	{
		if(argc<6)
		{
			fprintf(stderr,"XYZ_ContResi <protein_xyz> <ligand_xyz> <thres> <output_xyz> <RESI_or_ATOM> \n");
			fprintf(stderr,"[note]: RESI_or_ATOM -> 0 for residue_lvel and 1 for atom_level \n");
			exit(-1);
		}
		string protein_xyz=argv[1];
		string ligand_xyz=argv[2];
		double thres=atof(argv[3]);
		string output_xyz=argv[4];
		int RESI_or_ATOM=atoi(argv[5]);
		//proc
		XYZ_ContResi(protein_xyz,ligand_xyz,thres,output_xyz,RESI_or_ATOM);
		//exit
		exit(0);
	}
}


