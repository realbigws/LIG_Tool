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

//------- int to string --------//
template <typename T>
string NumberToString( T Number )
{
	ostringstream ss;
	ss << Number;
	return ss.str();
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
//----- check insert_code ------//
int Check_Ins(string &in)
{
	int i=(int)in.length()-1;
	if(in[i]>='0'&&in[i]<='9')return 0;
	else return 1;
}
//--------- load XYZ ----------//
void Load_XYZ_Label(string &fn,vector <vector <vector <double> > > &xyz,
	vector <vector <string> > &str, vector <vector <string> > &lab,
	vector <vector <string> > &remain, vector <string> &resi)
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
	lab.clear();
	remain.clear();
	resi.clear();
	vector <double> point(3);
	vector <vector <double> > xyz_tmp;
	vector <string> str_tmp;
	vector <string> lab_tmp;
	vector <string> remain_tmp;
	string prev="";
	string str_rec;
	int first=1;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		www>>temp;
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
			xyz.push_back(xyz_tmp);
			str.push_back(str_tmp);
			lab.push_back(lab_tmp);
			remain.push_back(remain_tmp);
			resi.push_back(prev);
			xyz_tmp.clear();
			str_tmp.clear();
			lab_tmp.clear();
			remain_tmp.clear();
			prev=temp;
		}
		www>>temp;
		str_tmp.push_back(temp);
		www>>point[0]>>point[1]>>point[2];
		xyz_tmp.push_back(point);
		//remain
		string remain_rec="";
		int count=0;
		string lab_value="0";
		for(;;)
		{
			if(! (www>>temp) )break;
			if(count>0)remain_rec=remain_rec+temp+" ";
			else lab_value=temp;
			count++;
		}
		lab_tmp.push_back(lab_value);
		remain_tmp.push_back(remain_rec);
	}
	//termi
	if(first==0)
	{
		xyz.push_back(xyz_tmp);
		str.push_back(str_tmp);
		lab.push_back(lab_tmp);
		remain.push_back(remain_tmp);
		resi.push_back(prev);
	}
}

//------- output XYZ ---------//
void Output_XYZ_ResiStart(FILE *fp,
	vector <vector <vector <double> > > &xyz,
	vector <vector <string> > &str, vector <vector <string> > &lab,
	vector <vector <string> > &remain, int resi_start)
{
	for(int i=0;i<(int)xyz.size();i++)
	{
		for(int j=0;j<(int)xyz[i].size();j++)
		{
			//output to file
			fprintf(fp,"%5d %3s %8.3f %8.3f %8.3f %s %s\n",
				i+resi_start,str[i][j].c_str(),
				xyz[i][j][0],xyz[i][j][1],xyz[i][j][2],
				lab[i][j].c_str(),remain[i][j].c_str());
		}
	}
}


//-------- main ---------//
int main(int argc,char **argv)
{
	//------ XYZ_ResiStart -------//
	{
		if(argc<4)
		{
			fprintf(stderr,"XYZ_ResiStart <xyz_in> <xyz_out> <resi_start> \n");
			fprintf(stderr,"[note]: set resi_start (1-base); \n");
			exit(-1);
		}
		string xyz_in=argv[1];
		string xyz_out=argv[2];
		int resi_start=atoi(argv[3]);
		//input
		vector <vector <vector <double> > > xyz;
		vector <vector <string> > str;
		vector <vector <string> > lab;
		vector <vector <string> > remain;
		vector <string> resi;
		Load_XYZ_Label(xyz_in,xyz,str,lab,remain,resi);
		//output
		FILE *fp=fopen(xyz_out.c_str(),"wb");
		Output_XYZ_ResiStart(fp,xyz,str,lab,remain,resi_start);
		//exit
		exit(0);
	}
}

