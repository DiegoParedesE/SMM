#include <bits/stdc++.h>
using namespace std;

typedef double db;
int n;//Number of Motors
const db pi=3.1415926536;
const db pi4=3.1415926536/4;
const db RU = 8314.462618; //Universal gas constant (J/(K*M))
const int sim_i=1000;
const int sim_i2=sim_i*0.1;
double sim[sim_i+sim_i2][40];
string prop;
//Static int Motor class counter

class Propellant{
public:
    string name; //Propellant name
    db p_den; //Ideal Density (g/cm3)
    db k2ph;//Ratio of specific heats, 2-ph
    db k;//Ratio of specific heats, mixture
    db p_M;//Effective molecular wt (g/mol)
    db p_chT;//Chamber temperature (K)
    db p_dr; //Real Density Ratio
    db p_ce;//Combustion Efficiency
    db pol[4][10];//Polynomial coefficients for Least Squares Fit of Kn as function of Pressure data (a,b,c,d,e,f,g,pMin,pMax)
    int br_n; //Number of burn rate intervals
    db br[7][6]; //Burn rate interval values: a, n, Min Pressure (MPa), Max Pressure (MPa)
    //calculated
    db p_rden; //Real density (g/cm3)
    db R; //Specific gas constant
    db ct;//real combustion temperature


    void SetProp(db p_den, db k2ph, db k, db p_M, db p_chT, db p_dr, db p_ce){
        this->p_den = p_den;
        this->k2ph = k2ph;
        this->k = k;
        this->p_M = p_M;
        this->p_chT = p_chT;
        this->p_dr = p_dr;
        this->p_ce = p_ce;
        p_rden=p_dr*p_den;
        R=RU/p_M;
        ct=p_chT*p_ce;
    }

    void SetBurn(int br_n,db b[7][6]){
        for(int i=0;i<7;i++){
            for(int i2=0;i2<6;i2++){
                br[i][i2]=b[i][i2];
            }
        }
        this->br_n = br_n;
    }

    void SetPol(db pol[4][10]){
        for(int i=0;i<5;i++){
            for(int i2=0;i2<11;i2++){
                this->pol[i][i2]=pol[i][i2];
            }
        }
    }
    double GetRP(db P){
        for(int i=2;i<7;i++){
            if(P<br[i][4]){
                return br[i][1]*pow(P,br[i][2]);
                break;
            }
        }
    }
};

class Motor:public Propellant{
public:


    db c_di; //Combustion inner length (mm)
    db c_l; //Combustion chamber length (mm)
    db g_do; //Grain Outer diameter (initial) (mm)
    db g_dc; //Grain Core diameter (initial) (mm)
    db g_sl; //Grain Segment length (initial) (mm)
    int g_n; //Grain Number of segments
    bool i_os;//Outer surface inhibition
    bool i_is;//Core surface inhibition
    bool i_es;//Ends surface inhibition
    db n_tdi;//Nozzle throat Diameter Initial
    //Calculated
    db c_LD; //Chamber Length to Diameter Ratio
    db g_l; //Grain length (initial)
    db c_v; //Combustion chamber volume
    db g_v;//Grain volume (initial)
    db g_m;//Grain mass (initial)
    db m_vlm;//Volumetric loading fraction
    db p_max;// Maximum allowed pressure (MPa)
    db br_ai;// Initial Burn Area (mm^2)
    db kn_max;//Maximum Kn
    db x_inc;//mm
    db c_cs;//Combustion chamber cross-section
    db n_ef;//Nozzle efficiency
    Propellant act;

    void SetData(db c_di, db c_l, db g_do, db g_dc, db g_sl, int g_n, bool i_os, bool i_is, bool i_es,db n_ef,Propellant act){
        this->c_l = c_l;
        this->c_di = c_di;
        this->g_do = g_do;
        this->g_dc = g_dc;
        this->g_sl = g_sl;
        this->g_n = g_n;
        this->i_os = i_os;
        this->i_is = i_is;
        this->i_es = i_es;
        this->n_ef = n_ef;
        this->act = act;
        g_l=g_sl*g_n;
        c_v=c_di*c_di*c_l*pi4/1000000000;
        g_v=g_do*g_do*2*g_l*pi4-g_dc*g_dc*2*g_l*pi4;
        g_m=g_v/1000*act.p_rden/1000;
        m_vlm=g_v/c_v;
        c_LD=c_l/c_di;
        c_cs=c_di*c_di*pi4;

    }

    void KnMax(db pmax){
         //Calculating Kn max
        for(int i=1;i<4;i++){
            if(pmax>=act.pol[i][7]&&pmax<act.pol[i][8]){
                kn_max=act.pol[i][0] + act.pol[i][1]*pmax + act.pol[i][2]*pow(pmax,2) + act.pol[i][3]*pow(pmax,3) + act.pol[i][4]*pow(pmax,4) + act.pol[i][5]*pow(pmax,5) + act.pol[i][6]*pow(pmax,6);
                ;
            }
        }
        this->p_max = pmax;
    }

    void SetInt(){
        if(i_is&&i_os){
            if(i_es&&(g_do-g_dc)/4>g_sl/2){
                x_inc=(g_sl/2)/sim_i;
            }
            else{
                x_inc=((g_do-g_dc)/4)/sim_i;
            }
        }
        else{
            if(i_es&&(g_do-g_dc)/2>g_sl/2){
                x_inc=(g_sl/2)/sim_i;
            }
            else{
                x_inc=((g_do-g_dc)/2)/sim_i;
            }
        }
    }

    void SimA(){
        db d=g_dc;
        db D=g_do;
        db L=g_l;
        db Ae;//Ends Area
        db br_amax=(i_es*(D*D-d*d)*2*g_n*pi4)+(i_is*d*pi*L)+(i_os*D*pi*L);
        for(int x=0;x<sim_i+1;x++){
            Ae=D*D-d*d;
            sim[x][0]=x*x_inc;//Regression
            sim[x][2]=(i_es*Ae*2*g_n*pi4)+(i_is*d*pi*L)+(i_os*D*pi*L);//Grain burn area
            sim[x][11]=Ae*L*pi4/1000000000; //Grain Volume
            sim[x][12]=c_cs-Ae*pi4;//Combustion Flow Area
            //cout<<sim[x][0]<<"\t"<<d<<"\t"<<D<<"\t"<<L<<"\t"<<sim[x][2]<<"\t"<<sim[x][11]<<"\t"<<sim[x][12]<<endl;
            d+=i_is*x_inc*2;
            D-=i_os*x_inc*2;
            L-=i_es*x_inc*g_n*2;
            if(sim[x][2]>br_amax){
                br_amax=sim[x][2];
            }
        }
        //cout<<kn_max<<endl;
        //cout<<br_amax<<endl;
        n_tdi=sqrt(br_amax/kn_max/pi)*2;
        //cout<<n_tdi<<endl;
        //cout<<endl<<endl<<endl<<endl<<endl<<endl<<endl;
    }

    double GetR(db P){
       return act.GetRP(P);
    }

    void SimB(){
        //Set initial conditions
        for(int i=0;i<11;i++){
            sim[0][i]=0;
        }

        //Full Sim
        db r;//regression rate
        db fv;//free chamber volume
        db mg;//mass generation
        db patm=0.101;//atmospheric pressure (MPa)
        db sden=act.p_rden*1000;//supplemental density (kg/mm^3)
        db na=n_tdi*n_tdi*pi4/1000000; //nozzle area (m^3)
        db msr;//mass storage rate
        db ms;//stored mass
        db ms2=0;//old stored mass
        db denpro=0; //density of products
        db cfa=1;
        r=GetR(sim[0][9]);
        sim[0][1]=n_tdi;//set initial nozzle diameter
        sim[0][2]=br_ai;//set initial Burn Area
        sim[0][4]=sim[0][11]*sden;//set initial Propellant Mass
        sim[0][5]=0;
        sim[0][9]=0.101;//set initial pressure
        sim[0][13]=1;//set initial expansion ratio
        //cout<<na<<endl;
        for(int x=1;x<sim_i+1;x++){
            r=GetR(sim[x-1][9]);//get regression rate
            sim[x][4]=sim[x][11]*sden;//get grain mass
            sim[x][5]=x_inc/r+sim[x-1][5];//get time
            fv=c_v-sim[x][11];//get free volume
            mg=(sim[x-1][4]-sim[x][4])/(sim[x][5]-sim[x-1][5]);//get mass generation
            //cout<<(sim[x-1][4]-sim[x][4])<<"\t"<<(sim[x][5]-sim[x-1][5])<<"\t";
            sim[x][6]=(sim[x-1][9]-patm)*1000000*na/sqrt(act.R*act.ct)*sqrt(act.k)*pow((2/(act.k+1)),((act.k+1)/2/(act.k-1)));//mass flow
            sim[x][7]=sim[x][6]/(sim[x][12]/1000000);//mass flux
            msr=mg-sim[x][6];//buildup rate
            ms=ms2+(sim[x-1][4]-sim[x][4])*msr;//buildup mass
            denpro=ms/fv;
            sim[x][9]=denpro*act.R*act.ct/1000000+patm;
            sim[x][13]=1/(pow(((act.k2ph+1)/2),(1/(act.k2ph-1)))*pow((patm/sim[x][9]),(1/act.k2ph))*sqrt((act.k2ph+1)/(act.k2ph-1)*(1-pow((patm/sim[x][9]),((act.k2ph-1)/act.k2ph)))));
            cfa+=sim[x][13];
            //end
            ms2=ms;
            //cout<<"\t p:"<<sim[x][9]<<"\t Ar:"<<sim[x][13]<<"\t cfa:"<<cfa<<"\t k2ph:"<<act.k2ph<<endl;
            //cout<<"\t r:"<<r<<"\t m:"<<sim[x][4]<<"\t t:"<<sim[x][5]<<"\t mg:"<<mg<<"\t mf:"<<sim[x][6]<<"\t mfx:"<<sim[x][7]<<"\t R:"<<act.R<<"\t ms:"<<ms<<"\t denpro:"<<denpro<<"\t p:"<<sim[x][9]<<"\t Ar:"<<sim[x][13]<<"\t cfa:"<<cfa<<endl;
        }
        int i=sim_i;
        /*while(sim[i][9]>patm+pmax*0.02){

        }*/
        cfa=cfa/sim_i;
        //cout<<cfa;
        db Pe=.101;
        db Ae=na*cfa;

        for(int x=1;x<sim_i+1;x++){
            sim[x][10]=(n_ef*sqrt(2*pow(act.k2ph,2)/(act.k2ph-1)*pow((2/(act.k2ph+1)),((act.k2ph+1)/(act.k2ph-1)))*(1-pow((Pe/sim[x][9]),((act.k2ph-1)/act.k2ph))))+(Pe-patm)/sim[x][9]*(Ae/na))*(sim[x][9]*na*1000000);
            cout<<sim[x][10]<<endl;
        }
    }
};



int main(){
//input
//cout<<"Input Number of Motors: ";
//cin>>n;
/**Read Propellant*/
/*  Input propellant
    Save propellant file
    Read Propellant file
*/
Propellant KNDX;
KNDX.SetProp(1.879, 1.043, 1.131, 42.39, 1710,0.95,0.95);
//Set Burn Rate Data
//(a, n, Min Pressure (MPa), Max Pressure (MPa))
db KNDX_br[7][6]={
    {0,0,0,0,0,0},
    {0,0,0,0,0,0},
    {0,8.875, 0.619, 0.100, 0.779,0},
    {0,7.553,-0.009, 0.779, 2.572,0},
    {0,3.841,	0.688, 2.572, 5.930,0},
    {0,17.20,	-0.148,5.930, 8.502,0},
    {0,4.775, 0.442, 8.502, 11.20,0}
};

KNDX.SetBurn(5,KNDX_br);
//Set KN prediction data
db KNDX_pol[4][10]={
    //(a,b,c,d,e,f,g,pMin,pMax)
    {43.500, 0.24168,50.484,-9.9115,0,0,0,1.03425,3.10275},
    {163.80, 22.224,-0.2240,0,0,0,0,3.10275,6.2055},
    {5095.0, -2460.4,453.48,-35.532,1.0175,0,0,6.2055,8.9635}
};
KNDX.SetPol(KNDX_pol);
/*  Input motor properties or read motor file
*/
Motor m;
m.SetData(32,192,32,12.8,64,3,0,1,1,0.85,KNDX);
/**Burn Sim**/
/*
Layout:
0   -Regression (mm)
1   -Nozzle diameter (mm)
2   -Burn Area (mm^2)
3   Kn----------------------
4   -Propellant Mass (kg)
5   -Time (s)
6   -Mass Flow (kg/s)
7   Mas Flux   (kg/s.m^2)
8   Nozzle Exit Pressure (MPa)
9   Chamber Pressure   (MPa)
10  Thrust  (N)
11  Grain volume    (mm^3)
12  Combustion Flow Area    (mm^2)
13  Optimum Nozzle expansion ratio

*/
db ri;//web regression interval
/**Bdb c_di, db c_l, db g_do, db g_dc, db g_sl, db g_n, bool i_os, bool i_is, bool i_es,Propellant act**/

db Patm=0.101; //Ambient pressure (MPa)
ri = (m.g_do-m.g_dc)/(2*sim_i);

//Calculating x_inc
m.SetInt();
//Calculating Kn max
m.KnMax(6.895);

//Area Sim
m.SimA();

m.SimB();

};

