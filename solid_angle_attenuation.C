{

// dintr-un punct se emit fotoni catre un patrat.
// coeficientul de atenuare al fotonilor la o anumita energie este mu
// probabilitatea a un foton sa NU interactioneze pe distanta x este e^(-mu*x)
// rezulta ca simularea parcursului liber mediu se face cu relatia l = - 1/(mu * ln r)


# define M_PI           3.14159265358979323846  /* pi */

double d = 10.0;     // distanta pana la patrat
double a = 100.0;     // latura patratului
double costeta, teta;
double fi;
double L;           // lungimea pana la plan
double u, v;        // cosinusurile directoare

double N = 1000000.0; // nr. total de aruncari 
double N_in = 0;    // punctele care pica in patrat

double x,y;
double omega;       // unghiul solid

double r1, r2, r3;
//double solid_angle = 4*asin( a*a / sqrt((a*a+4*d*d)*(a*a+4*d*d)) );
double solid_angle = 4*atan (a*a / (2*d * sqrt(4*d*d + a*a + a*a)));

double drum_lib_mediu;  // l = - 1/(mu * ln r)
double mu = 10.0;

TRandom3 *r = new TRandom3(0);

   TFile *f    = new TFile("attenuation.root","RECREATE");
   TTree *tree = new TTree("Tpart","particle tree");
   tree->Branch("x",&x,"x/D");
   tree->Branch("y",&y,"y/D");

for (int i=0; i<N; i++)
{
   r1 = r.Uniform();
   r2 = r.Uniform();
   r3 = r.Uniform();

    costeta = 2*r1-1; 
   //if (costeta <= (d / sqrt(d*d + a*a/2))) continue;

    fi = 2*M_PI*r2;
    teta = acos(costeta);
    
    u = sin(teta)*cos(fi);
    v = sin(teta)*sin(fi);

    L = d/costeta;
 
    drum_lib_mediu = - TMath::Log(r3) / mu;
    if (drum_lib_mediu < L) continue;

    x = L*u;
    if (TMath::Abs(x) > a/2) continue;
   
    y = L*v;
    if (TMath::Abs(y) > a/2) continue;

    N_in +=1; 

    tree->Fill();

}

f->Write();  

tree->Draw("x:y","","colz");

omega = 2*M_PI * N_in / N;

cout << omega << " " << solid_angle << endl;
}
