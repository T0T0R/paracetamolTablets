#include <iostream>
#include <cstdlib>
#include <map>
#include <vector>
#include <math.h>
#include <fstream>
#include <string>


/***
Calcul de pH par recherche de la racine d'un polynome.
Six especes presentes, chacune avec une concentration donnee.
Le programme calcule le pH pour chaque goutte de titrant ajoutee.
    Une fois cette burette 'videe', on recommence en changeant une des concentrations
des cinq autres especes, jusqu'a obtenir une courbe de pH se rapprochant suffisemment
de celle obtenue lors du titrage. On note alors les masses initiales de reactifs qui ont
engendre ce resultat.
    Ces concentrations sont calculees avec des masses initiales allant
de (m_min) mg a (m_max) mg par pas de (m_gap) mg.
***/


/***
Les seules especes dont on connait deja la concentration correcte sont:
    - Titrant : c0
    - Paracetamol (ParacH): c5
autres especes :
    - Acide tartrique (AH_3) : c1
    - CO_3 2- (2 Na+ probablement) : c2
    - HCO_3 - (Na+ probablement) : c3
    - PhCOO - (Na+ probablement) : c4
***/



/*** Functions designed to return concentrations of species in solution ***/
double getAH3(double const& h, double const& c1, std::vector<double> const& Ka){
	//Get AH3 concentration
	return c1/(1+Ka[1]/h*(1+Ka[2]/h*(1+Ka[3])));
}
double getAH2(double const& h, double const& c1, std::vector<double> const& Ka) {
	return c1/(1+h/Ka[1]+Ka[2]/h*(1+Ka[3]));
}
double getAH(double const& h, double const& c1, std::vector<double> const& Ka) {
	return c1/(1+h/Ka[2]*(h/Ka[1]+1)+Ka[3]/h);
}
double getA(double const& h, double const& c1, std::vector<double> const& Ka) {
	return c1/(1+h/Ka[3]*(h/Ka[2]*(h/Ka[1]+1)+1));
}
std::vector<double> getCO2(double const& h, double const& c2, double const& c3, std::vector<double> const& Ka, double const& pCO2air) {
	/* Returns a couple containing:
	table[0] : the concentration of CO2 diluted in solution
	table[1] : the 'concentration' of CO2 that cannot be disolved and therefore is being released
	*/
	double conc= (c2+c3)/(1+Ka[4]/h*(1+Ka[5]/h));
	std::vector<double> result (2);

	double k_Hcp (3.4e-2);	//Constante de Henry pour le CO2 a 25 °C
	double diff = conc-k_Hcp/101.3e3*pCO2air;
	if (diff>0){
		result[0] = k_Hcp/101.3e3*pCO2air;
		result[1] = diff;
	}else{
		result[0] = conc;
		result[1]=0;
	}

	return result;
}
double getHCO3(double const& h, double const& c2, double const& c3, std::vector<double> const& Ka, double const& pCO2air) {
	return getCO2(h, c2, c3, Ka, pCO2air)[0]*Ka[4]/h;
}
double getCO3(double const& h, double const& c2, double const& c3, std::vector<double> const& Ka, double const& pCO2air) {
	return getHCO3(h, c2, c3, Ka, pCO2air)*Ka[5]/h;
}
double getPhCOOH(double const& h, double const& c4, std::vector<double> const& Ka) {
	return c4/(1+Ka[6]/h);
}
double getPhCOO(double const& h, double const& c4, std::vector<double> const& Ka) {
	return c4/(1+h/Ka[6]);
}
double getParacH(double const& h, double const& c5, std::vector<double> const& Ka) {
	return c5/(1+Ka[7]/h);
}
double getParac(double const& h, double const& c5, std::vector<double> const& Ka) {
	return c5/(1+h/Ka[7]);
}



std::vector<std::vector<double>> readDatas(std::string filename){
	/** \brief Stocke les donnees experimentales dans un tableau de ph=f(VolumeAjoute)
	*
	* \param filename {Nom du fichier de donnees}
	* \return {Tableau de couples (V;pH)}
	*
	*/

	std::vector<std::vector<double>> datas;
	std::string delimiter = ";";
	std::ifstream myFile(filename);
	std::string line;

	if(myFile.is_open()){
		while (getline(myFile, line)){
			std::string volume = line.substr(0, line.find(delimiter));
			std::string pH = line.substr(line.find(delimiter)+1, line.find("\n"));
			datas.push_back(std::vector<double>{  std::atof(volume.c_str())/1000, std::atof(pH.c_str())  });
		}
		myFile.close();
		
	}else{
		std::cout<< "Unable to open file" <<std::endl;
	}

	return datas;
}



double polynomeB(double const& h, double const& c0, double const& c1, double const& c2, double const& c3, double const& c4, double const& c5,
				 std::vector<double> const& Ka, double const& pCO2air) {
	/** \brief Polynome dont la racine est la concentration en ions H+ de la solution pour le titrage par une base forte
	*
	* \param h {Indeterminee du polynome}
	* \param c {Concentrations du milieu}
	* \param Ka {Constantes de dissociation}
	* \return polynomeB(h)
	*
	*/

	return Ka[0]/h + getAH2(h, c1, Ka) + 2*getAH(h, c1, Ka) + 3*getA(h, c1, Ka) - getCO2(h, c2, c3, Ka, pCO2air)[0] - getCO2(h, c2, c3, Ka, pCO2air)[1] + getCO3(h, c2, c3, Ka, pCO2air) - c2
		- getPhCOOH(h, c4, Ka) + getParac(h, c5, Ka) - h - c0;
	//return Ka[0]/h  +  c1/(1+h/Ka[1]+Ka[2]/h*(1+Ka[3]/h))*(1+Ka[2]/h*(2+3*Ka[3]/h))  +  (c2+c3)/(1+h/Ka[5]*(1+h/Ka[4]))*(1-h*h/(Ka[4]*Ka[5]))  -  c2  +  c5/(1+h/Ka[7])  -  c4/(1+Ka[6]/h)  -  c0  -  h;
}



double pH(double const& c0, double const& c1, double const& c2, double const& c3, double const& c4, double const& c5, std::vector<double> const& Kas, double const& pCO2air){
	/** \brief Calcule le pH en cherchant la racine d'un polynome par dichotomie.
	*
	* \param c {Concentrations du milieu}
	* \param Kas {Constantes de dissociation}
	* \return pH du melange.
	*
	*/

	double a = pow(10, -14);	//pH minimal
	double b = pow(10, 0);	//pH maximal
	double eps = pow(10, -16);	//erreur accordee

	double x(0.0);

	if (polynomeB(a, c0, c1, c2, c3, c4, c5, Kas, pCO2air)*polynomeB(b, c0, c1, c2, c3, c4, c5, Kas, pCO2air)>0){
			return -1;
	}

	while (b-a > eps){
		x = (b+a)/2;
		if (polynomeB(x, c0, c1, c2, c3, c4, c5, Kas, pCO2air)*polynomeB(a, c0, c1, c2, c3, c4, c5, Kas, pCO2air)>0){
			a = x;
		}else{
			b = x;
		}
	}
	return -log10((a+b)/2);
}



double compares(double const& vaj, double const& pH, std::vector<std::vector<double>> const& table){
    /** \brief Compare deux tableaux de pH = f(V ajoute) avec une tolerance de eps.
     *
     * \param vaj {Volume ajoute actuel}
     * \param pH {pH a vaj}
     * \param table {Tableau de pH = f(V ajoute)}
     * \param eps {Seuil de pH}
     * \return true si les tableaux sont analogues, false sinon.
     *
     */

	double epsVolume (0.00025);	//Seuil de detection pour le volume

	for (std::vector<double> couple: table){	//Pour chaque couple VolumeAjoute/pH...
		if (abs(couple[0]-vaj)<epsVolume){	//... on cherche si c'est celui qui correspond au bon volume vaj ...
			return couple[1]-pH;	//... et dans ce cas, on retourne la difference.
		}else{	//Si le couple n'est pas celui avec le bon volume ajoute...
			continue;	//... on passe au couple suivant.
		}
	}

     return 0;	//Dans le cas ou le point recherche correspond a un volume qui n'est pas present dans les donnees experimentales, on l'accepte par defaut
}



double getExhaustedVol(double const& h, double const& c2, double const& c3, std::vector<double> const& Ka,
					double const& pAir, double const& xCO2, double const& volume, double const& T){
	//Retourne le volume de gaz CO2 degage
	return 8.134*T/pAir * getCO2(h, c2, c3, Ka, pAir*xCO2)[1]/volume;	//	V=RT/P*n
}



std::vector<double> getConc(double const& h, double const& c1, double const& c2, double const& c3, double const& c4, double const& c5, 
							std::vector<double> const& Ka, double const& pCO2air) {
	/** \brief Retourne les concentrations des diverses especes presentes en solution a l'equilibre
	*
	* \param h {Concentration en H+}
	* \param c {Concentrations initiales}
	* \param Ka {Constantes de dissociation}
	* \param pCO2air {Pression partielle du CO2 dans l'air}
	* \return concentrations
	*
	*/

	std::vector<double> conc;
	conc.push_back(c1/(1+Ka[1]/h*(1+Ka[2]/h*(1+Ka[3]))));	//AH3
	conc.push_back(c1/(1+h/Ka[1]+Ka[2]/h*(1+Ka[3])));	//AH2 -
	conc.push_back(c1/(1+h/Ka[2]*(h/Ka[1]+1)+Ka[3]/h));	//AH 2-
	conc.push_back(c1/(1+h/Ka[3]*(h/Ka[2]*(h/Ka[1]+1)+1)));	//A 3-

	conc.push_back(getCO2(h,c2,c3,Ka,pCO2air)[0]);	//CO2 dissous
	conc.push_back((c2+c3)/(1+h/Ka[4] +Ka[5]/h));	//HCO3 -
	conc.push_back((c2+c3)/(1+h/Ka[5]*(1+h/Ka[4])));	//CO3 2-

	conc.push_back(c4/(1+Ka[6]/h));	//PhCOOH
	conc.push_back(c4/(1+h/Ka[6]));	//PhCOO -

	conc.push_back(c5/(1+Ka[7]/h));	//ParacH
	conc.push_back(c5/(1+h/Ka[7]));	//Parac -

	return conc;
}



void titrage(double const& c_titrant, double const& v_burette, double const& c1, double const& c2, double const& c3, double const& c4, double const& c5,
			 double const& vi, std::vector<double> const& Kas, double const& pAir, double const& xCO2, double const& T, int const& nbCurve) {
	/** \brief Realise un titrage d'une solution
	*
	* \param c_titrant {Concentration de l'espece titrante}
	* \param v_burette {Volume total de titrant a verser}
	* \param c {Concentrations initiales}
	* \param Ka {Constantes de dissociation}
	* \param vi {Volume de solution a titrer (volume initial)}
	* \param pAir {Pression de l'air}
	* \param xCO2 {Fraction molaire de CO2 dans l'air}
	* \param T {Temperature}
	* \return concentrations
	*
	*/



	std::ofstream myFile("testpH_"+std::to_string(nbCurve)+".csv");

	myFile<<"Volume ajoute"<<";"<<"pH"
		<<";"<<"[AH_3]"<<";"<<"[AH_2 -]"<<";"<<"[AH 2-]"<<";"<<"[A 3-]"
		<<";"<<"[CO2]"<<";"<<"[HCO3 -]"<<";"<<"[CO3 2-]"
		<<";"<<"[PhCOOH]"<<";"<<"[PhCOO -]"
		<<";"<<"[ParacH]"<<";"<<"[Parac -]"
		<<";"<<"Vol degaze"<<std::endl;

	for (int vg(5); ((double)vg)/100000<=v_burette; vg += 5) {	//vg : centieme de mL
		double vaj = ((double)vg)/100000;
		double vtot = vaj+vi;
		double c0 = c_titrant*vaj/(vtot);	//Concentration de titrant dans le becher

		double pHnow = pH(c0, c1*vi/vtot, c2*vi/vtot, c3*vi/vtot, c4*vi/vtot, c5*vi/vtot, Kas, pAir*xCO2);
		std::vector<double> conc = getConc(pow(10, -pHnow), c1*vi/vtot, c2*vi/vtot, c3*vi/vtot, c4*vi/vtot, c5*vi/vtot, Kas, pAir*xCO2);

		myFile<<vaj*1000<<";"<<pHnow
			<<";"<<conc[0]<<";"<<conc[1]<<";"<<conc[2]<<";"<<conc[3]
			<<";"<<conc[4]<<";"<<conc[5]<<";"<<conc[6]
			<<";"<<conc[7]<<";"<<conc[8]
			<<";"<<conc[9]<<";"<<conc[10]
			<<";"<<getExhaustedVol(pow(10, -pHnow),c2,c3,Kas,pAir,xCO2,vtot,T)<<std::endl;	//Concentrations des especes dans le becher
	}
	myFile.close();
}




int main()
{
	//Constantes de dissociation
	double Ke = pow(10, -14);  //H2O <=> HO - + H+
	double Ka1 = pow(10, -3.13);  //AH_3 <=> AH_2 - + H+
	double Ka2 = pow(10, -4.76);  //AH_2 - <=> AH 2- + H+
	double Ka3 = pow(10, -6.39);  //AH 2- <=> A 3- + H+

	double Ka4 = pow(10, -6.35);  //CO_2 <=> HCO_3 - + H+
	double Ka5 = pow(10, -10.33);  //HCO_3 - <=> CO_3 2- + H+


	double Ka6 = pow(10, -4.20);  //PhCOOH <=> PhCOO - + H+
	double Ka7 = pow(10, -9.38);  //ParacH <=> Parac - + H+

	std::vector<double> Kas {Ke, Ka1, Ka2, Ka3, Ka4, Ka5, Ka6, Ka7};

	//Masses molaires (g/mol)
	double M_AH3 (192.12);	//Acide citrique
	double M_NaHCO3 (84.01);	//Bicarbonate de Sodium
	double M_Na2CO3 (105.99);	//Carbonate de Sodium
	double M_PhCOONa (144.10);	//Benzoate de Sodium
	double M_ParacH (151.16);	//Parcetamol


	//Donnees du titrage
	std::vector<std::vector<double>> expTitrage = readDatas("dataFile.csv");
	double vi (0.02);	//Volume de solution initiale dans le becher (L)
	double v_burette (0.04);	//Volume de la burette (L)

	double c_titrant (0.2);	//Concentration de titrant dans la burette (mol/L)
	double m_ParacH (500.0/1000);	//Masse de Paracetamol (g)


	double c5 = m_ParacH/(vi*M_ParacH);	//Concentration de Paracetamol dans le becher (mol/L)

	double pressure (101.3e3);	//Pression de l'air (Pa)
	double xCO2air (0.04);	//Fraction molaire de CO2 dans l'air
	double pCO2air = pressure*xCO2air;	//Pression partielle de CO2 dans le milieu exterieur (Pa)
	double temperature (293);	//Temperature du milieu (K)


	//Donnees de la simulation
	double eps (1.0);	//Ecart de pH maximal admissible
	double m_min (100.0);	//Masse minimale (mg)
	double m_max (5000.0);	//Masse maximale (mg)
	double m_gap (50.0);	//Pas entre deux masses (mg)


	//Initialisations de variables
	double pHnow (-1);
	bool stillCorrect (false);

	double prevSum (1E8);	//Variance initiale fixee volontairement tres grande
	double sum (0);
	int i (0);
	double gap (0);	//Ecart mesure entre simulation/experience
	int nbCurve (0);	//Numerote les courbes de sortie

	double _temp (0.0);


	/*** Entree utilisateur ***/
	std::cout<<"Volume becher (mL): ";
	{std::cin>>_temp;	vi = _temp/1000;}
	std::cout<<"Volume burette (mL): ";
	{std::cin>>_temp;	v_burette = _temp/1000;}
	std::cout<<"Concentration de l'espece titrante (mol/L): ";
	std::cin>> c_titrant;
	std::cout<<"Ecart de pH maximal: ";
	std::cin>> eps;
	std::cout<<"Masse minimale (mg): ";
	std::cin>>m_min;
	std::cout<<"Masse maximale (mg): ";
	std::cin>>m_max;
	std::cout<<"Pas entre deux masses successives (mg): ";
	std::cin>>m_gap;




	for (double m_PhCOONa = m_min; m_PhCOONa<= m_max; m_PhCOONa += m_gap) {	//Masse Benzoate de sodium (en mg)
		double c4 = m_PhCOONa/(1000*vi*M_PhCOONa);
		std::cout<<((double)m_PhCOONa)/((double)m_max)*100<<"\t%"<<std::endl;	//Affichage progression en %

		for (double m_NaHCO3 = m_min; m_NaHCO3<= m_max; m_NaHCO3 += m_gap) {	//Masse Bicarbonate de sodium (en mg)
			double c3 = m_NaHCO3/(1000*vi*M_NaHCO3);

			for (double m_Na2CO3 = m_min; m_Na2CO3<= m_max; m_Na2CO3 += m_gap) {	//Masse Carbonate de sodium (en mg)
				double c2 = m_Na2CO3/(1000*vi*M_Na2CO3);

				for (double m_AH3 = m_min; m_AH3<= m_max; m_AH3 += m_gap) {   //Masse Acide citrique (en mg)
					double c1 = m_AH3/(1000*vi*M_AH3);	//Concentration Acide citrique


					stillCorrect = true;	//Si la simulation d'un titrage est coherente avec l'experience jusque-la
					sum = 0;	//Somme des ecarts a l'experience en vue de calculer une variance
					i = 1;	//Nombre de points du titrage pour calculer la variance

					for (int vg(5); ((double)vg)/100000<=v_burette; vg+=50){	//vg : centieme de mL
						i++;
						double vaj = ((double)vg)/100000;	//Volume ajoute calcule a partir de vg (L)
						double vtot = vaj+vi;
						double c0 = c_titrant*vaj/(vtot);	//Concentration de titrant dans le becher


						pHnow = pH(c0, c1*vi/vtot, c2*vi/vtot, c3*vi/vtot, c4*vi/vtot, c5*vi/vtot, Kas, pCO2air);	//Concentrations dans le becher

						gap = compares(vaj, pHnow, expTitrage);

						if(abs(gap)>eps){
							stillCorrect = false;
							break;	//Sort du titrage si le pH ne correspond plus a partir d'un moment
						}else{
							sum += pow(gap, 2);
						}
					}


					sum /= i;

					if (stillCorrect && sum < prevSum){	//Si tous les points de la simulation sont acceptables et que la variance est plus faible que la precedente
						std::cout<<"s = "<< sqrt(sum)<<std::endl;	//Affichage ecart-type avec l'experience
						std::cout<<"m_AH3 = "<< m_AH3 << "mg."<<std::endl;
						std::cout<<"m_Na2CO3 = "<< m_Na2CO3 << "mg."<<std::endl;
						std::cout<<"m_NaHCO3 = "<< m_NaHCO3 << "mg."<<std::endl;
						std::cout<<"m_NaPhCOO = "<< m_PhCOONa << "mg.\n"<<std::endl;

						titrage(c_titrant, v_burette, c1, c2, c3, c4, c5, vi, Kas, pressure, xCO2air, temperature, nbCurve);
						prevSum = sum;	//Nouvelle variance a 'battre'
						nbCurve ++;
					}
				}
			}
		}
	}
	return EXIT_SUCCESS;
}