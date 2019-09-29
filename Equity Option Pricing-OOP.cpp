//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Project Name: European Equity Option (non dividend paying stock) Pricing and Risk Generation
//Version: 1.0
//Project Description: Calculates Equity Option Prices, Option Greeks - Delta, Vega, Gamma, Theta, Rho, does Vol Smile Fitting based on
//					   user inputted Volatility/ Strike grid
//Date: 28-Sep-2019
//Owner: Vivek Kumar Rauniyar
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Including C++ libraries
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <iomanip>
#include <bits/stdc++.h>


class Option{

    public:
    
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// This calculates d_j, for j in {1,2}. This term appears in the closed form solution for the European call or put price
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    
    double d_j(const int& j, const double& S, const double& K, const double& r, const double& v, const double& t) {
   		return (log(S/K) + (r + (pow(-1,j-1))*0.5*v*v)*t)/(v*(pow(t,0.5)));
	}	

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// An approximation to the cumulative distribution function for the standard normal distribution
	//This is a recursive function. I have used this piece of code from quantlib
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double norm_cdf(const double& x) {
    double k = 1.0/(1.0 + 0.2316419*x);
    double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));

	    if (x >= 0.0) {
	        return (1.0 - (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x) * k_sum);
	    } else {
	        return 1.0 - norm_cdf(-x);
	    }
    
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Greek: Gamma Calculator
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	double gamma (const double& S,const double& K,const double& r,const double& v,const double& t)
	{
		return ((pow(2.71828,(-pow(d_j(1,S,K,r,v,t),2)/2)))/pow(2*M_PI,0.5))/(S*v*pow(t,0.5));
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Greek: Vega Calculator
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	double vega (const double& S,const double& K,const double& r,const double& v,const double& t)
	{
		return S*((pow(2.71828,(-pow(d_j(1,S,K,r,v,t),2)/2)))/pow(2*M_PI,0.5))*pow(t,0.5)/100;
	}
	



	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Volatility Smile Fitting using user inputted Volatility Strike Grid
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double * Vol_Smile_Fitting ()
	{
	
				double Vol[20][20],SABR_param[10];
				double ATM_Vol,Volatility;
				int * loc;
				double ATM_Strike,Strike;
				double Strike_Increment,Strike_Range,Strike_Range_Low,Strike_Range_High,Strike_count;
				double A,B,C;
				double A_high=10,B_high=10,C_high=10;
				double A_low=-10.0,B_low=-10.0,C_low=-10.0;
				double A_Fitted,B_Fitted,C_Fitted;
				int ATM;
				
										
				std::cout << "Enter ATM Strike\n";
				std::cin >> ATM_Strike;
				
				//std::cout << "Enter Strike Increment for the Smile Grid. Please enter multiples of 5.\n";
				//std::cin >> Strike_Increment;
				
				Strike_Increment=5;
				
				std::cout << "Enter Lowest Strike for the Strike Range for the Smile Grid. Please enter multiples of 5.\n";
				std::cin >> Strike_Range_Low;
				
				std::cout << "Enter Highest Strike for the Strike Range for the Smile Grid. Please enter multiples of 5.\n";
				std::cin >> Strike_Range_High;
				
				
				Strike_count = (Strike_Range_High - Strike_Range_Low)/Strike_Increment+1;
				
				Strike = Strike_Range_Low;
				
				//std::cout << Strike_count << "      " << Strike << std::endl;
				
				for(int i=0;i<Strike_count;i++){
				
					Vol[0][i]=Strike/ATM_Strike;		
					std::cout << "\n";
					
					std::cout << "Enter Implied Vol in % for the Strike  : " << Strike << std::endl;
					std::cin >> Volatility;
					Vol[1][i]=Volatility/100;
					
					
						if(Strike==ATM_Strike){
							ATM=i;
							ATM_Vol=Volatility/100;
							}
								
					Strike = Strike + Strike_Increment;
					
				}
				
				
				
		/////////////////////////////////////////////////////////////////////////////////////////////////		
		//  Vol Smile Fitting using quadratic equation		
		////////////////////////////////////////////////////////////////////////////////////////////////
				
				double Min_Std_err=1, Std_err=0,Err=0, ATM_Diff=0.2,ATM_Lowest_Diff=0.5,A_ATM_lowest=0,B_ATM_lowest=0,C_ATM_lowest=0;
				double Min_Std_err_ATM_Fixed=1,A_Fitted_ATM_Fixed,B_Fitted_ATM_Fixed,C_Fitted_ATM_Fixed;
				
				for(A=A_low;A<=A_high;A+=0.01){
					for(B=B_low;B<=B_high;B+=0.1){
						for(C=C_low;C<=C_high;C+=0.1){
							
							Std_err=0;
							Err=0;
							
							for(int i=0;i<Strike_count;i++){
										
										if(Vol[0][i]!=1){
										Err=pow((Vol[1][i]-(A*pow(Vol[0][i],2.0)+B*Vol[0][i]+C)),2.0)*(1/sqrt(abs(ATM_Strike*Vol[0][i]-ATM_Strike)));
										}
										if(Vol[0][i]==1){
										Err=pow((Vol[1][i]-(A*pow(Vol[0][i],2.0)+B*Vol[0][i]+C)),2.0);
										}
										Std_err=Std_err+Err;
							}
							
							
							Std_err= pow(Std_err,0.5)/Strike_count;
							
													
							for(int i=0;i<Strike_count;i++){
								Vol[3][i]=A*pow(Vol[0][i],2.0)+B*Vol[0][i]+C;
								}
							
							if(Min_Std_err>Std_err){
								Min_Std_err=Std_err;
								A_Fitted = A;
								B_Fitted = B;
								C_Fitted = C;
							}
						}
					}
				}		
			
				for(int i=0;i<Strike_count;i++){
					Vol[2][i]=A_Fitted*pow(Vol[0][i],2.0)+B_Fitted*Vol[0][i]+C_Fitted;		
				}
				
	
				
				SABR_param[0]=A_Fitted;
				SABR_param[1]=B_Fitted;
				SABR_param[2]=C_Fitted;
				                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
				
				std::cout << "\n\n********************************************************\n";
				std::cout << "The best fit estimates for the parameters are the below\n";
				std::cout << "********************************************************\n";
				std::cout << "Convexity (a.k.a VoV)                 :" << A_Fitted << std::endl;
				std::cout << "Skew (a.k.a. Correlation)             :" << B_Fitted << std::endl;
				std::cout << "Constant                              :" << C_Fitted << std::endl;
				std::cout << "Minimum Standard error                :" << Min_Std_err << std::endl;
		
				std::cout << "\n*******************************************************************************************\n";
				std::cout << "\nRelative Strike        Implied Vol             Fitted_Vol            " << std::endl;
				std::cout << "**********************************************************************************************\n";
				
				for (int i=0;i<Strike_count;i++){
							std::cout << "     "<< std::setprecision(6) << Vol[0][i] << "                "<< std::setprecision(6)<< Vol[1][i] << "                      "<< std::setprecision(6) << Vol[2][i] << std::endl;
				}
				
				
					
				return SABR_param;
				
	}


      
      
}; 

class CallOption : public Option{
	
	public:
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//This calculates Call option prices based on user input
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	double Call_Opt_Pricer (const double& S,const double& K,const double& r,const double& v,const double& t)
	{
		double Call_Opt_Premium;
		
		Call_Opt_Premium = S*norm_cdf(d_j(1,S,K,r,v,t)) - K*exp(-r*t)*norm_cdf(d_j(2,S,K,r,v,t));
		
		return Call_Opt_Premium;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Greek: Delta Calculator
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	double delta (const double& S,const double& K,const double& r,const double& v,const double& t)
	{
		return norm_cdf(d_j(1,S,K,r,v,t));
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Greek: rho Calculator
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	double rho (const double& S,const double& K,const double& r,const double& v,const double& t)
	{
		return K*t*pow(2.71828,-r*t)*norm_cdf(d_j(2,S,K,r,v,t))/100;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Greek: theta Calculator
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	double theta(const double& S,const double& K,const double& r,const double& v,const double& t)
	{
		double theta_call;
		
		theta_call= (-(S*v*((pow(2.71828,(-pow(d_j(1,S,K,r,v,t),2)/2)))/pow(2*M_PI,0.5))/(2*pow(t,0.5)))-(r*K*pow(2.71828,-r*t)*norm_cdf(d_j(2,S,K,r,v,t))));
		
		return theta_call;
	}
		
};

class PutOption: public Option{
	
	public:

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//This calculates Put option prices based on user input
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	double Put_Opt_Pricer (const double& S,const double& K,const double& r,const double& v,const double& t)
    {	
		double Put_Opt_Premium;
		
		Put_Opt_Premium = K*exp(-r*t)*norm_cdf(-d_j(2,S,K,r,v,t)) - S*norm_cdf(-d_j(1,S,K,r,v,t));
		
		return Put_Opt_Premium;
    }

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Greek: delta Calculator
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	double delta (const double& S,const double& K,const double& r,const double& v,const double& t)
	{
		return norm_cdf(d_j(1,S,K,r,v,t))-1;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Greek: rho Calculator
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
		double rho (const double& S,const double& K,const double& r,const double& v,const double& t)
	{
		return -K*t*pow(2.71828,-r*t)*norm_cdf(-d_j(2,S,K,r,v,t))/100;
	}
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Greek: theta Calculator
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	double theta(const double& S,const double& K,const double& r,const double& v,const double& t)
	{
		double theta_put;
		
		theta_put= (-(S*v*((pow(2.71828,(-pow(d_j(1,S,K,r,v,t),2)/2)))/pow(2*M_PI,0.5))/(2*pow(t,0.5)))+(r*K*pow(2.71828,-r*t)*norm_cdf(-d_j(2,S,K,r,v,t))));
		
		return theta_put;
	}
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//MAIN Function: Takes in user inputs for Option Prices, calculates Option Greeks and prints various calculated outputs
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv)
{
	double S,K,r,v,t; //Declares Spot Price, Strike Price, Risk Free rate, Volatility and Time to expiry
	char Opt_type;	  //Can be either 'C' ( for Call) or 'P' (for Put)
	double *pVol;	  
	int Vol_Smile_flg=0;
	double Opt_Premium, Cal_Opt_Premium, Cal_Opt_Premium_ATM_Fixed,Opt_Premium_tweak, Cal_Vol,Cal_Vol_ATM_Fixed;
	double S_rel_positivetweak, S_rel_negativetweak,S_abs_positivetweak,S_abs_negativetweak;
	double delta,gamma,vega,theta,rho;
	double delta1,gamma1,gamma2,vega1,theta1,rho1;
	char Job_select;
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Taking User inputs for option pricing
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	std::cout << "What you would like to generate: Enter P for Option Pricing / Enter V for Vol Smile Fitting \n";
	std::cin >> Job_select;
	
		if(Job_select != 'P'){
			if(Job_select!='V'){
				std::cout << "Incorrect entry.......exiting \n";
				return 0;
			}
		}
		
		if(Job_select == 'V')
		{
			Option Obj;
			pVol=Obj.Vol_Smile_Fitting();
			Vol_Smile_flg=1;
			
		}
	
		
		if(Job_select == 'P')
		{
			std::cout<<"\n******************************************\n";
			std::cout<<"\n Enter Option Pricing Parameters\n";
			std::cout<<"\n******************************************\n\n";
			
			std::cout << "\nEnter Spot Price: \n";
			std::cin >> S;
			
			std::cout << "\nEnter Strike of Option: \n";
			std::cin >> K;
			
			std::cout << "\nEnter Risk free Interest Rate (in fraction): \n";
			std::cin >> r;
			
			std::cout << "\nEnter Volatility of Underlying (in fraction): \n";
			std::cin >> v;
			
			std::cout << "\nEnter Option Expiry (in years): \n";
			std::cin >> t;
			
			std::cout << "\nWhich Option do you want to price: Enter C for Call/ P for Put \n";
			std::cin >> Opt_type;
			
				if(Opt_type != 'C'){
					if(Opt_type!='P'){
					std::cout << "Incorrect Option type entry.......exiting \n";
					return 0;
					}
				}	
				
				
					if(Opt_type=='P'){
						
					    PutOption Put_O;
						Opt_Premium = Put_O.Put_Opt_Pricer(S,K,r,v,t);
						
						///////////////////////////////////////////////////////////////////////////////////////
						//Calculate Option Greeks
						//////////////////////////////////////////////////////////////////////////////////////
					 
						delta1=Put_O.delta(S,K,r,v,t);
					 	gamma1=Put_O.gamma(S,K,r,v,t);
						rho1=Put_O.rho(S,K,r,v,t);
						vega1=Put_O.vega(S,K,r,v,t);
						theta1=Put_O.theta(S,K,r,v,t);
						
						
									
					    					    
					}else if(Opt_type=='C'){ 
						
						CallOption Call_O;
						Opt_Premium = Call_O.Call_Opt_Pricer(S,K,r,v,t);
						
						///////////////////////////////////////////////////////////////////////////////////////
						//Calculate Option Greeks
						//////////////////////////////////////////////////////////////////////////////////////
											
						delta1=Call_O.delta(S,K,r,v,t);
					 	gamma1=Call_O.gamma(S,K,r,v,t);
						rho1=Call_O.rho(S,K,r,v,t);
						vega1=Call_O.vega(S,K,r,v,t);
						theta1=Call_O.theta(S,K,r,v,t);
											
					}else{
					
						std::cout << "Invalid Option Type \n";
						std::cout << "Closing Program \n";
						return 0; 
					}
		 
			
			system("cls");
			
			std::cout << "\n****************************************************" <<std::endl;
			std::cout << "Underlying Price:               " << S << std::endl;
		    std::cout << "Strike Price:                   " << K << std::endl;
		    std::cout << "Risk-Free Rate:                 " << r << std::endl;
		    std::cout << "Volatility:                     " << v << std::endl;
		    std::cout << "Maturity:                       " << t << std::endl; 
			std::cout << "Option Type:                    " << Opt_type << std::endl;
			std::cout << "Option Premium:                $" << Opt_Premium << std::endl;
			
								
			std::cout << "\n\nOption Greeks are as below" << std::endl;
			std::cout << "\n******************************" << std::endl;
			std::cout << "Delta:              " << delta1 << std::endl;
		    std::cout << "Gamma:              " << gamma1 << std::endl;
		    std::cout << "Vega:               " << vega1 << std::endl;
		    std::cout << "Rho:                " << rho1 << std::endl;
		    std::cout << "Theta:              " << theta1/365 << std::endl;
		}	
	return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//End of Equity Option Pricing and Risk Management module
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
