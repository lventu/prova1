// prova1.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include "box_muller.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <complex>
#include <cfloat>
#include <stdlib.h>
#include <stdio.h>
#include <limits>

using namespace std;

//costanti

#define M_idx 4
#define m_bit 2
#define DATA_LENGTH 1024*1024
#define SYM_LENGTH DATA_LENGTH/m_bit
#define Rb 10			// bit rate
#define Rs Rb/m_bit		// symbol rate
int EbNo_dB= -5;		// SNR
//int N0= 30;				// noise power


//variabili globali
int data_input[DATA_LENGTH];
int symbol_input[SYM_LENGTH];
float K = (float)sqrt(2.0); // è il nostro sqrt(E0)
//numeri complessi uscita
complex<float> out[SYM_LENGTH];
complex<float> out_noise[SYM_LENGTH];
complex<float> ricevuti[SYM_LENGTH];
//costanti
const int mapping[M_idx] = {-3,-1,3,1};
const float psk_I[M_idx] = {(float)cos(-3*PI/M_idx),(float)cos(-1*PI/M_idx),(float)cos(3*PI/M_idx),(float)cos(1*PI/M_idx)};
const float psk_Q[M_idx] = {(float)sin(-3*PI/M_idx),(float)sin(-1*PI/M_idx),(float)sin(3*PI/M_idx),(float)sin(1*PI/M_idx)};

float random(){
	srand(rand()*cos((double)rand())*(int)time(NULL));
	return (float)rand()/RAND_MAX;
}
//BOX MULLER ALGORITHM per ricavare una variabile con pdt gaussiana da 2 va uniformi
float awgn_noise(){
	static bool b_cached = false;
	float randNum,a,b,tmp,res;
	float sigmaI = (float)sqrt(2*2*(float)pow(10.0,(EbNo_dB/10)));
	do{
		randNum = random();
		a = randNum*2 - 1;
		randNum = random();
		b = randNum*2 - 1;
		randNum = a*a + b*b;
	}while(randNum>=1.0 || randNum==0);
	tmp = sqrt((-2*log(randNum))/randNum);
	if (b_cached){
		res = b*tmp;
	}else{
		res=a*tmp;
	}
	b_cached = !b_cached;
	return res/sigmaI;
}
// nota la costellazione e il segnale in ingresso, secondo il criterio ML ricaviamo il simbolo trasmesso
std::complex<float> ML_receiver(int idx){
	/*float angle_rec = (float)atan2(imag(out[idx]),real(out[idx]));
	if (angle_rec <= 0) angle_rec += 2*PI;
	return angle_rec*180/PI;
	*/
	float dSquareNew,dSquareOld;
	int idxSymTX=0;
	
	dSquareOld = (float)pow(real<float>(out_noise[idx])-psk_I[0],(int)2) +(float)pow(imag<float>(out_noise[idx])-psk_Q[0],(int)2);
	for(int h=1; h<M_idx; h++){
		dSquareNew = (float)pow(real<float>(out_noise[idx])-psk_I[h],(int)2) +(float)pow(imag<float>(out_noise[idx])-psk_Q[h],(int)2);
		if(dSquareNew<dSquareOld){ // vuol dire che l'h-esimo è più vicino
			idxSymTX=h;	dSquareOld = dSquareNew;
		}
	}
	complex<float> R (K*(float)psk_I[idxSymTX],K*(float)psk_Q[idxSymTX]);
	return R;
}

int _tmain(int argc, _TCHAR* argv[])
{
	bool stop_ok = true;
	do{
		do{
			//procedura generazione simboli complessi semplificata
			for(int k=0; k<SYM_LENGTH; k++){ // kernel A, proverò a lanciarlo con <<1,SYM_LENGTH>> e via via aumentare la capacità dei pacchetti
				// generazione coppia di bit casuale e conversione in base 10
				int sym_tmp = ((int)(2*random()))*2+((int)(2*random())); //la funzione random sarà una __device__
				//int sym_map = mapping[sym_tmp]; // l'array di mapping sarà nella constant memory
				out[k] = std::complex<float>(K*psk_I[sym_tmp], K*psk_Q[sym_tmp]);
				//out[k] = std::complex<float>(K*(float)cos(sym_map*PI/M_idx), K*(float)sin(sym_map*PI/M_idx)); // ritornerà un puntatore all'array
				//cout<<"symbol: "<<out[k]<<endl;
 			}

			//parte eventuale di aggiunta rumore AWGN
			for(int i=0; i<SYM_LENGTH;i++){
				float a_1 = (float)EbNo_dB/10;
				float a = (float)sqrt(1.0/(2.0*2.0*pow(10,a_1)));//pow(10.0,(EbNo_dB/10))));
				float g1 = (float)sqrt(2.0)*(float)randn_notrig(0.0,a);
				float g2 = (float)sqrt(2.0)*(float)randn_notrig(0.0,a);
				out_noise[i] = std::complex<float>(real(out[i])+g1, imag(out[i])+g2);
				//cout<<"sym+noise compl: "<<out_noise[i]<<endl;
			}

			//uscita definitiva g(t)*(coseno*cos(w0*t)+seno*sin(w0*t))
			// non abbiamo portante, cosideriamo la banda base
			//immagino la trasmissione di simboli complessi
			for(int i=0; i<SYM_LENGTH; i++){
				ricevuti[i] = ML_receiver(i);
				//cout<<i<<" symbol ric: "<<ricevuti[i]<<endl;
			}
			//error check
			int err=0;
			for(int i=0; i<SYM_LENGTH; i++){
				if((real<float>(ricevuti[i])!=real<float>(out[i]))||(imag<float>(ricevuti[i])!=imag<float>(out[i]))){
					err++;
				}
			}
			// calcolo della SER non la BER!!!
			double a =(((double)err)/(2*(double)SYM_LENGTH));
			cout<<"SNR	"<<EbNo_dB<<"			"<<err<<"/"<<SYM_LENGTH<<"		"<<a<<endl;
			EbNo_dB +=1;
	
		}while(EbNo_dB<11);
		EbNo_dB=-5;
	system("pause");
	}while(1);
	return(0);
}