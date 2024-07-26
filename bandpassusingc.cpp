#include <math.h>
#include <stdio.h>

#define PI 3.14159265358979323846

void bandpass_filter(double As, double Fs, double fs1, double fp1, double fp2, double fs2, int Kaiser, float *FIR_bandpass,int *N);

void convolution(float* FIR_bandpass,int N, float* input,int m,float* output){
	int length=m+N-1;
	int i,j;
	for(i=0; i<length;i++){
		output[i]=0;
	}
	for(i=0; i<N; i++){
		for(j=0; j<m; j++){
			output[i+j]+=FIR_bandpass[i]*input[j];
		}
	}
}

int main() {
    // Example usage
    int N,i;
    double As = 60.0;
    double Fs = 48000;
    double fs1 =7200 ;
    double fp1 = 9600;
    double fp2 = 12000;
    double fs2 = 14400;
    int Kaiser = 1;
	float input[100];
	float FIR_bandpass[1000];
	// 1 for Kaiser window, 0 for other windows
    float output[2000];
    // Calculate filter coefficients
    for(i=0;i<100; i++){
    	input[i]=sin(2*PI*10000*i/48000)+0.1*(sin(2*PI*800*i/48000))+0.5*(sin(2*PI*1500*i)/48000);
	}
    bandpass_filter(As, Fs, fs1, fp1, fp2, fs2, Kaiser, FIR_bandpass,&N);
    convolution(FIR_bandpass,N,input,100,output);
    
    printf("%i\n",N);
    
    for(i=0; i<100+N-1; i++){
    	printf("%f,",output[i]);
	}
	
	
    
    

    // Output the filter coefficients
    printf("FIR Bandpass Filter Coefficients:\n");
    for (int i = 0; i < N; ++i) {
        printf("%f\n", FIR_bandpass[i]);
    }

    return 0;
}

void bandpass_filter(double As, double Fs, double fs1, double fp1, double fp2, double fs2, int Kaiser, float *FIR_bandpass,int *N) {
    float fc1 = (fs1
	 + fp1) / 2;
    float fc2 = (fs2 + fp2) / 2;

    float Tb = 2 * PI * fmin((fs2 - fp2), (fp1 - fs1)) / Fs;

    float beta;
    if (Kaiser) {
        if (As > 50) {
            beta = 0.1102 * (As - 8.7);
        } else if (21 < As && As <= 50) {
            beta = 0.5842 * pow((As - 21), 0.4) + 0.07886 * (As - 21);
        } else {
            beta = 0;
        }
    }
      
    
    if (!Kaiser) {
        if (As <= 21) {
           *N = (int)ceil(1.8 * PI / Tb);
        } else if (As > 21 && As <= 26) {
            *N = (int)ceil(6.1 * PI / Tb);
        } else if (As > 26 && As <= 44) {
            *N = (int)ceil(6.2 * PI / Tb);
        } else if (As > 44 && As <= 53) {
            *N = (int)ceil(6.6 * PI / Tb);
        } else {
            *N = (int)ceil(11 * PI / Tb);
        }
    } else {
        *N= (int)ceil((As - 
		8) / (2.285 * Tb));
    }
    
    
    

    double w[*N];
    if (!Kaiser) {
        if (As <= 21) {
            for (int i = 0; i < *N; ++i) {
                w[i] = 1.0;
            }
        } else if (As > 21 && As <= 26) {
            for (int i = 0; i < *N; ++i) {
                w[i] = 1.0 - fabs(2.0 * i / (*N - 1) - 1.0);
            }
        } else if (As > 26 && As <= 44) {
            for (int i = 0; i < *N; ++i) {
                w[i] = 0.5 * (1.0 - cos(2.0 * PI * i / (*N - 1)));
            }
        } else if (As > 44 && As <= 53) {
            for (int i = 0; i < *N; ++i) {
                w[i] = 0.54 - 0.46 * cos(2.0 * PI * i / (*N - 1));
            }
        } else {
            for (int i = 0; i < *N; ++i) {
                w[i] = 0.42 - 0.5 * cos(2.0 * PI * i / (*N - 1)) + 0.08 * cos(4.0 * PI * i / (*N - 1));
            }
        }
    } else {
        for (int i = 0; i < *N; ++i) {
            w[i] = exp(beta * sqrt(1 - pow((2.0 * i - *N + 1) / (*N - 1), 2))) /
                   exp(beta);
        }
    }

    double alpha = *N / 2.0;
    for (int n = 0; n < *N; ++n) {
        if (n == alpha) {
            FIR_bandpass[n] = 2.0 * (fc2 / Fs - fc1 / Fs);
        } else {
            FIR_bandpass[n] = (sin(2.0 * PI * fc2 / Fs * (n - alpha)) -
                               sin(2.0 * PI * fc1 / Fs * (n - alpha))) /
                              (PI * (n - alpha));
        }
    }
   //  for (int i = 0; i < *N; ++i) {
    //    printf("%f\n",FIR_bandpass[i]);
  //  }

    for (int i = 0; i < *N; ++i) {
        FIR_bandpass[i] *= w[i];
    }
    //printf("%i\n",*N);
    //for(int i=0; i<*N; ++i){
    //	printf("%f\n",FIR_bandpass[i]);
	//}
}
