#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(){
	int num_alpha=5;
	int nr = 4;
	int ntrials = 6;

	double trials_2pcf[num_alpha][ntrials][nr];
	double alpha_2pcf[num_alpha][nr];
	double alpha_err[num_alpha][nr];
	double chi2;
	double nbody_2pcf[nr];

	int ir,j,k;

	for (ir=0;ir<nr;ir++){
		nbody_2pcf[ir] = (double) 4.0;
		//fprintf(stderr,"%e, ",nbody_2pcf[ir]);
		for(j=0;j<num_alpha;j++){
			for(k=0;k<ntrials;k++){
				if (k%2==0)
				trials_2pcf[j][k][ir] = (double) 4.0 + 0.1*(j-2) + 0.05*k;
				else
					trials_2pcf[j][k][ir] = (double) 4.0 + 0.1*(j-2) - 0.05*k;
				//fprintf(stderr,"\n%e, ",trials_2pcf[j][k][ir]);
			}
		}
	}

	for(j=0;j<num_alpha;j++){
		//Get mean and stdev of trials_2pcf
		chi2 = 0.0;
		for (ir=0;ir<nr;ir++){
			alpha_2pcf[j][ir] = 0.0;
			for(k=0;k<ntrials;k++){
				alpha_2pcf[j][ir] += trials_2pcf[j][k][ir];//+total_nr-nr
			}
			alpha_2pcf[j][ir] = alpha_2pcf[j][ir]/ntrials;
			alpha_err[j][ir] = 0.0;
			for(k=0;k<ntrials;k++){
				alpha_err[j][ir] += pow((trials_2pcf[j][k][ir]-alpha_2pcf[j][ir]),2);//+total_nr-nr
			}
			alpha_err[j][ir] = pow((alpha_err[j][ir]/(ntrials)),0.5);

			// Now get chi^2 values
			chi2 += pow(((alpha_2pcf[j][ir]-nbody_2pcf[ir])/alpha_err[j][ir]),2);

		}

		fprintf(stderr,"%d, %lf\n",j,chi2);

		//	gsl_vector_set(chi2_alpha,j,chi2/nr);
	//	gsl_vector_set(weights,j,nr/chi2);
	}

	fprintf(stderr,"MEAN:\n");
	for(j=0;j<num_alpha;j++){
		for(ir=0;ir<nr;ir++){
			fprintf(stderr,"%e, ",alpha_2pcf[j][ir]);
		}
		fprintf(stderr,"\n");
	}

	fprintf(stderr,"STD:\n");
		for(j=0;j<num_alpha;j++){
			for(ir=0;ir<nr;ir++){
				fprintf(stderr,"%e, ",alpha_err[j][ir]);
			}
			fprintf(stderr,"\n");
		}

}
