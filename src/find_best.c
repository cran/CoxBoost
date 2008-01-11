#include <math.h>
//#include <stdio.h>

void find_best(double *x, int *n, int *p, int *uncens, int *uncens_n, double *beta, double *risk, double *eta, 
               double *weights, double *weighted_risk, double *weighted_risk_sum, 
               double *penalty,
               int *min_index, double *min_deviance, double *min_beta_delta)
{
    int actual_x, i, j;
    double x_bar,V,ldenom;
    double U, I, beta_delta, deviance;
    
    double score, max_score, max_score_beta_delta;
    int max_score_index;
    
    for (actual_x = 0; actual_x < *p; actual_x++) {
        
        U = 0;
        I = 0;
        
        for (i = 0; i < *uncens_n; i++) {
            x_bar = 0;
            for (j = 0; j < *n; j++) x_bar += weighted_risk[(i * *n) + j] * x[(actual_x * *n) + j];
            x_bar /= weighted_risk_sum[i];
            U += x[(actual_x * *n) + uncens[i]] - x_bar; 
            
            V = 0;
            for (j = 0; j < *n; j++) V += weighted_risk[(i * *n) + j] * (x[(actual_x * *n) + j] - x_bar) * (x[(actual_x * *n) + j] - x_bar);
            V /= weighted_risk_sum[i]; 
            I += V;
        }
        
        beta_delta = U/(I + *penalty);
        
        score = U*U/(I + *penalty);

        if (actual_x == 0 || score > max_score) {
            max_score_index = actual_x + 1;
            max_score = score;
            max_score_beta_delta = beta_delta;
        }
        
        // deviance = 0;
        // 
        // for (i = 0; i < *uncens_n; i++) {
        //     ldenom = 0;
        //     for (j = 0; j < *n; j++) ldenom += weights[(i * *n) + j] * risk[j] * exp(x[(actual_x * *n) + j]*beta_delta);
        //     deviance += eta[uncens[i]] + x[(actual_x * *n) + uncens[i]]*beta_delta - log(ldenom); 
        // }
        // 
        // deviance = -2*deviance;
        //         
        // if (actual_x == 0 || deviance < *min_deviance) {
        //     *min_index = actual_x + 1;
        //     *min_deviance = deviance;
        //     *min_beta_delta = beta_delta;
        // }
    }

    *min_index = max_score_index;
    *min_beta_delta = max_score_beta_delta;
    
    *min_deviance = 0; 

    for (i = 0; i < *uncens_n; i++) {
        ldenom = 0;
        for (j = 0; j < *n; j++) ldenom += weights[(i * *n) + j] * risk[j] * exp(x[((*min_index - 1) * *n) + j]* *min_beta_delta);
        *min_deviance += eta[uncens[i]] + x[((*min_index - 1) * *n) + uncens[i]]* *min_beta_delta - log(ldenom); 
    }

    //printf("min %d %d\n",max_score_index,*min_index);
}
