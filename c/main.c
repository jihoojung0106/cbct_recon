#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>

#define PATH "."

float** interp2(float** trans_proj, float **pu, float **pv, int nx, int ny, int nu, int nv) {
    float **result = (float **)calloc(nx, sizeof(float *));
    for (int i = 0; i < ny; i++) {
        result[i] = (float *)calloc(ny, sizeof(float));
    }
    
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            float u = pu[i][j];
            float v = pv[i][j];
            if (u >= 1 && u <= nu && v >= 1 && v <= nv) {
                int iu = floor(u + 0.5) - 1;
                int iv = floor(v + 0.5) - 1;
                result[i][j] = trans_proj[iu][iv];
            } else {
                result[i][j] = 0.0;
            }
        }
    }

    return result;
}

void back_projection(float*** img, float*** proj, 
int nx, int ny, int nz, float sx, float sy, float sz, int nu, int nv, float su, float sv, float DSD, float DSO, 
float off_z, float off_u, float off_v, int num_angles) {
    
    float dx = sx / nx;
    float dy = sy / ny;
    float dz = sz / nz;
    float du = su / nu;
    float dv = sv / nv;

    float *xs = (float *)calloc(nx, sizeof(float));
    for (int i = 0; i < nx; i++) {
        xs[i] = (- (nx - 1) / 2.0 + i) * dx;
    }

    float *ys = (float *)calloc(ny, sizeof(float));
    for (int i = 0; i < ny; i++) {
        ys[i] = (- (ny - 1) / 2.0 + i) * dy;
    }

    float *zs = (float *)calloc(nz, sizeof(float));
    for (int i = 0; i < nz; i++) {
        zs[i] = (- (nz - 1) / 2.0 + i) * dz + off_z;
    }

    float *us = (float *)calloc(nu, sizeof(float));
    for (int i = 0; i < nu; i++) {
        us[i] = (- (nu - 1) / 2.0 + i) * du + off_u;
    }

    float *vs = (float *)calloc(nv, sizeof(float));
    for (int i = 0; i < nx; i++) {
        vs[i] = (- (nv - 1) / 2.0 + i) * dv + off_v;
    }

    float **xx = (float **)calloc(nx, sizeof(float *));
    for (int i = 0; i < nx; i++) {
        xx[i] = (float *)calloc(ny, sizeof(float));
    }
    
    float **yy = (float **)calloc(nx, sizeof(float *));
    for (int i = 0; i < nx; i++) {
        yy[i] = (float *)calloc(ny, sizeof(float));
    }

    float **rx = (float **)calloc(nx, sizeof(float *));
    for (int i = 0; i < nx; i++) {
        rx[i] = (float *)calloc(ny, sizeof(float));
    }

    float **ry = (float **)calloc(nx, sizeof(float *));
    for (int i = 0; i < nx; i++) {
        ry[i] = (float *)calloc(ny, sizeof(float));
    }

    float **pu = (float **)calloc(nx, sizeof(float *));
    for (int i = 0; i < nx; i++) {
        pu[i] = (float *)calloc(ny, sizeof(float));
    }

    float **pv = (float **)calloc(nx, sizeof(float *));
    for (int i = 0; i < nx; i++) {
        pv[i] = (float *)calloc(ny, sizeof(float));
    }

    float **ratio = (float **)calloc(nx, sizeof(float *));
    for (int i = 0; i < nx; i++) {
        ratio[i] = (float *)calloc(ny, sizeof(float));
    }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            xx[i][j] = xs[i];
            yy[i][j] = ys[j];
        }
    }

    for (int proj_idx = 0; proj_idx < num_angles; proj_idx++) {
        printf("proj_idx = %d\n",proj_idx);
        float **trans_proj = (float **)calloc(nu, sizeof(float *));
        for (int i = 0; i < nu; i++) {
            trans_proj[i] = (float *)calloc(nv, sizeof(float));
            for (int j = 0; j < nv; j++) {
                trans_proj[i][j] = proj[proj_idx][j][i];
            }
        }
        float theta = 2 * M_PI * proj_idx / num_angles;
        float sin_theta = sin(theta - M_PI / 2.0);
        float cos_theta = cos(theta - M_PI / 2.0);
        
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                rx[i][j] = xx[i][j] * cos_theta + yy[i][j] * sin_theta;
                ry[i][j] = - xx[i][j] * sin_theta + yy[i][j] * cos_theta;
            }
        }
        
        for (int j = 0; j < nx; j++) {
            for (int i = 0; i < ny; i++) {
                pu[j][i] = (rx[j][i] * DSD / (ry[j][i] + DSO) + us[0]) / (- du) + 1;
                ratio[i][j] = pow(DSO, 2.0) / pow(DSO + ry[j][i], 2.0);
            }
        }
        
        for (int slice = 0; slice < nz ; slice++){
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    pv[j][i] = (zs[slice] * DSD / (ry[i][j] + DSO) - vs[0]) / dv + 1;
                }
            }
            
            float **result = interp2(trans_proj, pu, pv, nx, ny, nu, nv);
            
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    img[slice][j][i] = img[slice][j][i] + ratio[i][j] * result[i][j];
                }
            }
            
            for (int i = 0; i < nx; i++) {
                free(result[i]);
            }
            free(result);
        }

        for (int i = 0; i < nu; i++) {
            free(trans_proj[i]);
        }
        free(trans_proj);
    }

    free(xs);
    free(ys);
    free(zs);
    free(us);
    free(vs);

    for (int i = 0; i < nx; i++) {
        free(xx[i]);
    }
    free(xx);
    for (int i = 0; i < nx; i++) {
        free(yy[i]);
    }
    free(yy);
    for (int i = 0; i < nx; i++) {
        free(rx[i]);
    }
    free(rx);
    for (int i = 0; i < nx; i++) {
        free(ry[i]);
    }
    free(ry);
    for (int i = 0; i < nx; i++) {
        free(pu[i]);
    }
    free(pu);
    for (int i = 0; i < nx; i++) {
        free(pv[i]);
    }
    free(pv);
    for (int i = 0; i < nx; i++) {
        free(ratio[i]);
    }
    free(ratio);
}

int nextpow2(int n) {
    int next = pow(2, ceil(log(n)/log(2)));
    int power = log(next)/log(2);
    return power;
}

float *ramp_flat(int n) {
    float *nn = (float *)calloc(n, sizeof(float));
    for (int i = 0; i < n; i++) {
        nn[i] = - (n / 2.0) + i;
    }

    float *h = (float *)calloc(n, sizeof(float));
    h[n / 2] = 1 / 4.0;
    for (int i = 0; i < n; i++) {
        float remainder = fmod(nn[i], 2.0);
        if (remainder < 0) remainder = 2.0 + remainder;
        if (remainder == 1.0) {
            h[i] = - 1 / pow(M_PI * nn[i], 2.0);
        }
    }

    free(nn);
    return h;
} 

float *Filter(int filter, float *kernel, int order, float d) {
    fftw_complex *array_fft_ini = (fftw_complex *)fftw_malloc(order * sizeof(fftw_complex));
    fftw_complex *array_fft_fin = (fftw_complex *)fftw_malloc(order * sizeof(fftw_complex));
    
    fftw_plan fft = fftw_plan_dft_1d(order, array_fft_ini, array_fft_fin, FFTW_FORWARD, FFTW_MEASURE);
    for (int i = 0; i < order; i++) {
        array_fft_ini[i][0] = kernel[i];
    }
    
    fftw_execute(fft); 
    
    for (int i = 0; i < order; i++) {
        kernel[i] = 2 * fabs(array_fft_fin[i][0]);
    }

    float *f = (float *)calloc(order / 2 + 1, sizeof(float));
    for (int i = 0; i < order / 2 + 1; i++) {
        f[i] = kernel[i];
    }

    float *w = (float *)calloc(order / 2 + 1, sizeof(float));
    for (int i = 0; i < order / 2 + 1; i++) {
        w[i] = 2 * M_PI * i / order;
    }

    switch (filter) {
        case 0: // ram-lak
            break;
        case 1: // shepp-logan
            break;
        case 2: // cosine
            break;
        case 3: // hamming
            break;
        case 4: // hann
            break;
        default:
            break;
    }

    for (int i = 0; i < order / 2 + 1; i++) {
        if (w[i] > M_PI * d) {
            f[i] = 0;
        }
    }

    float *filt = (float *)calloc(order, sizeof(float));
    for (int i = 0; i < order; i++) {
        if (i < order / 2 + 1) {
            filt[i] = f[i];
        }
        else {
            filt[i] = f[order - i];
        }
    }

    fftw_destroy_plan(fft);
    fftw_free(array_fft_ini);
    fftw_free(array_fft_fin);
    free(f);
    free(w);
    return filt;
}

void filtering(float*** proj, int nu, int nv, float su, float sv, float DSD, float DSO, float off_u, float off_v, int num_angles) {
    
    float du = su / nu;
    float dv = sv / nv;

    float or = 31.9;

    float *us = (float *)calloc(nu, sizeof(float));
    for (int i = 0; i < nu; i++) {
        us[i] = (- (nu - 1) / 2.0 + i) * du + off_u;
    }

    float *vs = (float *)calloc(nv, sizeof(float));
    for (int i = 0; i < nv; i++) {
        vs[i] = (- (nv - 1) / 2.0 + i) * dv + off_v;
    }

    float **uu = (float **)calloc(nv, sizeof(float *));
    for (int i = 0; i < nv; i++) {
        uu[i] = (float *)calloc(nu, sizeof(float));
        for (int j = 0; j < nu; j++) {
            uu[i][j] = us[j];
        } 
    }

    float **weight = (float **)calloc(nv, sizeof(float *));
    for (int i = 0; i < nv; i++) {
        weight[i] = (float *)calloc(nu, sizeof(float));
        for (int j = 0; j < nu; j++) {
            if (fabs(uu[i][j]) > or) {
                weight[i][j] = 1;
            }
            else {
                weight[i][j] = pow(sin(M_PI * (uu[i][j] - or) / (or * 2.0 * 2.0)), 2.0);
            }
        }
    }

    for (int proj_idx = 0; proj_idx < num_angles; proj_idx++) {
        for (int i = 0; i < nv ; i++){
            for (int j = 0; j < nu; j++) {
                proj[proj_idx][i][j] = proj[proj_idx][i][j] * weight[i][j];
            }
        }
    }
    
    int filt_len = 64 > pow(2, nextpow2(2 * nu)) ? 64 : pow(2, nextpow2(2 * nu));
    float *ramp_kernel = ramp_flat(filt_len);
    float d = 0.1;
    float *filt = Filter(0, ramp_kernel, filt_len, d);
    float **filt_proj = (float **)calloc(nv, sizeof(float *));
    for (int i = 0; i < nv; i++) {
        filt_proj[i] = (float *)calloc(filt_len, sizeof(float));
    }
        
    for (int proj_idx = 0; proj_idx < num_angles; proj_idx++) {
        for (int i = 0; i < nv ; i++) {
            fftw_complex *array_fft_ini = (fftw_complex *)fftw_malloc(filt_len * sizeof(fftw_complex));
            fftw_complex *array_fft_mid = (fftw_complex *)fftw_malloc(filt_len * sizeof(fftw_complex));
            fftw_complex *array_fft_fin = (fftw_complex *)fftw_malloc(filt_len * sizeof(fftw_complex));

            fftw_plan fft = fftw_plan_dft_1d(filt_len, array_fft_ini, array_fft_mid, FFTW_FORWARD, FFTW_MEASURE);
            fftw_plan ifft = fftw_plan_dft_1d(filt_len, array_fft_mid, array_fft_fin, FFTW_BACKWARD, FFTW_MEASURE);

            for (int j = 0; j < nu; j++) {
                filt_proj[i][filt_len / 2 - nu / 2 + j] = proj[proj_idx][i][j];
            }
            
            for (int j = 0; j < filt_len; j++) {
                array_fft_ini[j][0] = filt_proj[i][j];
                array_fft_ini[j][1] = 0;
            }
            fftw_execute(fft);
            
            for (int j = 0; j < filt_len; j++) {
                array_fft_mid[j][0] = array_fft_mid[j][0] * filt[j];
                array_fft_mid[j][1] = array_fft_mid[j][1] * filt[j];
            }
            fftw_execute(ifft);
            
            for (int j = 0; j < filt_len; j++) {
                filt_proj[i][j] = array_fft_fin[j][0] / filt_len;
            }
            
            fftw_destroy_plan(fft);
            fftw_destroy_plan(ifft);
            fftw_free(array_fft_ini);
            fftw_free(array_fft_mid);
            fftw_free(array_fft_fin);            
        }
        for (int i = 0; i < nv ; i++) {
            for (int j = 0; j < nu; j++) {
                proj[proj_idx][i][j] = filt_proj[i][filt_len / 2 - nu / 2 + j] / 2 / du * (2 * M_PI / num_angles) / 2 * (DSD / DSO);
            }
        }
    }

    free(us);
    free(vs);
    for (int i = 0; i < nv; i++) {
        free(uu[i]);
    }
    free(uu);
    for (int i = 0; i < nv; i++) {
        free(weight[i]);
    }
    free(weight);
    free(ramp_kernel);
    free(filt);
    for (int i = 0; i < nv; i++) {
        free(filt_proj[i]);
    }
    free(filt_proj);
}

int main(int argc, char *argv[]) {
    int num_projections = 706; // number of projections
    int index_projection = 43; // starting point of projections

    int nx = 750; // width voxels of image
    int ny = 750; // height voxels of image
    int nz = 450; // number of images

    float sx = 150; // real width of image [mm]
    float sy = 150; // real height of image [mm]
    float sz = 90; // real width of image [mm]
    
    int nu = 1500; // width of projection
    int nv = 1628; // height of projection

    float su = 147; // real width of projection [mm]
    float sv = 159.544; // real height of projection [mm]

    float DSD = 658.45;
    float DSO = 409.70;

    float off_z = -40;
    float off_u = -41.633;
    float off_v = -74.662;

    float ***projection = (float ***)calloc(num_projections, sizeof(float **));
    for (int i = 0; i < num_projections; i++) {
        projection[i] = (float **)calloc(nv, sizeof(float *));
        for (int j = 0; j < nv; j++) {
            projection[i][j] = (float *)calloc(nu, sizeof(float));
        }
    } // allocate memory of 3D array for projections

    float ***image = (float ***)calloc(nz, sizeof(float **));
    for (int i = 0; i < nz; i++) {
        image[i] = (float **)calloc(ny, sizeof(float *));
        for (int j = 0; j < ny; j++) {
            image[i][j] = (float *)calloc(nx, sizeof(float));
        }
    } // allocate memory of 3D array for images

    for (int i = 0; i < num_projections; i++) {
        char filename[100];
        sprintf(filename, "%s/input/%04d.raw", PATH, i + index_projection);
        FILE *fp = fopen(filename, "rb"); // load file
        if (fp == NULL) {
            fputs("File error\n", stderr);
            exit(1);
        } // check if file is loaded
        unsigned short *buffer = (unsigned short *)calloc(nu * nv, sizeof(unsigned short)); // allocate memory to buffer
        fread(buffer, sizeof(unsigned short), nu * nv, fp); // read file to buffer
        for (int j = 0; j < nv; j++) {
            for (int k = 0; k < nu; k++) {
                projection[i][j][k] = buffer[nu * j + k] / 28415.0;                
                if (projection[i][j][k] == 0) projection[i][j][k] = 1;
                projection[i][j][k] = -log(projection[i][j][k]); // copy nomalized value of buffer to projections
            }         
        }
        fclose(fp);
        free(buffer);
    } // open input files and store to memory
    printf("Projections Loaded.\n");

    filtering(projection, nu, nv, su, sv, DSD, DSO, off_u, off_v, num_projections);
    printf("Filtering Done.\n");

    back_projection(image, projection, nx, ny, nz, sx, sy, sz, nu, nv, su, sv, DSD, DSO, off_z, off_u, off_v, num_projections);
    printf("Backprojection Done.\n");

    for (int i = 0; i < nz; i++) {
        char filename[100];
        sprintf(filename, "%s/output/%04d.raw", PATH, i);
        FILE *fp = fopen(filename, "wb"); // open file to write
        if (fp == NULL) {
            fputs("File error\n", stderr);
            exit(1);
        } // check if file is loaded
        short *buffer = (short *)calloc(nx * ny, sizeof(short)); // allocate memory to buffer
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nx; k++) {
                if (image[i][j][k] < 0.022802) {
                    image[i][j][k] = image[i][j][k] * 60775.4 + 1385.80;
                }
                else {
                    image[i][j][k] = image[i][j][k] * 68036.46 + 1551.36;
                }
                buffer[nx * j + k] = (short)image[i][j][k]; // copy value of images to buffer
            }
        }
        fwrite(buffer, sizeof(short), nx * ny, fp); // read file to buffer
        fclose(fp);
        free(buffer);
    } // save output images to files
    printf("Images Saved.\n");

    for (int i = 0; i < num_projections; i++) {
        for (int j = 0; j < nv; j++) {
            free(projection[i][j]);
        }
        free(projection[i]);
    }
    free(projection);

    for (int i = 0; i < nz; i++) {
        for (int j = 0; j < ny; j++) {
            free(image[i][j]);
        }
        free(image[i]);
    }
    free(image);

    return 0;
}