
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define PATH "C:/Users/jungu/Downloads"

__global__ void make_rx_ry(float* rx, float* ry, const float* xx, const float* yy, const float theta, const int nx, const int ny, const float M_PI)
{
    int i = threadIdx.x;
    int j = blockIdx.x;
    rx[j * nx + i] = xx[i * ny + j] * cos(theta - M_PI / 2.0) + yy[i * ny + j] * sin(theta - M_PI / 2.0);
    ry[j * nx + i] = -xx[i * ny + j] * sin(theta - M_PI / 2.0) + yy[i * ny + j] * cos(theta - M_PI / 2.0);

    // change xx and yy to 1-d array first, and after this function, change rx and ry into 1-d array, threadIdx.x and blockIdx.x must be 750
}

__global__ void make_pu_ratio(float* pu, float* Ratio, const float* rx, const float* ry, const float DSD, const float DSO, const float us, const float du, const int nx, const int ny)
{
    int i = threadIdx.x;
    int j = blockIdx.x;
    pu[j * ny + i] = (((rx[j * nx] * DSD / (ry[j * ny + i] + DSO))) + us) / (-du) + 1;
    Ratio[j * ny + i] = (DSO * DSO) / ((DSO + ry[j * ny + i]) * (DSO + ry[j * ny + i]));
}

__global__ void make_pv(float* pv, const float* ry, const float* zs, const float DSD, const float DSO, const float dv, const int ny, const float vs) {

    int i = threadIdx.x;
    int j = blockIdx.x;
    int k = blockIdx.y;
    pv[k * ny * ny + i * ny + j] = ((zs[k] * DSD) / (ry[i * ny + j] + DSO) - vs) / dv + 1;

}

__global__ void interp2(float* proj, float* pu, float* pv, float* result, const int proj_idx) {
    int i = threadIdx.x;
    int j = blockIdx.x;
    int k = blockIdx.y;
    //int *x, *y;

    //x[i * 750 + j] = (0 < pu[i * 750 + j] < 1628) ? ceil(pu[i * 750 + j]) : 1;
    //y[i * 750 + j] = (0 < pv[i * 750 + j] < 1500) ? ceil(pv[i * 750 + j]) : 1;
    result[k * 750 * 750 + i * 750 + j] = proj[1500 * 1628 * proj_idx + ((int)ceil(pu[i * 750 + j]) - 1) * 1628 + (int)ceil(pv[k * 750 * 750 + i * 750 + j]) - 1];
}

__global__ void image_final(float* img, const float* Ratio, const float* result, const int nx, const int ny) {

    int j = threadIdx.x;
    int k = blockIdx.x;
    int i = blockIdx.y;
    img[i * nx * ny + nx * j + k] += Ratio[j * ny + k] * result[i * nx * ny + j * ny + k];

}



void back_projection(float* final_img, float* proj,
    int nx, int ny, int nz, float sx, float sy, float sz, int nu, int nv, float su, float sv, float DSD, float DSO,
    float off_z, float off_u, float off_v, int num_angles, float* xs, float* ys, float* zs, float* us, float* vs) {

    //clock_t start1 = clock();
    float dx = sx / nx;
    float dy = sy / ny;
    float dz = sz / nz;
    float du = 0.098;
    float dv = 0.098;
    float M_PI = 3.141592;
    int i, j;

    float* xx = (float*)malloc(nx * ny * sizeof(float));
    float* yy = (float*)malloc(nx * ny * sizeof(float));

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            xx[i * ny + j] = xs[i];
            yy[i * ny + j] = ys[j];

        }
    }
    //clock_t end1 = clock();
    //printf("pre setting: %lf", (double)(end1 - start1));
    //printf("%f%f", xx[10000], yy[10000]);
    float* d_img, * d_zs, * d_xx, * d_yy, * d_proj;
    cudaMalloc(&d_img, nx * ny * nz * sizeof(float));
    cudaMalloc(&d_zs, nz * sizeof(float));
    cudaMalloc(&d_xx, nx * ny * sizeof(float));
    cudaMalloc(&d_yy, nx * ny * sizeof(float));
    cudaMalloc(&d_proj, nu * nv * num_angles * sizeof(float));
    cudaMemcpy(d_xx, xx, nx * ny * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_yy, yy, nx * ny * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_zs, zs, nz * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_proj, proj, nu * nv * num_angles * sizeof(float), cudaMemcpyHostToDevice);

    for (int proj_idx = 0; proj_idx < num_angles; proj_idx++) {

        float* d_rx, * d_ry, * d_pu, * d_ratio;
        //clock_t start2 = clock();
        cudaMalloc(&d_rx, nx * ny * sizeof(float));
        cudaMalloc(&d_ry, nx * ny * sizeof(float));
        cudaMalloc(&d_pu, nx * ny * sizeof(float));
        cudaMalloc(&d_ratio, nx * ny * sizeof(float));


        float theta = 2 * M_PI * proj_idx / num_angles;
        //clock_t end2 = clock();
        //printf("Memory copying: %lf", (double)(end2 - start2));
        //clock_t start3 = clock();
        make_rx_ry << < 750, 750 >> > (d_rx, d_ry, d_xx, d_yy, theta, nx, ny, M_PI);
        //clock_t end3 = clock();
        //float* rx = (float*)calloc(nu * nv, sizeof(float));
        //float* ry = (float*)calloc(nu * nv, sizeof(float));
        //cudaMemcpy(rx, d_rx, nx * ny * sizeof(float), cudaMemcpyDeviceToHost);
        //cudaMemcpy(ry, d_ry, nx * ny * sizeof(float), cudaMemcpyDeviceToHost);
        //printf("%f%f", rx[10000], ry[10000]);
        //printf("make rx ry: %lf", (double)(end3 - start3));
        //clock_t start4 = clock();
        make_pu_ratio << < 750, 750 >> > (d_pu, d_ratio, d_rx, d_ry, DSD, DSO, us[0], du, nx, ny);
        //clock_t end4 = clock();
        //printf("make pu ratio: %lf", (double)(end4 - start4));
        //float* pu = (float*)calloc(nx * ny, sizeof(float));
        //float* ratio = (float*)calloc(nx * ny, sizeof(float));
        //cudaMemcpy(pu, d_pu, nx * ny * sizeof(float), cudaMemcpyDeviceToHost);
        //cudaMemcpy(ratio, d_ratio, nx * ny * sizeof(float), cudaMemcpyDeviceToHost);
        //printf("%f%f", pu[11000], ratio[11000]);

        //change iz value


        float* d_pv, * d_result;

        cudaMalloc(&d_pv, nx * ny * nz * sizeof(float));
        cudaMalloc(&d_result, nx * ny * nz * sizeof(float));
        //clock_t start5 = clock();
        dim3 blocks1(750, 450);
        make_pv << < blocks1, 750 >> > (d_pv, d_ry, d_zs, DSD, DSO, dv, ny, vs[0]);
        //clock_t end5 = clock();
        //printf("make pv: %lf", (double)(end5 - start5));
            //float* pv = (float*)calloc(nx * ny, sizeof(float));
            //cudaMemcpy(pv, d_pv, nx * ny * sizeof(float), cudaMemcpyDeviceToHost);
            //printf("%f", pv[11000]);
        //clock_t start6 = clock();
        dim3 blocks2(750, 450);
        interp2 << < blocks2, 750 >> > (d_proj, d_pu, d_pv, d_result, proj_idx);
        //float* result = (float*)calloc(nx * ny, sizeof(float));
        //cudaMemcpy(result, d_result, nx * ny * sizeof(float), cudaMemcpyDeviceToHost);
        //printf("%f", result[300000]);
    //clock_t end6 = clock();
    //printf("interp2: %lf", (double)(end6 - start6));
    // start7 = clock();
        dim3 blocks3(750, 450);
        image_final << < blocks3, 750 >> > (d_img, d_ratio, d_result, nx, ny);
        //float* img = (float*)calloc(nx * ny * nz, sizeof(float));
        //cudaMemcpy(img, d_img, nx * ny * nz * sizeof(float), cudaMemcpyDeviceToHost);
        //printf("%f", img[112800000]);
    //clock_t end7 = clock();
    //printf("image final: %lf", (double)(end7 - start7));
        cudaFree(d_pv);
        cudaFree(d_result);
        //printf("%f", final_img[112800000]);
        cudaFree(d_rx);
        cudaFree(d_ry);
        cudaFree(d_pu);
        cudaFree(d_ratio);
        printf("%d", proj_idx);
    }
    cudaMemcpy(final_img, d_img, nx * ny * nz * sizeof(float), cudaMemcpyDeviceToHost);
    //printf("%f", final_img[112800000]);
    cudaFree(d_xx);
    cudaFree(d_yy);
    cudaFree(d_img);
    cudaFree(d_zs);
}

int main(int argc, char* argv[]) {
    int num_projections = 353; // number of projections
    int index_projection = 1;     // starting point of projections

    int nx = 750; // width voxels of image
    int ny = 750; // height voxels of image
    int nz = 450; // number of images

    float sx = 150; // real width of image [mm]
    float sy = 150; // real height of image [mm]
    float sz = 90;  // real width of image [mm]

    int nu = 1628; // width of projection
    int nv = 1500; // height of projection

    float su = 147;     // real width of projection [mm]
    float sv = 159.544; // real height of projection [mm]

    float DSD = 658.45;
    float DSO = 409.70;

    float off_z = -40;
    float off_u = -41.633;
    float off_v = -74.662;
    /**array param**/

    float dx = sx / nx;
    float dy = sy / ny;
    float dz = sz / nz;
    float du = su / nu;
    float dv = sv / nv;

    int i;
    float* xs = (float*)malloc(nx * sizeof(float));
    float* ys = (float*)malloc(ny * sizeof(float));
    float* zs = (float*)malloc(nz * sizeof(float));
    float* us = (float*)malloc(nu * sizeof(float));
    float* vs = (float*)malloc(nv * sizeof(float));

    // Generate xs array
    for (i = 0; i < nx; i++) {
        xs[i] = ((-(nx - 1) / 2.0) + i) * dx;
    }

    // Generate ys array
    for (i = 0; i < ny; i++) {
        ys[i] = ((-(ny - 1) / 2.0) + i) * dy;
    }

    // Generate zs array
    for (i = 0; i < nz; i++) {
        zs[i] = ((-(nz - 1) / 2.0) + i) * dz + off_z;
    }

    // Generate us array
    for (i = 0; i < nu; i++) {
        us[i] = ((-(nu - 1) / 2.0) + i) * du + off_u;
    }

    // Generate vs array
    for (i = 0; i < nv; i++) {
        vs[i] = ((-(nv - 1) / 2.0) + i) * dv + off_v;
    }

    printf("1, finished");
    float* projection =
        (float*)calloc(num_projections * nu * nv,
            sizeof(float)); // allocate memory to projections
    float* image =
        (float*)calloc(nx * ny * nz, sizeof(float)); // allocate memory to images

    for (int i = 0; i < num_projections; i++) {
        char filename[300];
        sprintf(filename, "%s/my_input/input_%03d.raw", PATH, i + index_projection);
        FILE* fp = fopen(filename, "rb"); // load file
        if (fp == NULL) {
            fputs("File error\n", stderr);
            exit(1);
        } // check if file is loaded

        float* buffer =
            (float*)calloc(nu * nv, sizeof(float)); // allocate memory to buffer
        fread(buffer, sizeof(float), nu * nv, fp);   // read file to buffer
        fclose(fp);
        for (int j = 0; j < nu * nv; j++) {
            //            buffer[j] = buffer[j] / 28415;
            //            if (buffer[j] == 0) buffer[j] = 1;
            //            buffer[j] = -log(buffer[j]);
            projection[nu * nv * i + j] =
                buffer[j]; // copy value of buffer to projections
        }
        free(buffer);
        if (i % 10 == 0)
            printf("reading %d\n", i + index_projection);
    }
    // open input files and store to memory
    printf("done reading the file");
    //    printf("%u\n", projection[32546]);
    //    printf("Projections Loaded.\n");
    //
    back_projection(image, projection, nx, ny, nz, sx, sy, sz, nu, nv, su, sv,
        DSD, DSO, off_z, off_u, off_v, num_projections, xs, ys, zs,
        us, vs);
    printf("Backprojection Done.\n");

    for (int i = 0; i < 450; i++) {
        char filename[100];
        sprintf(filename, "%s/output/cc_%04d.raw", PATH, i);
        FILE* fp = fopen(filename, "wb"); // open file to write
        if (fp == NULL) {
            fputs("File error\n", stderr);
            exit(1);
        } // check if file is loaded
        float* buffer =
            (float*)calloc(nx * ny, sizeof(float)); // allocate memory to buffer
        for (int j = 0; j < nx * ny; j++) {
            buffer[j] = image[nx * ny * i + j]; // copy value of images to buffer
        }
        fwrite(buffer, sizeof(float), nx * ny, fp); // read file to buffer
        fclose(fp);
        free(buffer);
    } // save output images to files

    printf("Images Saved.\n");

    free(projection);
    free(image);

    return 0;
}